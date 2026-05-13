#!/usr/bin/env python3
import argparse
import gzip
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Harmonize cohort VCF REF/ALT to match a reference VCF by CHROM:POS")
    parser.add_argument("--cohort-vcf", required=True, help="Input cohort VCF (.vcf.gz)")
    parser.add_argument("--reference-vcf", required=True, help="Input reference VCF (.vcf.gz)")
    parser.add_argument("--output-vcf", required=True, help="Output harmonized VCF (.vcf)")
    parser.add_argument("--report", required=True, help="Output text report path")
    return parser.parse_args()


def load_reference_map(reference_vcf: str) -> dict[tuple[str, str], tuple[str, str]]:
    cmd = [
        "bcftools",
        "query",
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT\n",
        reference_vcf,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"bcftools query failed on reference VCF: {proc.stderr.strip()}")

    mapping: dict[tuple[str, str], tuple[str, str]] = {}
    for raw in proc.stdout.splitlines():
        if not raw:
            continue
        chrom, pos, ref, alt = raw.split("\t")
        if "," in alt:
            continue
        mapping[(chrom.replace("chr", ""), pos)] = (ref.upper(), alt.upper())
    return mapping


def flip_gt(gt_field: str) -> str:
    if gt_field in {".", "./.", ".|."}:
        return gt_field

    sep = "/"
    if "|" in gt_field:
        sep = "|"

    alleles = gt_field.split(sep)
    flipped = []
    for allele in alleles:
        if allele == ".":
            flipped.append(".")
        elif allele == "0":
            flipped.append("1")
        elif allele == "1":
            flipped.append("0")
        else:
            # Unexpected multi-allelic/haplotype code, preserve to avoid corrupting data.
            flipped.append(allele)
    return sep.join(flipped)


def harmonize_sample(sample_col: str) -> str:
    if sample_col == ".":
        return sample_col
    parts = sample_col.split(":")
    if not parts:
        return sample_col
    parts[0] = flip_gt(parts[0])
    return ":".join(parts)


def main():
    args = parse_args()

    ref_map = load_reference_map(args.reference_vcf)
    if not ref_map:
        raise RuntimeError("No usable biallelic reference variants were loaded")

    output_path = Path(args.output_vcf)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    Path(args.report).parent.mkdir(parents=True, exist_ok=True)

    total = 0
    kept_same = 0
    kept_swapped = 0
    dropped_not_in_ref = 0
    dropped_mismatch = 0

    with gzip.open(args.cohort_vcf, "rt", encoding="utf-8", errors="ignore") as fin, output_path.open(
        "w", encoding="utf-8"
    ) as fout:
        for raw in fin:
            if raw.startswith("#"):
                fout.write(raw)
                continue

            total += 1
            line = raw.rstrip("\n")
            parts = line.split("\t")
            if len(parts) < 8:
                continue

            chrom = parts[0].replace("chr", "")
            pos = parts[1]
            ref = parts[3].upper()
            alt = parts[4].upper()

            if "," in alt:
                dropped_mismatch += 1
                continue

            ref_entry = ref_map.get((chrom, pos))
            if ref_entry is None:
                dropped_not_in_ref += 1
                continue

            ref_ref, ref_alt = ref_entry

            if ref == ref_ref and alt == ref_alt:
                fout.write(raw)
                kept_same += 1
                continue

            if ref == ref_alt and alt == ref_ref:
                parts[3] = ref_ref
                parts[4] = ref_alt
                if len(parts) > 9:
                    for i in range(9, len(parts)):
                        parts[i] = harmonize_sample(parts[i])
                fout.write("\t".join(parts) + "\n")
                kept_swapped += 1
                continue

            dropped_mismatch += 1

    with Path(args.report).open("w", encoding="utf-8") as rep:
        rep.write(f"total_variants={total}\n")
        rep.write(f"kept_same={kept_same}\n")
        rep.write(f"kept_swapped={kept_swapped}\n")
        rep.write(f"dropped_not_in_ref={dropped_not_in_ref}\n")
        rep.write(f"dropped_mismatch={dropped_mismatch}\n")

    print(
        "Harmonization complete: "
        f"same={kept_same}, swapped={kept_swapped}, "
        f"not_in_ref={dropped_not_in_ref}, mismatch={dropped_mismatch}"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        raise

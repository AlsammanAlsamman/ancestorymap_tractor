#!/usr/bin/env python3
import argparse
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path


GT_SPLIT_RE = re.compile(r"[\/|]")


def parse_args():
    parser = argparse.ArgumentParser(description="Calculate SNP MAF from a VCF using optional sample subsets")
    parser.add_argument("--vcf", required=True, help="Input VCF/BCF file")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--label", required=True, help="Source label to write in the output")
    parser.add_argument("--sample-file", default="", help="Optional file with one sample ID per line")
    parser.add_argument("--regions-file", default="", help="Optional loci file with columns including chr/start/end")
    return parser.parse_args()


def parse_gt(gt_value: str) -> tuple[int, int]:
    gt = gt_value.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return 0, 0

    alt_count = 0
    allele_number = 0
    for allele in GT_SPLIT_RE.split(gt):
        if allele == "." or allele == "":
            continue
        if allele not in {"0", "1"}:
            raise ValueError(f"Non-biallelic allele '{allele}' found in genotype '{gt_value}'")
        allele_number += 1
        if allele == "1":
            alt_count += 1
    return alt_count, allele_number


def build_regions_bed(regions_file: str) -> str:
    temp_handle = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False, encoding="utf-8")
    with open(regions_file, "r", encoding="utf-8") as handle:
        next(handle, None)
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            chrom = str(parts[1]).replace("chr", "")
            start = int(parts[2])
            end = int(parts[3])
            temp_handle.write(f"{chrom}\t{start}\t{end}\n")
    temp_handle.close()
    return temp_handle.name


def stream_variants(vcf_path: str, sample_file: str, regions_file: str):
    view_cmd = ["bcftools", "view", "-m2", "-M2", "-v", "snps"]
    temp_regions_bed = ""
    if regions_file:
        temp_regions_bed = build_regions_bed(regions_file)
        view_cmd.extend(["-R", temp_regions_bed])
    if sample_file:
        view_cmd.extend(["-S", sample_file])
    view_cmd.extend(["-Ou", vcf_path])

    query_cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n"]

    with subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as view_proc:
        with subprocess.Popen(
            query_cmd,
            stdin=view_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        ) as query_proc:
            assert view_proc.stdout is not None
            view_proc.stdout.close()

            assert query_proc.stdout is not None
            for line in query_proc.stdout:
                yield line.rstrip("\n")

            query_stderr = query_proc.stderr.read() if query_proc.stderr is not None else ""
            query_returncode = query_proc.wait()
            view_stderr = view_proc.stderr.read().decode() if view_proc.stderr is not None else ""
            view_returncode = view_proc.wait()

    if temp_regions_bed:
        Path(temp_regions_bed).unlink(missing_ok=True)

    if view_returncode != 0:
        raise RuntimeError(f"bcftools view failed for {vcf_path}:\n{view_stderr.strip()}")
    if query_returncode != 0:
        raise RuntimeError(f"bcftools query failed for {vcf_path}:\n{query_stderr.strip()}")


def main():
    args = parse_args()
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows_written = 0
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chr", "pos", "rsid", "ref", "alt", "maf", "maf_allele", "alt_af", "alt_ac", "an", "n_samples_called", "source"])

        for record in stream_variants(args.vcf, args.sample_file, args.regions_file):
            if not record:
                continue
            parts = record.split("\t")
            if len(parts) < 5:
                continue

            chrom, pos, rsid, ref, alt = parts[:5]
            genotypes = parts[5:]

            alt_ac = 0
            an = 0
            n_samples_called = 0
            for gt in genotypes:
                sample_alt, sample_an = parse_gt(gt)
                alt_ac += sample_alt
                an += sample_an
                if sample_an > 0:
                    n_samples_called += 1

            if an == 0:
                alt_af = ""
                maf = ""
                maf_allele = ""
            else:
                alt_af_value = alt_ac / an
                maf_value = min(alt_af_value, 1.0 - alt_af_value)
                maf_allele = alt if alt_af_value <= 0.5 else ref
                alt_af = f"{alt_af_value:.8f}"
                maf = f"{maf_value:.8f}"

            writer.writerow([chrom.replace("chr", ""), pos, rsid, ref, alt, maf, maf_allele, alt_af, alt_ac, an, n_samples_called, args.label])
            rows_written += 1

    print(f"Wrote {rows_written} SNP rows to {output_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        raise

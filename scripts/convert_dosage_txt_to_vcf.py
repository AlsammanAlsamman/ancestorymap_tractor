#!/usr/bin/env python3
import argparse
import gzip


def to_gt(value):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "./."
    if x < 0:
        return "./."
    dosage = int(round(x))
    if dosage <= 0:
        return "0/0"
    if dosage == 1:
        return "0/1"
    if dosage >= 2:
        return "1/1"
    return "./."


def main():
    parser = argparse.ArgumentParser(description="Convert ancestry dosage table to a simple GT VCF")
    parser.add_argument("--dosage", required=True, help="Input dosage table from extract_tracts.py")
    parser.add_argument("--output-vcf-gz", required=True, help="Output VCF.GZ path")
    parser.add_argument("--source", default="tractor_dosage", help="Source tag for VCF header")
    args = parser.parse_args()

    with open(args.dosage, "r", encoding="utf-8") as src:
        header = src.readline().rstrip("\n").split("\t")
        if len(header) < 6:
            raise ValueError(f"Unexpected dosage header format in {args.dosage}")
        sample_ids = header[5:]

        with gzip.open(args.output_vcf_gz, "wt", encoding="utf-8") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write(f"##source={args.source}\n")
            out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Converted genotype from ancestry dosage\">\n")
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
            out.write("\t".join(sample_ids) + "\n")

            for line in src:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue

                chrom = parts[0]
                pos = parts[1]
                rsid = parts[2] if parts[2] else "."
                ref = parts[3] if parts[3] else "N"
                alt = parts[4] if parts[4] else "<ALT>"

                gts = [to_gt(v) for v in parts[5:]]
                out.write(
                    f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n"
                )


if __name__ == "__main__":
    main()

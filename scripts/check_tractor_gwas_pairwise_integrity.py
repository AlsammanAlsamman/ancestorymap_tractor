#!/usr/bin/env python3
import argparse
import csv
import os
import sys


def read_first_row(path: str):
    with open(path, "r", encoding="utf-8") as handle:
      reader = csv.DictReader(handle, delimiter="\t")
      return reader.fieldnames, next(reader, None)


def main():
    parser = argparse.ArgumentParser(description="Validate Tractor pairwise GWAS outputs")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--pairs", nargs="+", required=True)
    parser.add_argument("--chromosomes", nargs="+", required=True)
    args = parser.parse_args()

    for pair in args.pairs:
        if "_" not in pair:
            raise ValueError(f"Invalid pair name: {pair}")
        label_i, label_j = pair.split("_", 1)

        merged_path = os.path.join(args.output_dir, f"merged.model_{pair}.gwas.tsv")
        png_path = os.path.join(args.output_dir, f"merged.model_{pair}.manhattan.png")
        pdf_path = os.path.join(args.output_dir, f"merged.model_{pair}.manhattan.pdf")

        for path in (merged_path, png_path, pdf_path):
            if not os.path.exists(path):
                raise FileNotFoundError(f"Missing expected pairwise output: {path}")

        fields, first_row = read_first_row(merged_path)
        required = {
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "p_dosage_i",
            "p_dosage_j",
            "ancestry_i",
            "ancestry_j",
        }
        missing = required.difference(fields or [])
        if missing:
            raise ValueError(f"Merged GWAS output missing columns for {pair}: {sorted(missing)}")
        if first_row is None:
            raise ValueError(f"Merged GWAS output has no rows: {merged_path}")
        if str(first_row["ancestry_i"]) != label_i or str(first_row["ancestry_j"]) != label_j:
            raise ValueError(
                f"Merged GWAS ancestry labels do not match pair {pair}: "
                f"found ancestry_i={first_row['ancestry_i']}, ancestry_j={first_row['ancestry_j']}"
            )

        for chrom in args.chromosomes:
            per_chr_path = os.path.join(args.output_dir, f"chr{chrom}.model_{pair}.gwas.tsv")
            if not os.path.exists(per_chr_path):
                raise FileNotFoundError(f"Missing per-chromosome GWAS output: {per_chr_path}")
            per_fields, per_first_row = read_first_row(per_chr_path)
            missing = required.difference(per_fields or [])
            if missing:
                raise ValueError(f"Per-chromosome GWAS output missing columns for {pair}, chr{chrom}: {sorted(missing)}")
            if per_first_row is not None:
                if str(per_first_row["ancestry_i"]) != label_i or str(per_first_row["ancestry_j"]) != label_j:
                    raise ValueError(
                        f"Per-chromosome GWAS ancestry labels do not match pair {pair}, chr{chrom}: "
                        f"found ancestry_i={per_first_row['ancestry_i']}, ancestry_j={per_first_row['ancestry_j']}"
                    )

    print("Tractor pairwise GWAS ancestry integrity check passed.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"tractor_gwas_pairwise integrity check failed: {exc}", file=sys.stderr)
        sys.exit(1)
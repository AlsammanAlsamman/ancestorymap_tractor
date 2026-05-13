#!/usr/bin/env python3
import argparse
import csv
import sys


def main():
    parser = argparse.ArgumentParser(description="Validate merged MAF + GWAS summary output")
    parser.add_argument("--summary", required=True)
    parser.add_argument("--ancestries", nargs="+", required=True)
    parser.add_argument("--pairs", nargs="+", required=True)
    args = parser.parse_args()

    with open(args.summary, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        first_row = next(reader, None)

    if first_row is None:
      raise ValueError(f"Merged summary is empty: {args.summary}")

    required_base = [
        "chr",
        "pos",
        "rsid",
        "locus_name",
        "gene",
        "gene_distance_bp",
        "ref",
        "alt",
        "cohort_maf",
    ]
    required = set(required_base)
    required.update(f"{ancestry}_maf" for ancestry in args.ancestries)
    required.update(f"OR_local_ancestry_{pair}" for pair in args.pairs)
    required.update(f"p_local_ancestry_{pair}" for pair in args.pairs)
    required.update(f"p_dosage_{ancestry}" for ancestry in args.ancestries)
    required.update(f"OR_dosage_{ancestry}" for ancestry in args.ancestries)
    required.update(f"se_dosage_{ancestry}" for ancestry in args.ancestries)

    missing = sorted(required.difference(fields))
    if missing:
        raise ValueError(f"Merged summary missing required columns: {missing}")

    for column in required_base:
        if str(first_row.get(column, "")).strip() == "":
            raise ValueError(f"Merged summary first row has empty required field: {column}")

    print("Merged MAF + GWAS summary ancestry integrity check passed.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"merge_maf_gwas_summary integrity check failed: {exc}", file=sys.stderr)
        sys.exit(1)
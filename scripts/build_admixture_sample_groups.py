#!/usr/bin/env python3
import argparse
import csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build sample group table for ADMIXTURE plotting (reference populations + cohort)."
    )
    parser.add_argument("--fam", required=True, help="Merged PLINK .fam file")
    parser.add_argument("--reference-sample-map", required=True, help="sample_id ancestry_label table")
    parser.add_argument("--cohort-label", default="Cohort", help="Label to assign non-reference samples")
    parser.add_argument("--output", required=True, help="Output TSV path")
    return parser.parse_args()


def load_reference_map(path: str) -> dict:
    mapping = {}
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            mapping[parts[0]] = parts[1]
    return mapping


def main() -> None:
    args = parse_args()
    ref_map = load_reference_map(args.reference_sample_map)

    with open(args.fam, "r", encoding="utf-8") as fam_handle, open(
        args.output, "w", encoding="utf-8", newline=""
    ) as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["IID", "GROUP", "SOURCE"])

        for raw in fam_handle:
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            iid = parts[1]
            if iid in ref_map:
                writer.writerow([iid, ref_map[iid], "reference"])
            else:
                writer.writerow([iid, args.cohort_label, "cohort"])


if __name__ == "__main__":
    main()

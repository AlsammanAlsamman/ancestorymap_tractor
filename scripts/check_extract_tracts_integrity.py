#!/usr/bin/env python3
import argparse
import os
import sys


def parse_label_code_map(raw: str):
    mapping = {}
    for pair in raw.split(","):
        item = pair.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"Invalid label/code entry: {item}")
        label, code = item.split(":", 1)
        label = label.strip()
        code = code.strip()
        if not label or not code:
            raise ValueError(f"Empty label/code entry: {item}")
        if not code.isdigit():
            raise ValueError(f"Non-integer code in mapping: {item}")
        code_int = int(code)
        if label in mapping:
            raise ValueError(f"Duplicate label in mapping: {label}")
        if code_int in mapping.values():
            raise ValueError(f"Duplicate code in mapping: {code_int}")
        mapping[label] = code_int
    if not mapping:
        raise ValueError("No ancestry mappings provided")
    return mapping


def read_header(path: str):
    with open(path, "r", encoding="utf-8") as handle:
        first = handle.readline().strip()
    if not first:
        raise ValueError(f"Empty output file: {path}")
    return first.split("\t")


def main():
    parser = argparse.ArgumentParser(description="Validate Tractor extract_tracts outputs")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--num-ancs", type=int, required=True)
    parser.add_argument("--label-code-map", required=True)
    parser.add_argument("--chromosomes", nargs="+", required=True)
    args = parser.parse_args()

    label_to_code = parse_label_code_map(args.label_code_map)
    expected_codes = sorted(label_to_code.values())
    expected_sequence = list(range(args.num_ancs))
    if expected_codes != expected_sequence:
        raise ValueError(
            f"Expected contiguous ancestry codes 0..{args.num_ancs - 1}, found {expected_codes}"
        )

    for chrom in args.chromosomes:
        dosage_headers = []
        hapcount_headers = []
        for anc in expected_sequence:
            dosage = os.path.join(args.output_dir, f"chr{chrom}.phased.anc{anc}.dosage.txt")
            hapcount = os.path.join(args.output_dir, f"chr{chrom}.phased.anc{anc}.hapcount.txt")
            if not os.path.exists(dosage):
                raise FileNotFoundError(f"Missing dosage output: {dosage}")
            if not os.path.exists(hapcount):
                raise FileNotFoundError(f"Missing hapcount output: {hapcount}")
            dosage_headers.append(read_header(dosage))
            hapcount_headers.append(read_header(hapcount))

        first_dosage = dosage_headers[0]
        first_hapcount = hapcount_headers[0]
        for header in dosage_headers[1:]:
            if header != first_dosage:
                raise ValueError(f"Dosage headers do not match for chr{chrom}")
        for header in hapcount_headers[1:]:
            if header != first_hapcount:
                raise ValueError(f"Hapcount headers do not match for chr{chrom}")
        if first_dosage != first_hapcount:
            raise ValueError(f"Dosage/hapcount headers do not match for chr{chrom}")

    print("Tractor extract_tracts ancestry integrity check passed.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"extract_tracts integrity check failed: {exc}", file=sys.stderr)
        sys.exit(1)

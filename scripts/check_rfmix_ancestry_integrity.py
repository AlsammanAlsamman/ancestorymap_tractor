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
            raise ValueError(f"Invalid --label-code-map entry (expected LABEL:CODE): {item}")
        label, code = item.split(":", 1)
        label = label.strip()
        code = code.strip()
        if not label or not code:
            raise ValueError(f"Invalid --label-code-map entry (empty label or code): {item}")
        if label in mapping:
            raise ValueError(f"Duplicate label in --label-code-map: {label}")
        if code in mapping.values():
            raise ValueError(f"Duplicate code in --label-code-map: {code}")
        if not code.isdigit():
            raise ValueError(f"Non-integer code in --label-code-map: {code}")
        mapping[label] = int(code)

    if not mapping:
        raise ValueError("--label-code-map produced no ancestry mappings")
    return mapping


def read_sample_map_labels(path: str):
    labels = set()
    with open(path, "r", encoding="utf-8") as src:
        for line in src:
            row = line.strip()
            if not row:
                continue
            parts = row.split()
            if len(parts) < 2:
                continue
            labels.add(parts[1])
    return labels


def scan_msp_codes(paths, expected_codes):
    seen_codes = set()
    for msp in paths:
        if not os.path.exists(msp):
            raise FileNotFoundError(f"Missing MSP file: {msp}")
        with open(msp, "r", encoding="utf-8") as src:
            for line in src:
                row = line.strip()
                if not row or row.startswith("#"):
                    continue
                parts = row.split()
                for token in parts[6:]:
                    if token.lstrip("-").isdigit():
                        code = int(token)
                        if code not in expected_codes:
                            raise ValueError(
                                f"Observed code {code} not in expected codes {sorted(expected_codes)} (file: {msp})"
                            )
                        seen_codes.add(code)
    return seen_codes


def main():
    parser = argparse.ArgumentParser(description="Validate ancestry label/code integrity for RFMix outputs")
    parser.add_argument("--sample-map", required=True)
    parser.add_argument("--label-code-map", required=True)
    parser.add_argument("--msp-files", nargs="+", required=True)
    args = parser.parse_args()

    label_to_code = parse_label_code_map(args.label_code_map)
    expected_labels = set(label_to_code.keys())
    expected_codes = set(label_to_code.values())

    observed_labels = read_sample_map_labels(args.sample_map)
    if observed_labels != expected_labels:
        missing = sorted(expected_labels - observed_labels)
        extra = sorted(observed_labels - expected_labels)
        raise ValueError(
            "Sample map labels mismatch. "
            f"Missing labels: {missing}; Unexpected labels: {extra}"
        )

    seen_codes = scan_msp_codes(args.msp_files, expected_codes)
    if seen_codes != expected_codes:
        missing_codes = sorted(expected_codes - seen_codes)
        raise ValueError(
            "Not all expected ancestry codes were observed in MSP outputs. "
            f"Missing codes: {missing_codes}; Seen codes: {sorted(seen_codes)}"
        )

    print("RFMix ancestry integrity check passed.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"RFMix ancestry check failed: {exc}", file=sys.stderr)
        sys.exit(1)

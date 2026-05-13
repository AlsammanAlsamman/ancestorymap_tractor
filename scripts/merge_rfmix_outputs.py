#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sys
from typing import Dict, List


CHR_RE = re.compile(r"chr(?P<chr>[^./]+)\.deconvoluted")
ANC_TOKEN_RE = re.compile(r"^(?:anc)?(?P<code>\d+)$", re.IGNORECASE)


def parse_label_code_map(raw: str):
    label_to_code = {}
    for pair in raw.split(","):
        item = pair.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"Invalid label-code map entry: {item}")
        label, code = item.split(":", 1)
        label = label.strip()
        code = code.strip()
        if not label or not code:
            raise ValueError(f"Invalid label-code map entry: {item}")
        if not code.isdigit():
            raise ValueError(f"Non-numeric ancestry code in label-code map: {item}")
        if label in label_to_code:
            raise ValueError(f"Duplicate ancestry label in label-code map: {label}")
        if code in label_to_code.values():
            raise ValueError(f"Duplicate ancestry code in label-code map: {code}")
        label_to_code[label] = code
    if not label_to_code:
        raise ValueError("Label-code map is empty")
    return label_to_code


def infer_chr(path: str) -> str:
    match = CHR_RE.search(os.path.basename(path))
    if not match:
        raise ValueError(f"Unable to infer chromosome from filename: {path}")
    return match.group("chr")


def ancestry_code_to_label(token: str, code_to_label: Dict[str, str]) -> str:
    match = ANC_TOKEN_RE.match(token.strip())
    if not match:
        return token
    code = match.group("code")
    if code not in code_to_label:
        raise ValueError(f"Observed ancestry token '{token}' is not present in configured mapping")
    return code_to_label[code]


def merge_q_files(q_files: List[str], code_to_label: Dict[str, str], output_path: str) -> None:
    expected_labels = [code_to_label[str(code)] for code in sorted(int(code) for code in code_to_label)]
    wrote_header = False

    with open(output_path, "w", encoding="utf-8", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")

        for path in q_files:
            chromosome = infer_chr(path)
            header = None

            with open(path, "r", encoding="utf-8") as in_handle:
                for raw_line in in_handle:
                    line = raw_line.strip()
                    if not line:
                        continue
                    if line.startswith("#sample"):
                        parts = line.lstrip("#").split()
                        if len(parts) < 2:
                            raise ValueError(f"Malformed Q header in {path}: {line}")
                        mapped_labels = [ancestry_code_to_label(token, code_to_label) for token in parts[1:]]
                        if mapped_labels != expected_labels:
                            raise ValueError(
                                f"Q header ancestry order mismatch in {path}. "
                                f"Observed {mapped_labels}; expected {expected_labels}"
                            )
                        header = ["chr", "sample", *mapped_labels]
                        if not wrote_header:
                            writer.writerow(header)
                            wrote_header = True
                        continue
                    if line.startswith("#"):
                        continue
                    if header is None:
                        raise ValueError(f"Missing '#sample' header before data in {path}")
                    parts = line.split()
                    if len(parts) != len(header) - 1:
                        raise ValueError(
                            f"Q row column count mismatch in {path}. "
                            f"Expected {len(header) - 1}; observed {len(parts)}"
                        )
                    writer.writerow([chromosome, *parts])

    if not wrote_header:
        raise ValueError("No Q records were written to merged output")


def build_msp_header(raw_header: str) -> List[str]:
    parts = raw_header.lstrip("#").split()
    if len(parts) < 8:
        raise ValueError(f"Malformed MSP header: {raw_header}")
    first_cols = ["chm", "spos", "epos", "sgpos", "egpos", "n_snps"]
    sample_cols = parts[7:]
    return ["chr", *first_cols, *sample_cols]


def merge_msp_files(msp_files: List[str], code_to_label: Dict[str, str], output_path: str) -> None:
    header = None

    with open(output_path, "w", encoding="utf-8", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")

        for path in msp_files:
            chromosome = infer_chr(path)
            local_header = None
            with open(path, "r", encoding="utf-8") as in_handle:
                for raw_line in in_handle:
                    line = raw_line.strip()
                    if not line:
                        continue
                    if line.startswith("#chm"):
                        local_header = build_msp_header(line)
                        if header is None:
                            header = local_header
                            writer.writerow(header)
                        elif local_header != header:
                            raise ValueError(
                                f"MSP header mismatch in {path}. "
                                f"Observed {local_header[:10]}...; expected {header[:10]}..."
                            )
                        continue
                    if line.startswith("#"):
                        continue
                    if local_header is None:
                        raise ValueError(f"Missing '#chm' header before MSP data in {path}")

                    parts = line.split()
                    expected_columns = len(local_header) - 1
                    if len(parts) != expected_columns:
                        raise ValueError(
                            f"MSP row column count mismatch in {path}. "
                            f"Expected {expected_columns}; observed {len(parts)}"
                        )

                    fixed = parts[:6]
                    ancestry_calls = [ancestry_code_to_label(token, code_to_label) for token in parts[6:]]
                    writer.writerow([chromosome, *fixed, *ancestry_calls])

    if header is None:
        raise ValueError("No MSP records were written to merged output")


def validate_output(path: str) -> None:
    with open(path, "r", encoding="utf-8") as handle:
        first = handle.readline().strip()
        if not first:
            raise ValueError(f"Merged output is empty: {path}")
        if first.startswith("#"):
            raise ValueError(f"Merged output still contains a comment header: {path}")
        if first.split("\t", 1)[0] != "chr":
            raise ValueError(f"Merged output is missing the leading chr column: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge per-chromosome RFMix Q and MSP outputs")
    parser.add_argument("--q-files", nargs="+", required=True)
    parser.add_argument("--msp-files", nargs="+", required=True)
    parser.add_argument("--label-code-map", required=True)
    parser.add_argument("--output-q", required=True)
    parser.add_argument("--output-msp", required=True)
    args = parser.parse_args()

    label_to_code = parse_label_code_map(args.label_code_map)
    code_to_label = {code: label for label, code in label_to_code.items()}

    merge_q_files(args.q_files, code_to_label, args.output_q)
    merge_msp_files(args.msp_files, code_to_label, args.output_msp)
    validate_output(args.output_q)
    validate_output(args.output_msp)

    print(f"Merged {len(args.q_files)} Q files into {args.output_q}")
    print(f"Merged {len(args.msp_files)} MSP files into {args.output_msp}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"merge_rfmix_outputs failed: {exc}", file=sys.stderr)
        sys.exit(1)
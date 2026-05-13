#!/usr/bin/env python3
import argparse
import csv
import re


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Select best ADMIXTURE K using CV error.")
    parser.add_argument("--k-logs", required=True, help="Comma-separated list of K=path_to_log")
    parser.add_argument("--output-summary", required=True, help="Output summary TSV")
    parser.add_argument("--output-best-k", required=True, help="Output file containing best K")
    return parser.parse_args()


def parse_cv_error(log_path: str, k: int) -> float:
    cv_value = None
    pattern = re.compile(r"CV\s+error\s*\(K\s*=\s*([0-9]+)\)\s*:\s*([0-9eE+\-.]+)")
    with open(log_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            match = pattern.search(line)
            if match:
                k_found = int(match.group(1))
                if k_found == k:
                    cv_value = float(match.group(2))
    if cv_value is None:
        raise ValueError(f"Could not parse CV error for K={k} from log: {log_path}")
    return cv_value


def main() -> None:
    args = parse_args()

    entries = []
    for token in args.k_logs.split(","):
        token = token.strip()
        if not token:
            continue
        k_str, path = token.split("=", 1)
        k = int(k_str)
        cv = parse_cv_error(path, k)
        entries.append((k, cv, path))

    if not entries:
        raise ValueError("No K logs were provided")

    entries.sort(key=lambda x: x[0])
    best_k, best_cv, _ = min(entries, key=lambda x: x[1])

    with open(args.output_summary, "w", encoding="utf-8", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["K", "CV_ERROR", "LOG_FILE"])
        for k, cv, path in entries:
            writer.writerow([k, cv, path])

    with open(args.output_best_k, "w", encoding="utf-8") as out_handle:
        out_handle.write(str(best_k) + "\n")

    print(f"Best K: {best_k} (CV={best_cv})")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create population sample lists and RFMix sample map from 1000G panel"
    )
    parser.add_argument("--panel", required=True)
    parser.add_argument(
        "--population",
        action="append",
        required=True,
        metavar="CODE:LABEL:OUT",
        help="Repeat for each population (e.g. --population CHB:EAS:results/CHB.samples.txt)",
    )
    parser.add_argument("--combined", required=True)
    parser.add_argument("--sample-map", dest="sample_map", required=True)
    return parser.parse_args()


def load_panel(panel_path: Path):
    with panel_path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            raise ValueError("Panel file is empty")
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            yield parts[0], parts[1]


def write_lines(path: Path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(f"{line}\n")


def main():
    args = parse_args()

    pop_configs = []
    for spec in args.population:
        parts = str(spec).split(":", 2)
        if len(parts) != 3:
            raise ValueError(
                f"Invalid --population spec '{spec}'. Expected CODE:LABEL:OUT"
            )
        code, label, out_path = parts
        code = code.strip()
        label = label.strip()
        out_path = out_path.strip()
        if not code or not label or not out_path:
            raise ValueError(
                f"Invalid --population spec '{spec}'. Empty CODE, LABEL, or OUT"
            )
        pop_configs.append((code, label, Path(out_path)))

    if not pop_configs:
        raise ValueError("At least one --population is required")

    pop_stacks: dict[str, list[str]] = {code: [] for code, _, _ in pop_configs}

    panel_path = Path(args.panel)
    if not panel_path.exists():
        raise FileNotFoundError(f"Panel file not found: {panel_path}")

    for sample_id, pop in load_panel(panel_path):
        if pop in pop_stacks:
            pop_stacks[pop].append(sample_id)

    for code, label, out_path in pop_configs:
        samples = pop_stacks[code]
        if not samples:
            raise ValueError(f"No {code} samples found in panel")
        write_lines(out_path, samples)
        print(f"{code} ({label}) count: {len(samples)}")

    combined = sorted(set(s for samples in pop_stacks.values() for s in samples))
    write_lines(Path(args.combined), combined)
    print(f"Combined unique count: {len(combined)}")

    sample_map_lines = []
    for code, label, _ in pop_configs:
        for sample in pop_stacks[code]:
            sample_map_lines.append(f"{sample} {label}")
    write_lines(Path(args.sample_map), sample_map_lines)


if __name__ == "__main__":
    main()


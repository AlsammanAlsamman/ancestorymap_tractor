#!/usr/bin/env python3
import argparse
from bisect import bisect_left
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def chr_sort_key(value: str) -> tuple[int, object]:
    value = str(value).replace("chr", "")
    return (0, int(value)) if value.isdigit() else (1, value)


def parse_args():
    parser = argparse.ArgumentParser(description="Merge cohort/reference MAF tables and annotate SNPs")
    parser.add_argument(
        "--source",
        action="append",
        required=True,
        help="Source table in the form LABEL=/path/to/file.tsv (repeat for each source)",
    )
    parser.add_argument("--annotation", required=True, help="Gene annotation file from analysis.yml")
    parser.add_argument("--loci-file", required=True, help="Loci definition file from analysis.yml")
    parser.add_argument("--output", required=True, help="Output merged annotated TSV")
    return parser.parse_args()


def first_non_missing(values) -> str:
    for value in values:
        if pd.notna(value):
            text = str(value).strip()
            if text and text not in {".", "nan", "None"}:
                return text
    return "."


def _process_chunk(chunk: "pd.DataFrame", label: str, filter_keys: "set | None") -> "pd.DataFrame":
    """Normalise one chunk and optionally filter to a set of position keys."""
    if "rsid" not in chunk.columns:
        chunk = chunk.copy()
        chunk["rsid"] = "."

    required = {"chr", "pos", "ref", "alt", "alt_af", "alt_ac", "an"}
    missing = required.difference(chunk.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    chunk = chunk.copy()
    chunk["chr"] = chunk["chr"].astype(str).str.replace("chr", "", regex=False)
    chunk["pos"] = pd.to_numeric(chunk["pos"], errors="raise")
    chunk["ref"] = chunk["ref"].astype(str).str.upper()
    chunk["alt"] = chunk["alt"].astype(str).str.upper()
    chunk["allele_a"] = chunk[["ref", "alt"]].min(axis=1)
    chunk["allele_b"] = chunk[["ref", "alt"]].max(axis=1)

    if filter_keys is not None:
        key_col = chunk["chr"] + "_" + chunk["pos"].astype(str) + "_" + chunk["allele_a"] + "_" + chunk["allele_b"]
        chunk = chunk[key_col.isin(filter_keys)]

    if chunk.empty:
        return chunk

    keep = chunk[["chr", "pos", "allele_a", "allele_b"]].copy()
    keep[f"{label}_rsid"] = chunk["rsid"].fillna(".").astype(str)
    keep[f"{label}_source_ref"] = chunk["ref"].values
    keep[f"{label}_source_alt"] = chunk["alt"].values
    keep[f"{label}_source_alt_af"] = pd.to_numeric(chunk["alt_af"], errors="coerce").values
    keep[f"{label}_source_alt_ac"] = pd.to_numeric(chunk["alt_ac"], errors="coerce").values
    keep[f"{label}_an"] = pd.to_numeric(chunk["an"], errors="coerce").values
    if "n_samples_called" in chunk.columns:
        keep[f"{label}_n_samples_called"] = pd.to_numeric(chunk["n_samples_called"], errors="coerce").values
    else:
        keep[f"{label}_n_samples_called"] = pd.NA
    return keep


def load_source_table(label: str, path: str, filter_keys: "set | None" = None) -> pd.DataFrame:
    """Load a MAF table in chunks, optionally filtering rows to *filter_keys*.

    *filter_keys* is a set of strings of the form ``chr_pos_allele_a_allele_b``.
    When provided (for large reference-panel files) only matching rows are kept,
    dramatically reducing peak memory usage.
    """
    chunks = []
    for chunk in pd.read_csv(path, sep="\t", dtype={"chr": str}, chunksize=500_000):
        processed = _process_chunk(chunk, label, filter_keys)
        if not processed.empty:
            chunks.append(processed)

    if not chunks:
        # Return a correctly-typed empty DataFrame so downstream merges work.
        cols = ["chr", "pos", "allele_a", "allele_b",
                f"{label}_rsid", f"{label}_source_ref", f"{label}_source_alt",
                f"{label}_source_alt_af", f"{label}_source_alt_ac",
                f"{label}_an", f"{label}_n_samples_called"]
        return pd.DataFrame(columns=cols)

    df = pd.concat(chunks, ignore_index=True)
    return df.drop_duplicates(subset=["chr", "pos", "allele_a", "allele_b"], keep="first")


def load_loci(loci_path: str) -> Dict[str, List[Tuple[int, int, str]]]:
    loci_df = pd.read_csv(loci_path, sep=r"\s+", engine="python")
    rename_map = {col.lower(): col for col in loci_df.columns}
    required = {"locus", "chr", "start", "end"}
    if not required.issubset(rename_map):
        raise ValueError(f"Loci file must contain columns: {sorted(required)}")

    loci_by_chr: Dict[str, List[Tuple[int, int, str]]] = {}
    for _, row in loci_df.iterrows():
        chrom = str(row[rename_map["chr"]]).replace("chr", "")
        start = int(row[rename_map["start"]])
        end = int(row[rename_map["end"]])
        locus = str(row[rename_map["locus"]])
        loci_by_chr.setdefault(chrom, []).append((start, end, locus))
    return loci_by_chr


def load_genes(annotation_path: str) -> Dict[str, List[Tuple[int, int, str]]]:
    genes_df = pd.read_csv(annotation_path, sep=r"\s+", header=None, engine="python", comment="#")
    if genes_df.shape[1] < 6:
        raise ValueError(f"Annotation file must have at least 6 columns: {annotation_path}")

    genes_by_chr: Dict[str, List[Tuple[int, int, str]]] = {}
    for _, row in genes_df.iterrows():
        gene_name = str(row.iloc[5]).strip()
        if not gene_name or "-" in gene_name:
            continue
        chrom = str(row.iloc[1]).replace("chr", "")
        start = int(row.iloc[2])
        end = int(row.iloc[3])
        genes_by_chr.setdefault(chrom, []).append((start, end, gene_name))

    for chrom in genes_by_chr:
        genes_by_chr[chrom].sort(key=lambda item: (item[0], item[1], item[2]))
    return genes_by_chr


def assign_locus(chrom: str, pos: int, loci_by_chr: Dict[str, List[Tuple[int, int, str]]]) -> str:
    for start, end, locus_name in loci_by_chr.get(chrom, []):
        if start <= pos <= end:
            return locus_name
    return ""


def nearest_gene(chrom: str, pos: int, genes_by_chr: Dict[str, List[Tuple[int, int, str]]]) -> Tuple[str, int]:
    genes = genes_by_chr.get(chrom, [])
    if not genes:
        return "", -1

    starts = [item[0] for item in genes]
    idx = bisect_left(starts, pos)

    candidate_indices = set()
    for offset in range(-3, 4):
        candidate_idx = idx + offset
        if 0 <= candidate_idx < len(genes):
            candidate_indices.add(candidate_idx)

    scan_idx = idx - 1
    while scan_idx >= 0 and genes[scan_idx][1] >= pos:
        candidate_indices.add(scan_idx)
        scan_idx -= 1

    best_gene = ""
    best_distance = -1
    for candidate_idx in sorted(candidate_indices):
        start, end, gene_name = genes[candidate_idx]
        if start <= pos <= end:
            distance = 0
        else:
            distance = min(abs(pos - start), abs(pos - end))
        if best_distance == -1 or distance < best_distance or (distance == best_distance and gene_name < best_gene):
            best_gene = gene_name
            best_distance = distance
    return best_gene, best_distance


def apply_harmonized_alt_frequency(merged: pd.DataFrame, label: str) -> None:
    ref_col = f"{label}_source_ref"
    alt_col = f"{label}_source_alt"
    af_col = f"{label}_source_alt_af"

    source_af = pd.to_numeric(merged[af_col], errors="coerce")

    same_mask = merged[ref_col].eq(merged["ref"]) & merged[alt_col].eq(merged["alt"])
    flip_mask = merged[ref_col].eq(merged["alt"]) & merged[alt_col].eq(merged["ref"])

    harmonized_af = pd.Series(float("nan"), index=merged.index, dtype="float64")
    harmonized_af.loc[same_mask] = source_af.loc[same_mask]
    harmonized_af.loc[flip_mask] = 1.0 - source_af.loc[flip_mask]

    merged[f"{label}_maf"] = harmonized_af.where(harmonized_af <= 0.5, 1.0 - harmonized_af).round(8)


def main():
    args = parse_args()

    # Parse all source specs first so we can load cohort before reference panels.
    source_specs = []
    for source_spec in args.source:
        if "=" not in source_spec:
            raise ValueError(f"Invalid --source value '{source_spec}'. Use LABEL=/path/to/file.tsv")
        label, path = source_spec.split("=", 1)
        source_specs.append((label.strip(), path.strip()))

    # ------------------------------------------------------------------
    # Load cohort table first (small: loci-region variants only).
    # Build a set of position keys used to filter the large reference-
    # panel files so we never load all 77 M rows per file into RAM.
    # ------------------------------------------------------------------
    if not source_specs:
        raise ValueError("No --source arguments were provided.")

    cohort_label, cohort_path = source_specs[0]
    print(f"Loading cohort source: {cohort_label} from {cohort_path}", flush=True)
    cohort_df = load_source_table(cohort_label, cohort_path, filter_keys=None)
    cohort_keys: set = set(
        cohort_df["chr"].astype(str)
        + "_"
        + cohort_df["pos"].astype(str)
        + "_"
        + cohort_df["allele_a"].astype(str)
        + "_"
        + cohort_df["allele_b"].astype(str)
    )
    print(f"  → {len(cohort_df):,} cohort variants; built position-key filter set.", flush=True)

    merged = cohort_df
    source_labels = [cohort_label]

    for label, path in source_specs[1:]:
        print(f"Loading reference-panel source: {label} from {path} (filtered to cohort positions)", flush=True)
        source_df = load_source_table(label, path, filter_keys=cohort_keys)
        print(f"  → {len(source_df):,} variants retained after position filter.", flush=True)
        merged = merged.merge(source_df, on=["chr", "pos", "allele_a", "allele_b"], how="left")
        source_labels.append(label)

    if merged.empty:
        raise ValueError("No MAF tables were provided for merging")

    # Keep source order as provided by the workflow (cohort first, then configured ancestries).
    source_labels = list(dict.fromkeys(source_labels))

    merged["chr"] = merged["chr"].astype(str).str.replace("chr", "", regex=False)
    merged["pos"] = pd.to_numeric(merged["pos"], errors="raise")

    merged["ref"] = pd.NA
    merged["alt"] = pd.NA
    for label in source_labels:
        ref_col = f"{label}_source_ref"
        alt_col = f"{label}_source_alt"
        if ref_col in merged.columns and alt_col in merged.columns:
            mask = merged["ref"].isna() & merged[ref_col].notna() & merged[alt_col].notna()
            merged.loc[mask, "ref"] = merged.loc[mask, ref_col]
            merged.loc[mask, "alt"] = merged.loc[mask, alt_col]

    merged["ref"] = merged["ref"].fillna(merged["allele_a"])
    merged["alt"] = merged["alt"].fillna(merged["allele_b"])

    rsid_cols = [f"{label}_rsid" for label in source_labels if f"{label}_rsid" in merged.columns]
    merged["rsid"] = merged[rsid_cols].apply(first_non_missing, axis=1) if rsid_cols else "."
    merged = merged[merged["rsid"].astype(str).str.strip().ne(".") & merged["rsid"].astype(str).str.strip().ne("")].copy()

    for label in source_labels:
        apply_harmonized_alt_frequency(merged, label)

    loci_by_chr = load_loci(args.loci_file)
    genes_by_chr = load_genes(args.annotation)

    merged["locus_name"] = [assign_locus(chrom, pos, loci_by_chr) for chrom, pos in zip(merged["chr"], merged["pos"])]
    nearest = [nearest_gene(chrom, pos, genes_by_chr) for chrom, pos in zip(merged["chr"], merged["pos"])]
    merged["gene"] = [item[0] for item in nearest]
    merged["gene_distance_bp"] = [item[1] for item in nearest]

    ordered_columns = ["chr", "pos", "rsid", "locus_name", "gene", "gene_distance_bp", "ref", "alt"]
    for label in source_labels:
        ordered_columns.append(f"{label}_maf")

    merged = merged[ordered_columns].sort_values(
        by=["chr", "pos"],
        key=lambda col: col.map(chr_sort_key) if col.name == "chr" else col,
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote merged annotated MAF table to {output_path} with {len(merged)} SNPs")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Merge harmonized MAF table with merged Tractor GWAS pairwise tables")
    parser.add_argument("--maf", required=True, help="Merged MAF table path")
    parser.add_argument(
        "--gwas",
        action="append",
        required=True,
        help="GWAS input in the form PAIR=/path/to/merged.model_PAIR.gwas.tsv",
    )
    parser.add_argument("--output", required=True, help="Output merged summary TSV")
    return parser.parse_args()


def load_maf_table(path: str, ancestry_order: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"rsid": str, "ref": str, "alt": str})
    required = [
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
    required.extend([f"{ancestry}_maf" for ancestry in ancestry_order])
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"MAF table missing required columns: {missing}")

    df = df[required].copy()
    df["chr"] = df["chr"].astype(str).str.replace("chr", "", regex=False)
    df["pos"] = pd.to_numeric(df["pos"], errors="raise")
    df["rsid"] = df["rsid"].astype(str).str.strip()
    df = df[df["rsid"].ne("") & df["rsid"].ne(".")].copy()
    df["ref"] = df["ref"].astype(str).str.upper()
    df["alt"] = df["alt"].astype(str).str.upper()
    return df.drop_duplicates(subset=["chr", "pos", "rsid", "ref", "alt"], keep="first")


def load_gwas_table(pair_name: str, path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"ID": str, "REF": str, "ALT": str, "ancestry_i": str, "ancestry_j": str})
    required = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "beta_local_ancestry",
        "beta_dosage_i",
        "beta_dosage_j",
        "se_dosage_i",
        "se_dosage_j",
        "p_local_ancestry",
        "p_dosage_i",
        "p_dosage_j",
        "ancestry_i",
        "ancestry_j",
    ]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"GWAS table for {pair_name} missing required columns: {missing}")

    df = df[required].copy()
    df = df.rename(columns={"CHROM": "chr", "POS": "pos", "ID": "rsid", "REF": "gwas_ref", "ALT": "gwas_alt"})
    df["chr"] = df["chr"].astype(str).str.replace("chr", "", regex=False)
    df["pos"] = pd.to_numeric(df["pos"], errors="raise")
    df["rsid"] = df["rsid"].astype(str).str.strip()
    df = df[df["rsid"].ne("") & df["rsid"].ne(".")].copy()
    df["gwas_ref"] = df["gwas_ref"].astype(str).str.upper()
    df["gwas_alt"] = df["gwas_alt"].astype(str).str.upper()
    return df.drop_duplicates(subset=["chr", "pos", "rsid", "gwas_ref", "gwas_alt"], keep="first")


def add_pair_metrics(summary: pd.DataFrame, pair_name: str, gwas_df: pd.DataFrame, ancestry_order: list[str]) -> pd.DataFrame:
    merged = summary.merge(gwas_df, on=["chr", "pos", "rsid"], how="left")

    same_mask = merged["gwas_ref"].eq(merged["ref"]) & merged["gwas_alt"].eq(merged["alt"])
    flip_mask = merged["gwas_ref"].eq(merged["alt"]) & merged["gwas_alt"].eq(merged["ref"])
    orientation_mask = same_mask | flip_mask

    beta_local = pd.to_numeric(merged["beta_local_ancestry"], errors="coerce")
    harmonized_local_beta = pd.Series(np.nan, index=merged.index, dtype="float64")
    harmonized_local_beta.loc[same_mask] = beta_local.loc[same_mask]
    harmonized_local_beta.loc[flip_mask] = -beta_local.loc[flip_mask]
    merged[f"OR_local_ancestry_{pair_name}"] = np.exp(harmonized_local_beta)

    local_p = pd.Series(np.nan, index=merged.index, dtype="float64")
    local_p.loc[orientation_mask] = pd.to_numeric(merged.loc[orientation_mask, "p_local_ancestry"], errors="coerce")
    merged[f"p_local_ancestry_{pair_name}"] = local_p

    for role in ["i", "j"]:
        ancestry_col = f"ancestry_{role}"
        beta_col = f"beta_dosage_{role}"
        p_col = f"p_dosage_{role}"
        se_col = f"se_dosage_{role}"

        beta = pd.to_numeric(merged[beta_col], errors="coerce")
        harmonized_beta = pd.Series(np.nan, index=merged.index, dtype="float64")
        harmonized_beta.loc[same_mask] = beta.loc[same_mask]
        harmonized_beta.loc[flip_mask] = -beta.loc[flip_mask]
        odds_ratio = np.exp(harmonized_beta)

        for ancestry in ancestry_order:
            mask = merged[ancestry_col].astype(str).eq(ancestry) & orientation_mask
            merged.loc[mask, f"__p_dosage_{ancestry}_{pair_name}"] = pd.to_numeric(merged.loc[mask, p_col], errors="coerce")
            merged.loc[mask, f"__beta_dosage_{ancestry}_{pair_name}"] = harmonized_beta.loc[mask]
            merged.loc[mask, f"__OR_dosage_{ancestry}_{pair_name}"] = odds_ratio.loc[mask]
            merged.loc[mask, f"__se_dosage_{ancestry}_{pair_name}"] = pd.to_numeric(merged.loc[mask, se_col], errors="coerce")

    drop_cols = [
        "gwas_ref",
        "gwas_alt",
        "beta_local_ancestry",
        "beta_dosage_i",
        "beta_dosage_j",
        "se_dosage_i",
        "se_dosage_j",
        "p_local_ancestry",
        "p_dosage_i",
        "p_dosage_j",
        "ancestry_i",
        "ancestry_j",
    ]
    return merged.drop(columns=[col for col in drop_cols if col in merged.columns])


def mean_from_columns(df: pd.DataFrame, columns: list[str]) -> pd.Series:
    available = [col for col in columns if col in df.columns]
    if not available:
        return pd.Series(np.nan, index=df.index, dtype="float64")
    return df[available].apply(pd.to_numeric, errors="coerce").mean(axis=1, skipna=True)


def main():
    args = parse_args()

    gwas_specs = {}
    for item in args.gwas:
        if "=" not in item:
            raise ValueError(f"Invalid --gwas value '{item}'. Use PAIR=/path/to/file.tsv")
        pair_name, path = item.split("=", 1)
        pair_name = pair_name.strip()
        gwas_specs[pair_name] = path.strip()

    pair_order = list(gwas_specs.keys())
    ancestry_order = []
    for pair_name in pair_order:
        parts = pair_name.split("_")
        if len(parts) != 2:
            raise ValueError(f"GWAS pair '{pair_name}' must be in the form ANC1_ANC2")
        for ancestry in parts:
            if ancestry not in ancestry_order:
                ancestry_order.append(ancestry)

    summary = load_maf_table(args.maf, ancestry_order)

    for pair_name in pair_order:
        gwas_df = load_gwas_table(pair_name, gwas_specs[pair_name])
        summary = add_pair_metrics(summary, pair_name, gwas_df, ancestry_order)

    for ancestry in ancestry_order:
        summary[f"p_dosage_{ancestry}"] = mean_from_columns(
            summary,
            [f"__p_dosage_{ancestry}_{pair_name}" for pair_name in pair_order],
        )
        mean_beta = mean_from_columns(
            summary,
            [f"__beta_dosage_{ancestry}_{pair_name}" for pair_name in pair_order],
        )
        summary[f"OR_dosage_{ancestry}"] = np.exp(mean_beta)
        summary[f"se_dosage_{ancestry}"] = mean_from_columns(
            summary,
            [f"__se_dosage_{ancestry}_{pair_name}" for pair_name in pair_order],
        )

    final_columns = [
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
    final_columns.extend([f"{ancestry}_maf" for ancestry in ancestry_order])
    final_columns.extend([f"OR_local_ancestry_{pair_name}" for pair_name in pair_order])
    final_columns.extend([f"p_local_ancestry_{pair_name}" for pair_name in pair_order])
    final_columns.extend([f"p_dosage_{ancestry}" for ancestry in ancestry_order])
    final_columns.extend([f"OR_dosage_{ancestry}" for ancestry in ancestry_order])
    final_columns.extend([f"se_dosage_{ancestry}" for ancestry in ancestry_order])

    required_metric_columns = []
    required_metric_columns.extend([f"OR_local_ancestry_{pair_name}" for pair_name in pair_order])
    required_metric_columns.extend([f"p_local_ancestry_{pair_name}" for pair_name in pair_order])
    required_metric_columns.extend([f"p_dosage_{ancestry}" for ancestry in ancestry_order])
    required_metric_columns.extend([f"OR_dosage_{ancestry}" for ancestry in ancestry_order])
    # Keep partially informative SNPs; only remove rows where every required metric is missing.
    summary = summary.dropna(subset=required_metric_columns, how="all")
    summary = summary[final_columns].sort_values(by=["chr", "pos", "rsid"])

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_path, sep="\t", index=False, na_rep="")
    print(f"Wrote merged MAF+GWAS summary to {output_path} with {len(summary)} rows")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import glob
import os
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def extract_base_iid(sample_id):
    x = str(sample_id)
    underscores = [i for i, ch in enumerate(x) if ch == "_"]
    for idx in underscores:
        left = x[:idx]
        right = x[idx + 1 :]
        if left == right:
            return left
    return x


def run_cmd(cmd, cwd=None):
    result = subprocess.run(cmd, cwd=cwd, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(cmd)}\n"
            f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    return result


def chr_sort_key(ch):
    s = str(ch).replace("chr", "")
    return (0, int(s)) if s.isdigit() else (1, s)


def load_psam_iids(psam_path):
    psam = pd.read_csv(psam_path, sep="\t", dtype=str)
    iid_col = "IID" if "IID" in psam.columns else "#IID"
    if iid_col not in psam.columns:
        raise ValueError(f"Could not find IID column in {psam_path}")
    return psam[iid_col].astype(str).tolist()


def build_sample_map(sample_ids):
    mapping = {}
    for sid in sample_ids:
        mapping[extract_base_iid(sid)] = sid
    return mapping


def prepare_phenotype(args, sample_map, out_path):
    if args.phenotype_format == "fam":
        fam = pd.read_csv(args.phenotype, sep=r"\s+", header=None, dtype=str, engine="python")
        if fam.shape[1] < 6:
            raise ValueError(f"FAM must have >=6 columns: {args.phenotype}")
        raw = pd.DataFrame({"IID_BASE": fam.iloc[:, 1].astype(str), "PHENO": fam.iloc[:, 5].astype(str)})
    else:
        tab = pd.read_csv(args.phenotype, sep=None, dtype=str, engine="python")
        if args.table_iid_column not in tab.columns or args.table_phenotype_column not in tab.columns:
            raise ValueError("Phenotype table missing IID/PHENO configured columns")
        raw = pd.DataFrame(
            {
                "IID_BASE": tab[args.table_iid_column].astype(str),
                "PHENO": tab[args.table_phenotype_column].astype(str),
            }
        )

    raw["IID"] = raw["IID_BASE"].map(sample_map)
    pheno = raw.dropna(subset=["IID", "PHENO"]).copy()
    pheno["FID"] = pheno["IID"]
    pheno = pheno[["FID", "IID", "PHENO"]].drop_duplicates(subset=["IID"])
    pheno.to_csv(out_path, sep="\t", index=False)
    if pheno.empty:
        raise ValueError("No phenotype rows matched merged PLINK sample IDs")


def prepare_covariates(args, sample_map, out_path):
    if args.covariates_format != "plink_eigenvec":
        raise ValueError(f"Unsupported covariates format: {args.covariates_format}")

    has_header = parse_bool(args.covariates_has_header)
    if has_header:
        cov = pd.read_csv(args.covariates_file, sep=r"\s+", dtype=str, engine="python")
    else:
        cov = pd.read_csv(args.covariates_file, sep=r"\s+", header=None, dtype=str, engine="python")
        if cov.shape[1] < 3:
            raise ValueError("Covariates file needs at least FID IID PC1")
        cov.columns = [args.covariates_fid_column, args.covariates_iid_column] + [
            f"PC{i}" for i in range(1, cov.shape[1] - 1)
        ]

    pc_cols = [c for c in cov.columns if str(c).upper().startswith("PC")][: args.n_pcs]
    if len(pc_cols) < args.n_pcs:
        raise ValueError(f"Requested {args.n_pcs} PCs; found {len(pc_cols)}")

    out = pd.DataFrame()
    out["IID_BASE"] = cov[args.covariates_iid_column].astype(str)
    out["IID"] = out["IID_BASE"].map(sample_map)
    out = out.dropna(subset=["IID"]).copy()
    out["FID"] = out["IID"]
    for col in pc_cols:
        out[col] = pd.to_numeric(cov.loc[out.index, col], errors="coerce")

    out = out[["FID", "IID"] + pc_cols].drop_duplicates(subset=["IID"])
    out.to_csv(out_path, sep="\t", index=False)
    if out.empty:
        raise ValueError("No covariate rows matched merged PLINK sample IDs")
    return pc_cols


def choose_glm_file(prefix):
    candidates = sorted(glob.glob(prefix + ".PHENO1.glm.*"))
    if not candidates:
        raise FileNotFoundError(f"No PLINK GLM output found for prefix: {prefix}")
    for path in candidates:
        if ".logistic" in path:
            return path
    return candidates[0]


def make_gwas_table(glm_path, out_tsv):
    g = pd.read_csv(glm_path, sep=r"\s+", dtype=str, engine="python")
    if "TEST" in g.columns:
        g = g[g["TEST"] == "ADD"].copy()

    chr_col = "#CHROM" if "#CHROM" in g.columns else ("CHROM" if "CHROM" in g.columns else None)
    if chr_col is None or "POS" not in g.columns or "ID" not in g.columns or "P" not in g.columns:
        raise ValueError(f"Unexpected GLM columns in {glm_path}")

    g["P"] = pd.to_numeric(g["P"], errors="coerce")
    g = g[g["P"].notna() & (g["P"] > 0)].copy()
    g["CHR"] = g[chr_col].astype(str).str.replace("chr", "", regex=False)
    g["POS"] = pd.to_numeric(g["POS"], errors="coerce")
    g = g[g["POS"].notna()].copy()

    keep = ["CHR", "POS", "ID", "P"]
    for extra in ["A1", "A1_FREQ", "BETA", "SE", "OR", "Z_STAT", "T_STAT", "ERRCODE"]:
        if extra in g.columns:
            keep.append(extra)

    out = g[keep].copy()
    out.to_csv(out_tsv, sep="\t", index=False)
    return out


def plot_manhattan(gwas_df, png_path, title):
    if gwas_df.empty:
        plt.figure(figsize=(12, 5))
        plt.text(0.5, 0.5, "No GWAS points to plot", ha="center", va="center")
        plt.axis("off")
        plt.savefig(png_path, dpi=180, bbox_inches="tight")
        plt.close()
        return

    df = gwas_df.copy()
    df["CHR"] = df["CHR"].astype(str)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df.dropna(subset=["POS", "P"])
    df = df[df["P"] > 0]
    if df.empty:
        plt.figure(figsize=(12, 5))
        plt.text(0.5, 0.5, "No valid GWAS p-values to plot", ha="center", va="center")
        plt.axis("off")
        plt.savefig(png_path, dpi=180, bbox_inches="tight")
        plt.close()
        return

    chr_order = sorted(df["CHR"].unique(), key=chr_sort_key)
    offsets = {}
    running = 0
    centers = []
    for ch in chr_order:
        ch_max = df.loc[df["CHR"] == ch, "POS"].max()
        offsets[ch] = running
        centers.append((ch, running + ch_max / 2.0))
        running += ch_max + 1_000_000

    df["x"] = df.apply(lambda r: r["POS"] + offsets[r["CHR"]], axis=1)
    df["y"] = -np.log10(df["P"].astype(float))

    plt.figure(figsize=(14, 6))
    colors = ["#1f77b4", "#ff7f0e"]
    for i, ch in enumerate(chr_order):
        sub = df[df["CHR"] == ch]
        plt.scatter(sub["x"], sub["y"], s=8, c=colors[i % 2], alpha=0.75, linewidths=0)

    plt.axhline(-np.log10(5e-8), color="#cc2f2f", linestyle="--", linewidth=1)
    plt.title(title)
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(P)")
    plt.xticks([c for _, c in centers], [ch for ch, _ in centers], rotation=0)
    plt.tight_layout()
    plt.savefig(png_path, dpi=180)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Run ancestry-specific merged PLINK GWAS and Manhattan plot")
    parser.add_argument("--ancestry", required=True)
    parser.add_argument("--chromosomes", required=True)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--merged-prefix", required=True)
    parser.add_argument("--gwas-prefix", required=True)
    parser.add_argument("--gwas-tsv", required=True)
    parser.add_argument("--manhattan-png", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--phenotype-format", choices=["fam", "table"], required=True)
    parser.add_argument("--fam-case-code", default="2")
    parser.add_argument("--fam-control-code", default="1")
    parser.add_argument("--table-iid-column", default="IID")
    parser.add_argument("--table-phenotype-column", default="PHENO")
    parser.add_argument("--covariates-file", required=True)
    parser.add_argument("--covariates-format", default="plink_eigenvec")
    parser.add_argument("--covariates-has-header", default="false")
    parser.add_argument("--covariates-fid-column", default="FID")
    parser.add_argument("--covariates-iid-column", default="IID")
    parser.add_argument("--n-pcs", type=int, default=5)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    chr_list = [c for c in str(args.chromosomes).split(",") if c]

    pfile_triples = []
    for ch in chr_list:
        prefix = os.path.join(args.input_dir, f"chr{ch}.dosage")
        pgen = prefix + ".pgen"
        pvar = prefix + ".pvar"
        psam = prefix + ".psam"
        if os.path.exists(pgen) and os.path.exists(pvar) and os.path.exists(psam):
            pfile_triples.append((pgen, pvar, psam, prefix))
    if not pfile_triples:
        raise FileNotFoundError(f"No chromosome PLINK dosage files found in {args.input_dir}")

    first_prefix = pfile_triples[0][3]
    if len(pfile_triples) == 1:
        run_cmd(
            [
                "plink2",
                "--pfile",
                first_prefix,
                "--allow-extra-chr",
                "--make-pgen",
                "--out",
                args.merged_prefix,
            ]
        )
    else:
        merge_list = os.path.join(args.output_dir, "merge_list.txt")
        with open(merge_list, "w", encoding="utf-8") as handle:
            for pgen, pvar, psam, _ in pfile_triples[1:]:
                handle.write(f"{pgen} {pvar} {psam}\n")

        run_cmd(
            [
                "plink2",
                "--pfile",
                first_prefix,
                "--pmerge-list",
                merge_list,
                "--allow-extra-chr",
                "--make-pgen",
                "--out",
                args.merged_prefix,
            ]
        )

    merged_psam = args.merged_prefix + ".psam"
    sample_ids = load_psam_iids(merged_psam)
    sample_map = build_sample_map(sample_ids)

    pheno_file = os.path.join(args.output_dir, "phenotype_for_plink.tsv")
    covar_file = os.path.join(args.output_dir, "covariates_for_plink.tsv")
    prepare_phenotype(args, sample_map, pheno_file)
    pc_cols = prepare_covariates(args, sample_map, covar_file)

    run_cmd(
        [
            "plink2",
            "--pfile",
            args.merged_prefix,
            "--allow-extra-chr",
            "--pheno",
            pheno_file,
            "--pheno-name",
            "PHENO",
            "--covar",
            covar_file,
            "--covar-name",
            ",".join(pc_cols),
            "--glm",
            "hide-covar",
            "--out",
            args.gwas_prefix,
        ]
    )

    glm_file = choose_glm_file(args.gwas_prefix)
    gwas_df = make_gwas_table(glm_file, args.gwas_tsv)
    plot_manhattan(gwas_df, args.manhattan_png, f"{args.ancestry} ancestry dosage GWAS")


if __name__ == "__main__":
    main()

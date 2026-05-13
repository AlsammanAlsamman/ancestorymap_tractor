#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def parse_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def parse_ancestry_map(value):
    mapping = {}
    for part in str(value).split(","):
        idx, label = part.split(":", 1)
        mapping[str(idx)] = str(label)
    return mapping


def add_sample_id_aliases_multi(df, iid_col, value_columns, fid_col=None):
    value_columns = list(value_columns)
    iid = df[iid_col].astype(str)

    base = df[value_columns].copy()
    base.insert(0, "sample_id", iid)
    alias_frames = [base]

    iid_iid = df[value_columns].copy()
    iid_iid.insert(0, "sample_id", iid + "_" + iid)
    alias_frames.append(iid_iid)

    if fid_col is not None and fid_col in df.columns:
        fid = df[fid_col].astype(str)
        fid_iid = df[value_columns].copy()
        fid_iid.insert(0, "sample_id", fid + "_" + iid)
        alias_frames.append(fid_iid)

    alias_df = pd.concat(alias_frames, ignore_index=True)
    alias_df = alias_df.drop_duplicates(subset=["sample_id"], keep="first")
    return alias_df


def load_phenotype(args):
    if args.phenotype_format == "fam":
        fam = pd.read_csv(args.phenotype, sep=r"\s+", header=None, dtype=str, engine="python")
        if fam.shape[1] < 6:
            raise ValueError(f"FAM file must have at least 6 columns; found {fam.shape[1]} in {args.phenotype}")

        phe_df = pd.DataFrame(
            {
                "FID": fam.iloc[:, 0].astype(str),
                "IID": fam.iloc[:, 1].astype(str),
                "SEX": fam.iloc[:, 4].astype(str),
                "PHENO": fam.iloc[:, 5].astype(str),
            }
        )
        phe_df["y"] = phe_df["PHENO"].map(
            {
                str(args.fam_case_code): 1.0,
                str(args.fam_control_code): 0.0,
            }
        )
        phe_df["sex_female"] = phe_df["SEX"].map({"2": 1.0, "1": 0.0})
        keep_cols = ["FID", "IID", "y"]
        if parse_bool(args.include_sex):
            keep_cols.append("sex_female")
        phe_alias = add_sample_id_aliases_multi(phe_df, iid_col="IID", fid_col="FID", value_columns=keep_cols[2:])
        phe_alias["FID"] = phe_alias["sample_id"].str.replace(r"_.*$", "", regex=True)
        phe_alias["IID"] = phe_alias["sample_id"].str.replace(r".*_", "", regex=True)
        phe_alias = phe_alias.drop_duplicates(subset=["sample_id"], keep="first")
        return phe_alias

    phe_df = pd.read_csv(args.phenotype, sep=None, engine="python", dtype=str)
    if args.table_iid_column not in phe_df.columns:
        raise ValueError(
            f"IID column '{args.table_iid_column}' not found in phenotype table {args.phenotype}. "
            f"Available columns: {list(phe_df.columns)}"
        )
    if args.table_phenotype_column not in phe_df.columns:
        raise ValueError(
            f"Phenotype column '{args.table_phenotype_column}' not found in phenotype table {args.phenotype}. "
            f"Available columns: {list(phe_df.columns)}"
        )

    phe_df["y"] = pd.to_numeric(phe_df[args.table_phenotype_column], errors="coerce")
    keep_cols = ["y"]
    if parse_bool(args.include_sex) and "SEX" in phe_df.columns:
        phe_df["sex_female"] = pd.to_numeric(phe_df["SEX"], errors="coerce").map({2.0: 1.0, 1.0: 0.0})
        keep_cols.append("sex_female")
    alias = add_sample_id_aliases_multi(
        phe_df,
        iid_col=args.table_iid_column,
        fid_col="FID" if "FID" in phe_df.columns else None,
        value_columns=keep_cols,
    )
    alias = alias.drop_duplicates(subset=["sample_id"], keep="first")
    alias["IID"] = alias["sample_id"].str.replace(r".*_", "", regex=True)
    return alias


def load_covariates(args):
    if not args.covariates_file:
        return None, []

    has_header = parse_bool(args.covariates_has_header)
    if args.covariates_format != "plink_eigenvec":
        raise ValueError(f"Unsupported covariate format: {args.covariates_format}")

    if has_header:
        cov_df = pd.read_csv(args.covariates_file, sep=r"\s+", dtype=str, engine="python")
    else:
        cov_df = pd.read_csv(args.covariates_file, sep=r"\s+", header=None, dtype=str, engine="python")
        if cov_df.shape[1] < 3:
            raise ValueError(
                f"Covariates file must contain at least FID, IID, and one PC column; found {cov_df.shape[1]} columns in {args.covariates_file}"
            )
        cov_df.columns = [args.covariates_fid_column, args.covariates_iid_column] + [f"PC{i}" for i in range(1, cov_df.shape[1] - 1)]

    pc_columns = [col for col in cov_df.columns if str(col).upper().startswith("PC")][: args.n_pcs]
    if len(pc_columns) < args.n_pcs:
        raise ValueError(
            f"Requested {args.n_pcs} PCs, but only found {len(pc_columns)} in {args.covariates_file}."
        )

    for column in pc_columns:
        cov_df[column] = pd.to_numeric(cov_df[column], errors="coerce")

    alias = add_sample_id_aliases_multi(
        cov_df,
        iid_col=args.covariates_iid_column,
        fid_col=args.covariates_fid_column if args.covariates_fid_column in cov_df.columns else None,
        value_columns=pc_columns,
    )
    alias = alias.drop_duplicates(subset=["sample_id"], keep="first")
    return alias, pc_columns


def load_age(args):
    if not parse_bool(args.include_age) or not args.age_file:
        return None

    age_path = Path(args.age_file)
    if not age_path.exists():
        print(f"Age file was requested but not found: {args.age_file}; age will be skipped.")
        return None

    age_df = pd.read_csv(age_path, sep=None, engine="python", dtype=str)
    if args.age_iid_column not in age_df.columns or args.age_column not in age_df.columns:
        raise ValueError(
            f"Age file must contain columns '{args.age_iid_column}' and '{args.age_column}'. Available: {list(age_df.columns)}"
        )
    age_df["age"] = pd.to_numeric(age_df[args.age_column], errors="coerce")
    alias = add_sample_id_aliases_multi(
        age_df,
        iid_col=args.age_iid_column,
        fid_col="FID" if "FID" in age_df.columns else None,
        value_columns=["age"],
    )
    alias = alias.drop_duplicates(subset=["sample_id"], keep="first")
    return alias


def load_global_ancestry(args, ancestry_map, chromosomes):
    frames = []
    for chr_name in chromosomes:
        path = Path(args.rfmix_q_pattern.replace("{chr}", str(chr_name)))
        q_df = pd.read_csv(path, sep=r"\s+", engine="python")
        q_df.columns = [str(col) for col in q_df.columns]
        sample_col = "#sample" if "#sample" in q_df.columns else q_df.columns[0]
        q_df = q_df.rename(columns={sample_col: "sample_id"})

        keep_cols = ["sample_id"]
        rename_map = {}
        for idx, label in ancestry_map.items():
            idx_col = str(idx)
            if idx_col in q_df.columns:
                rename_map[idx_col] = f"global_{label}_prop"
                keep_cols.append(idx_col)
        q_df = q_df[keep_cols].rename(columns=rename_map)
        frames.append(q_df)

    global_df = pd.concat(frames, ignore_index=True)
    value_cols = [col for col in global_df.columns if col.startswith("global_")]
    for col in value_cols:
        global_df[col] = pd.to_numeric(global_df[col], errors="coerce")
    global_df = global_df.groupby("sample_id", as_index=False)[value_cols].mean()
    return global_df


def select_snps(summary_tsv, ancestry_labels, p_threshold, max_snps_per_ancestry, use_top_snp_per_locus):
    summary = pd.read_csv(summary_tsv, sep="\t")
    summary["chr"] = summary["chr"].astype(str)
    summary["pos"] = pd.to_numeric(summary["pos"], errors="coerce").astype("Int64")

    if "locus_name" not in summary.columns:
        summary["locus_name"] = np.nan
    if "gene" not in summary.columns:
        summary["gene"] = np.nan

    locus_fallback = (
        summary["chr"].astype(str)
        + ":"
        + summary["pos"].astype(str)
        + ":"
        + summary["rsid"].astype(str)
    )
    summary["locus_name"] = summary["locus_name"].astype(str)
    summary["gene"] = summary["gene"].astype(str)
    summary["locus_name"] = summary["locus_name"].where(
        summary["locus_name"].notna() & (summary["locus_name"] != "") & (summary["locus_name"] != "nan"),
        summary["gene"],
    )
    summary["locus_name"] = summary["locus_name"].where(
        summary["locus_name"].notna() & (summary["locus_name"] != "") & (summary["locus_name"] != "nan"),
        locus_fallback,
    )
    selected_frames = []

    for label in ancestry_labels:
        p_col = f"p_dosage_{label}"
        or_col = f"OR_dosage_{label}"
        if p_col not in summary.columns or or_col not in summary.columns:
            raise ValueError(f"Required columns missing from summary table: {p_col}, {or_col}")

        sub = summary[["chr", "pos", "rsid", "locus_name", "gene", p_col, or_col]].copy()
        sub = sub.rename(columns={p_col: "p_value", or_col: "odds_ratio"})
        sub["p_value"] = pd.to_numeric(sub["p_value"], errors="coerce")
        sub["odds_ratio"] = pd.to_numeric(sub["odds_ratio"], errors="coerce")
        sub = sub.replace([np.inf, -np.inf], np.nan).dropna(subset=["pos", "p_value", "odds_ratio"])
        sub = sub[sub["odds_ratio"] > 0].sort_values(["p_value", "pos"])

        selection_mode = "threshold"
        filtered = sub[sub["p_value"] <= p_threshold].copy()
        if filtered.empty:
            filtered = sub.copy()
            selection_mode = "fallback_top_hits"

        if use_top_snp_per_locus and "locus_name" in filtered.columns:
            filtered = filtered.drop_duplicates(subset=["locus_name"], keep="first")
        else:
            filtered = filtered.drop_duplicates(subset=["chr", "pos", "rsid"], keep="first")

        filtered = filtered.head(max_snps_per_ancestry).copy()
        filtered["ancestry"] = label
        filtered["selection_mode"] = selection_mode
        filtered["log_or_weight"] = np.log(filtered["odds_ratio"])
        filtered["variant_key"] = (
            filtered["chr"].astype(str) + ":" + filtered["pos"].astype(str) + ":" + filtered["rsid"].astype(str)
        )
        filtered["found_in_dosage"] = False
        selected_frames.append(filtered)

    selected = pd.concat(selected_frames, ignore_index=True)
    return selected


def load_sample_ids_from_first_dosage(args, ancestry_map, chromosomes):
    first_anc = sorted(ancestry_map.keys(), key=lambda x: int(x))[0]
    first_chr = chromosomes[0]
    path = Path(
        args.tractor_dosage_pattern.replace("{chr}", str(first_chr)).replace("{anc}", str(first_anc))
    )
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
    return header[5:]


def compute_dosage_scores(args, selected, ancestry_map, chromosomes):
    sample_ids = load_sample_ids_from_first_dosage(args, ancestry_map, chromosomes)
    result = pd.DataFrame({"sample_id": sample_ids})
    found_keys = set()

    ancestry_by_label = {label: idx for idx, label in ancestry_map.items()}
    for label in ancestry_by_label:
        feature_name = f"hap_dosage_{label}"
        accum = np.zeros(len(sample_ids), dtype=float)
        label_selected = selected[selected["ancestry"] == label].copy()
        if label_selected.empty:
            result[feature_name] = accum
            continue

        for chr_name in chromosomes:
            chr_selected = label_selected[label_selected["chr"] == str(chr_name)]
            if chr_selected.empty:
                continue

            target_weights = {
                (str(row.chr), str(int(row.pos)), str(row.rsid)): float(row.log_or_weight)
                for row in chr_selected.itertuples(index=False)
            }
            dosage_path = Path(
                args.tractor_dosage_pattern.replace("{chr}", str(chr_name)).replace("{anc}", str(ancestry_by_label[label]))
            )

            with open(dosage_path) as handle:
                header = handle.readline().rstrip("\n").split("\t")
                if header[5:] != sample_ids:
                    raise ValueError(f"Sample order mismatch detected in dosage file: {dosage_path}")
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    key = (parts[0], parts[1], parts[2])
                    weight = target_weights.get(key)
                    if weight is None:
                        continue
                    values = np.asarray(parts[5:], dtype=float)
                    accum += values * weight
                    found_keys.add(f"{key[0]}:{key[1]}:{key[2]}:{label}")

        result[feature_name] = accum

    selected["found_in_dosage"] = selected.apply(
        lambda row: f"{row['chr']}:{int(row['pos'])}:{row['rsid']}:{row['ancestry']}" in found_keys,
        axis=1,
    )
    return result, selected


def build_feature_table(phenotype_df, global_df, covariates_df, age_df, dosage_df):
    feature_df = phenotype_df.merge(global_df, on="sample_id", how="inner")
    feature_df = feature_df.merge(dosage_df, on="sample_id", how="left")

    if covariates_df is not None:
        feature_df = feature_df.merge(covariates_df, on="sample_id", how="left")
    if age_df is not None:
        feature_df = feature_df.merge(age_df, on="sample_id", how="left")

    feature_df = feature_df.drop_duplicates(subset=["sample_id"], keep="first").reset_index(drop=True)
    return feature_df


def evaluate_model(feature_df, model_name, feature_cols, cv_folds, random_state):
    model_df = feature_df[["sample_id", "IID", "y"] + feature_cols].copy()
    for col in feature_cols:
        model_df[col] = pd.to_numeric(model_df[col], errors="coerce")

    model_df = model_df.dropna(subset=["y"]).reset_index(drop=True)
    X = model_df[feature_cols]
    y = model_df["y"].astype(int).to_numpy()

    n_cases = int((y == 1).sum())
    n_controls = int((y == 0).sum())
    n_splits = min(cv_folds, n_cases, n_controls)
    if n_splits < 2:
        raise ValueError(f"Not enough cases/controls for cross-validation in model {model_name}.")

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    probs = np.full(len(model_df), np.nan)
    folds = np.zeros(len(model_df), dtype=int)

    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
        pipeline = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                ("clf", LogisticRegression(max_iter=2000, solver="liblinear")),
            ]
        )
        pipeline.fit(X.iloc[train_idx], y[train_idx])
        probs[test_idx] = pipeline.predict_proba(X.iloc[test_idx])[:, 1]
        folds[test_idx] = fold_idx

    preds = (probs >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, preds, labels=[0, 1]).ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) else np.nan
    specificity = tn / (tn + fp) if (tn + fp) else np.nan

    metrics = pd.DataFrame(
        [
            {
                "model": model_name,
                "n_samples": len(model_df),
                "n_cases": n_cases,
                "n_controls": n_controls,
                "n_features": len(feature_cols),
                "cv_folds": n_splits,
                "auc": roc_auc_score(y, probs),
                "pr_auc": average_precision_score(y, probs),
                "brier_score": brier_score_loss(y, probs),
                "accuracy": accuracy_score(y, preds),
                "balanced_accuracy": balanced_accuracy_score(y, preds),
                "sensitivity": sensitivity,
                "specificity": specificity,
            }
        ]
    )

    predictions = model_df[["sample_id", "IID", "y"]].copy()
    predictions["model"] = model_name
    predictions["fold"] = folds
    predictions["predicted_probability"] = probs
    predictions["predicted_label"] = preds

    full_pipeline = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(max_iter=2000, solver="liblinear")),
        ]
    )
    full_pipeline.fit(X, y)
    clf = full_pipeline.named_steps["clf"]
    coefs = pd.DataFrame(
        {
            "model": model_name,
            "feature": feature_cols,
            "coefficient": clf.coef_[0],
        }
    )
    intercept_row = pd.DataFrame(
        [{"model": model_name, "feature": "intercept", "coefficient": clf.intercept_[0]}]
    )
    coefs = pd.concat([intercept_row, coefs], ignore_index=True)

    return predictions, metrics, coefs


def make_plots(predictions_df, metrics_df, roc_png, calibration_png):
    plt.figure(figsize=(7, 6))
    for model_name, sub in predictions_df.groupby("model"):
        y_true = sub["y"].astype(int).to_numpy()
        y_prob = sub["predicted_probability"].astype(float).to_numpy()
        fpr, tpr, _ = roc_curve(y_true, y_prob)
        auc = metrics_df.loc[metrics_df["model"] == model_name, "auc"].iloc[0]
        plt.plot(fpr, tpr, lw=2, label=f"{model_name} (AUC={auc:.3f})")
    plt.plot([0, 1], [0, 1], linestyle="--", color="gray", lw=1)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title("Admixture-validation ROC curves")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(roc_png, dpi=200)
    plt.close()

    plt.figure(figsize=(7, 6))
    for model_name, sub in predictions_df.groupby("model"):
        y_true = sub["y"].astype(int).to_numpy()
        y_prob = sub["predicted_probability"].astype(float).to_numpy()
        frac_pos, mean_pred = calibration_curve(y_true, y_prob, n_bins=10, strategy="quantile")
        plt.plot(mean_pred, frac_pos, marker="o", lw=1.5, label=model_name)
    plt.plot([0, 1], [0, 1], linestyle="--", color="gray", lw=1)
    plt.xlabel("Mean predicted probability")
    plt.ylabel("Observed case fraction")
    plt.title("Admixture-validation calibration")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(calibration_png, dpi=200)
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser(description="Validate admixture-informed disease prediction using RFMix and Tractor outputs")
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--rfmix-q-pattern", required=True)
    parser.add_argument("--tractor-dosage-pattern", required=True)
    parser.add_argument("--chromosomes", required=True)
    parser.add_argument("--ancestry-map", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--phenotype-format", required=True, choices=["fam", "table"])
    parser.add_argument("--fam-case-code", default="2")
    parser.add_argument("--fam-control-code", default="1")
    parser.add_argument("--table-iid-column", default="IID")
    parser.add_argument("--table-phenotype-column", default="PHENO")
    parser.add_argument("--covariates-file", default="")
    parser.add_argument("--covariates-format", default="plink_eigenvec")
    parser.add_argument("--covariates-has-header", default="false")
    parser.add_argument("--covariates-fid-column", default="FID")
    parser.add_argument("--covariates-iid-column", default="IID")
    parser.add_argument("--n-pcs", type=int, default=5)
    parser.add_argument("--include-sex", default="true")
    parser.add_argument("--include-age", default="false")
    parser.add_argument("--age-file", default="")
    parser.add_argument("--age-iid-column", default="IID")
    parser.add_argument("--age-column", default="AGE")
    parser.add_argument("--global-reference-ancestry", default="AMR")
    parser.add_argument("--snp-pvalue-threshold", type=float, default=0.01)
    parser.add_argument("--max-snps-per-ancestry", type=int, default=50)
    parser.add_argument("--use-top-snp-per-locus", default="true")
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--sample-features-out", required=True)
    parser.add_argument("--selected-snps-out", required=True)
    parser.add_argument("--cv-predictions-out", required=True)
    parser.add_argument("--metrics-out", required=True)
    parser.add_argument("--coefficients-out", required=True)
    parser.add_argument("--roc-png", required=True)
    parser.add_argument("--calibration-png", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    chromosomes = [item for item in str(args.chromosomes).split(",") if item]
    ancestry_map = parse_ancestry_map(args.ancestry_map)
    ancestry_labels = [ancestry_map[idx] for idx in sorted(ancestry_map.keys(), key=lambda x: int(x))]

    phenotype_df = load_phenotype(args)
    covariates_df, pc_columns = load_covariates(args)
    age_df = load_age(args)
    global_df = load_global_ancestry(args, ancestry_map, chromosomes)

    selected = select_snps(
        args.summary_tsv,
        ancestry_labels=ancestry_labels,
        p_threshold=args.snp_pvalue_threshold,
        max_snps_per_ancestry=args.max_snps_per_ancestry,
        use_top_snp_per_locus=parse_bool(args.use_top_snp_per_locus),
    )
    dosage_df, selected = compute_dosage_scores(args, selected, ancestry_map, chromosomes)

    feature_df = build_feature_table(phenotype_df, global_df, covariates_df, age_df, dosage_df)

    global_cols = [f"global_{label}_prop" for label in ancestry_labels if label != args.global_reference_ancestry]
    dosage_cols = [f"hap_dosage_{label}" for label in ancestry_labels]
    optional_cols = []
    if parse_bool(args.include_sex) and "sex_female" in feature_df.columns:
        optional_cols.append("sex_female")
    if parse_bool(args.include_age) and "age" in feature_df.columns and feature_df["age"].notna().any():
        optional_cols.append("age")

    model_definitions = {
        "pcs_only": pc_columns + optional_cols,
        "pcs_plus_global": pc_columns + optional_cols + global_cols,
        "full_admixture": pc_columns + optional_cols + global_cols + dosage_cols,
    }

    prediction_frames = []
    metrics_frames = []
    coef_frames = []

    for model_name, feature_cols in model_definitions.items():
        feature_cols = [col for col in feature_cols if col in feature_df.columns]
        preds, metrics, coefs = evaluate_model(
            feature_df=feature_df,
            model_name=model_name,
            feature_cols=feature_cols,
            cv_folds=args.cv_folds,
            random_state=args.random_state,
        )
        prediction_frames.append(preds)
        metrics_frames.append(metrics)
        coef_frames.append(coefs)
        print(f"Evaluated {model_name} with {len(feature_cols)} features: AUC={metrics['auc'].iloc[0]:.4f}")

    predictions_df = pd.concat(prediction_frames, ignore_index=True)
    metrics_df = pd.concat(metrics_frames, ignore_index=True)
    coefficients_df = pd.concat(coef_frames, ignore_index=True)

    Path(args.sample_features_out).parent.mkdir(parents=True, exist_ok=True)
    feature_df.to_csv(args.sample_features_out, sep="\t", index=False)
    selected.to_csv(args.selected_snps_out, sep="\t", index=False)
    predictions_df.to_csv(args.cv_predictions_out, sep="\t", index=False)
    metrics_df.to_csv(args.metrics_out, sep="\t", index=False)
    coefficients_df.to_csv(args.coefficients_out, sep="\t", index=False)
    make_plots(predictions_df, metrics_df, args.roc_png, args.calibration_png)

    print(
        f"Wrote sample features to {args.sample_features_out}, metrics to {args.metrics_out}, "
        f"and plots to {args.roc_png} / {args.calibration_png}"
    )


if __name__ == "__main__":
    main()

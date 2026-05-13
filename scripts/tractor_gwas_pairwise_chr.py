#!/usr/bin/env python3
import argparse
import hail as hl
import pandas as pd


def add_sample_id_aliases(df, iid_col, outcome_col, fid_col=None):
    iid = df[iid_col].astype(str)
    alias_frames = [pd.DataFrame({"sample_id": iid, outcome_col: df[outcome_col]})]
    alias_frames.append(pd.DataFrame({"sample_id": iid + "_" + iid, outcome_col: df[outcome_col]}))

    if fid_col is not None and fid_col in df.columns:
        fid = df[fid_col].astype(str)
        alias_frames.append(pd.DataFrame({"sample_id": fid + "_" + iid, outcome_col: df[outcome_col]}))

    alias_df = pd.concat(alias_frames, ignore_index=True)
    alias_df = alias_df.drop_duplicates(subset=["sample_id"], keep="first")
    return alias_df


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


def parse_args():
    parser = argparse.ArgumentParser(description="Run per-chromosome pairwise Tractor GWAS with Hail")
    parser.add_argument("--chr", required=True)
    parser.add_argument("--anc-hapcount", required=True)
    parser.add_argument("--anc-i-dosage", required=True)
    parser.add_argument("--anc-j-dosage", required=True)
    parser.add_argument("--label-i", required=True)
    parser.add_argument("--label-j", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--phenotype-format", required=True, choices=["fam", "table"])
    parser.add_argument("--fam-case-code", default="2")
    parser.add_argument("--fam-control-code", default="1")
    parser.add_argument("--table-iid-column", default="IID")
    parser.add_argument("--table-phenotype-column", default="PHENO")
    parser.add_argument("--pca-enabled", default="false")
    parser.add_argument("--pca-file", default="")
    parser.add_argument("--pca-format", default="plink_eigenvec")
    parser.add_argument("--pca-has-header", default="false")
    parser.add_argument("--pca-fid-column", default="FID")
    parser.add_argument("--pca-iid-column", default="IID")
    parser.add_argument("--pca-n-pcs", type=int, default=10)
    parser.add_argument("--min-partitions", type=int, default=32)
    parser.add_argument("--output-tsv", required=True)
    return parser.parse_args()


def import_phenotype(args):
    phe_tmp = args.output_tsv + ".phenotype.tsv"

    if args.phenotype_format == "fam":
        fam = pd.read_csv(args.phenotype, sep=r"\s+", header=None, dtype=str, engine="python")
        if fam.shape[1] < 6:
            raise ValueError(f"FAM file must have at least 6 columns; found {fam.shape[1]} in {args.phenotype}")

        phe_df = pd.DataFrame(
            {
                "FID": fam.iloc[:, 0].astype(str),
                "IID": fam.iloc[:, 1].astype(str),
                "PHENO": fam.iloc[:, 5].astype(str),
            }
        )
        phe_df["y"] = phe_df["PHENO"].map(
            {
                str(args.fam_case_code): 1.0,
                str(args.fam_control_code): 0.0,
            }
        )
        phe_alias_df = add_sample_id_aliases(phe_df, iid_col="IID", fid_col="FID", outcome_col="y")
        phe_alias_df.to_csv(phe_tmp, sep="\t", index=False)
        print(f"Loaded {len(phe_df)} phenotype rows and {len(phe_alias_df)} sample-ID aliases from {args.phenotype}")

        phe = hl.import_table(phe_tmp, types={"sample_id": hl.tstr, "y": hl.tfloat64}).key_by("sample_id")
        return phe

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
    phe_alias_df = add_sample_id_aliases(
        phe_df,
        iid_col=args.table_iid_column,
        fid_col="FID" if "FID" in phe_df.columns else None,
        outcome_col="y",
    )
    phe_alias_df.to_csv(phe_tmp, sep="\t", index=False)
    print(f"Loaded {len(phe_df)} phenotype rows and {len(phe_alias_df)} sample-ID aliases from {args.phenotype}")

    phe = hl.import_table(phe_tmp, types={"sample_id": hl.tstr, "y": hl.tfloat64}).key_by("sample_id")
    return phe


def import_mt(path, min_partitions):
    row_fields = {"CHROM": hl.tstr, "POS": hl.tint, "ID": hl.tstr, "REF": hl.tstr, "ALT": hl.tstr}
    mt = hl.import_matrix_table(path, row_fields=row_fields, row_key=[], min_partitions=min_partitions)
    mt = mt.key_rows_by().drop("row_id")
    mt = mt.key_rows_by(locus=hl.locus(mt.CHROM, mt.POS))
    return mt


def import_pca_covariates(args):
    enabled = str(args.pca_enabled).strip().lower() in {"1", "true", "yes", "y"}
    if not enabled:
        return None, []

    if not args.pca_file:
        raise ValueError("PCA covariates are enabled but no --pca-file was provided.")

    if args.pca_format != "plink_eigenvec":
        raise ValueError(f"Unsupported PCA format: {args.pca_format}. Expected 'plink_eigenvec'.")

    has_header = str(args.pca_has_header).strip().lower() in {"1", "true", "yes", "y"}
    if has_header:
        pca_df = pd.read_csv(args.pca_file, sep=r"\s+", dtype=str, engine="python")
    else:
        pca_df = pd.read_csv(args.pca_file, sep=r"\s+", header=None, dtype=str, engine="python")
        if pca_df.shape[1] < 3:
            raise ValueError(
                f"PCA file must contain at least FID, IID, and one PC column; found {pca_df.shape[1]} columns in {args.pca_file}"
            )
        pca_df.columns = [args.pca_fid_column, args.pca_iid_column] + [f"PC{i}" for i in range(1, pca_df.shape[1] - 1)]

    if args.pca_iid_column not in pca_df.columns:
        raise ValueError(
            f"PCA IID column '{args.pca_iid_column}' was not found in {args.pca_file}. Available columns: {list(pca_df.columns)}"
        )

    pc_columns = [col for col in pca_df.columns if str(col).upper().startswith("PC")]
    if len(pc_columns) < args.pca_n_pcs:
        raise ValueError(
            f"Requested {args.pca_n_pcs} PCs, but only found {len(pc_columns)} in {args.pca_file}."
        )
    pc_columns = pc_columns[: args.pca_n_pcs]

    for column in pc_columns:
        pca_df[column] = pd.to_numeric(pca_df[column], errors="coerce")

    cov_tmp = args.output_tsv + ".pca.tsv"
    cov_alias_df = add_sample_id_aliases_multi(
        pca_df,
        iid_col=args.pca_iid_column,
        fid_col=args.pca_fid_column if args.pca_fid_column in pca_df.columns else None,
        value_columns=pc_columns,
    )
    cov_alias_df = cov_alias_df.dropna(subset=pc_columns, how="any")
    cov_alias_df.to_csv(cov_tmp, sep="\t", index=False)

    cov_types = {"sample_id": hl.tstr}
    cov_types.update({column: hl.tfloat64 for column in pc_columns})
    cov = hl.import_table(cov_tmp, types=cov_types).key_by("sample_id")
    print(f"Loaded PCA covariates from {args.pca_file}: using {len(pc_columns)} PCs ({', '.join(pc_columns)})")
    return cov, pc_columns


def main():
    args = parse_args()
    hl.init()

    anc_hap = import_mt(args.anc_hapcount, args.min_partitions)
    anc_i_dose = import_mt(args.anc_i_dosage, args.min_partitions)
    anc_j_dose = import_mt(args.anc_j_dosage, args.min_partitions)

    phe = import_phenotype(args)
    cov, pc_columns = import_pca_covariates(args)

    total_cols = anc_hap.count_cols()
    mt = anc_hap.annotate_cols(Pheno=phe[anc_hap.col_id])
    if cov is not None:
        mt = mt.annotate_cols(Cov=cov[mt.col_id])

    matched_cols = mt.aggregate_cols(hl.agg.count_where(hl.is_defined(mt.Pheno.y)))
    print(f"Phenotype-matched samples: {matched_cols}/{total_cols}")
    if matched_cols == 0:
        example_ids = [row.sample_id for row in mt.cols().select(sample_id=mt.col_id).take(5)]
        raise ValueError(
            "No overlapping sample IDs were found between the phenotype file and Tractor matrices. "
            f"Example Tractor IDs: {example_ids}"
        )

    if cov is not None:
        matched_cov = mt.aggregate_cols(hl.agg.count_where(hl.is_defined(mt.Cov)))
        print(f"PCA-matched samples: {matched_cov}/{total_cols}")
        if matched_cov == 0:
            example_ids = [row.sample_id for row in mt.cols().select(sample_id=mt.col_id).take(5)]
            raise ValueError(
                "No overlapping sample IDs were found between the PCA covariate file and Tractor matrices. "
                f"Example Tractor IDs: {example_ids}"
            )

    keep_expr = hl.is_defined(mt.Pheno.y)
    if cov is not None:
        keep_expr = keep_expr & hl.is_defined(mt.Cov)
        for column in pc_columns:
            keep_expr = keep_expr & hl.is_defined(mt.Cov[column])

    mt = mt.filter_cols(keep_expr)
    mt = mt.annotate_entries(
        anc_i_dose=anc_i_dose[mt.locus, mt.col_id],
        anc_j_dose=anc_j_dose[mt.locus, mt.col_id],
    )

    regression_terms = [1.0, mt.x, mt.anc_i_dose.x, mt.anc_j_dose.x]
    if cov is not None:
        regression_terms.extend([mt.Cov[column] for column in pc_columns])

    mt = mt.annotate_rows(
        lm=hl.agg.linreg(
            mt.Pheno.y,
            regression_terms,
        )
    )

    mt = mt.select_rows(
        CHROM=mt.locus.contig,
        POS=mt.locus.position,
        ID=mt.ID,
        REF=mt.REF,
        ALT=mt.ALT,
        n=mt.lm.n,
        beta_local_ancestry=mt.lm.beta[1],
        beta_dosage_i=mt.lm.beta[2],
        beta_dosage_j=mt.lm.beta[3],
        se_local_ancestry=mt.lm.standard_error[1],
        se_dosage_i=mt.lm.standard_error[2],
        se_dosage_j=mt.lm.standard_error[3],
        p_local_ancestry=mt.lm.p_value[1],
        p_dosage_i=mt.lm.p_value[2],
        p_dosage_j=mt.lm.p_value[3],
        ancestry_i=hl.literal(args.label_i),
        ancestry_j=hl.literal(args.label_j),
    )

    mt.rows().export(args.output_tsv)
    hl.stop()


if __name__ == "__main__":
    main()

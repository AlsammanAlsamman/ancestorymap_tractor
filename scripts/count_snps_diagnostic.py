#!/usr/bin/env python3
"""
Count SNPs per chromosome per ancestry from per-ancestry PLINK pvar files,
and compare to the original cohort PLINK bim/fam. Outputs a TSV table.
Python 3.7 compatible.
"""

import argparse
import os
import sys


def count_pvar_snps(pvar_path):
    """Count variant lines (non-header) in a .pvar file."""
    count = 0
    with open(pvar_path, "r") as fh:
        for line in fh:
            if not line.startswith("#"):
                count += 1
    return count


def count_psam_samples(psam_path):
    """Count sample lines (non-header) in a .psam file."""
    count = 0
    with open(psam_path, "r") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                count += 1
    return count


def count_bim_snps_by_chr(bim_path):
    """Return dict of {chr_str -> snp_count} from a .bim file."""
    counts = {}
    with open(bim_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom = line.split()[0]
            counts[chrom] = counts.get(chrom, 0) + 1
    return counts


def count_fam_samples(fam_path):
    """Count samples in a .fam file."""
    count = 0
    with open(fam_path, "r") as fh:
        for line in fh:
            if line.strip():
                count += 1
    return count


def main():
    parser = argparse.ArgumentParser(
        description="Count SNPs/samples per ancestry per chromosome vs original PLINK data"
    )
    parser.add_argument("--ancestry-plink-dir", required=True,
                        help="Directory containing per-ancestry PLINK folders")
    parser.add_argument("--ancestries", required=True,
                        help="Comma-separated ancestry labels, e.g. AFR,AMR,EUR,EAS")
    parser.add_argument("--chromosomes", required=True,
                        help="Comma-separated chromosome numbers, e.g. 1,2,...,22")
    parser.add_argument("--original-bim", required=True,
                        help="Original cohort .bim file path")
    parser.add_argument("--original-fam", required=True,
                        help="Original cohort .fam file path")
    parser.add_argument("--output-table", required=True,
                        help="Output TSV file path")
    args = parser.parse_args()

    ancestries = [a.strip() for a in args.ancestries.split(",") if a.strip()]
    chromosomes = [c.strip() for c in args.chromosomes.split(",") if c.strip()]

    # --- Original cohort counts ---
    print("Counting original BIM SNPs per chromosome ...", file=sys.stderr)
    orig_snps_by_chr = count_bim_snps_by_chr(args.original_bim)
    print("Counting original FAM samples ...", file=sys.stderr)
    orig_samples = count_fam_samples(args.original_fam)

    # snp_data[chrom][ancestry] = snp_count
    # sample_data[chrom][ancestry] = sample_count
    snp_data = {}
    sample_data = {}

    for chrom in chromosomes:
        snp_data[chrom] = {}
        sample_data[chrom] = {}
        snp_data[chrom]["original"] = orig_snps_by_chr.get(str(chrom), 0)
        sample_data[chrom]["original"] = orig_samples

    # --- Per-ancestry per-chromosome counts ---
    for ancestry in ancestries:
        print("Counting ancestry: {} ...".format(ancestry), file=sys.stderr)
        for chrom in chromosomes:
            pvar = os.path.join(args.ancestry_plink_dir, ancestry, "chr{}.dosage.pvar".format(chrom))
            psam = os.path.join(args.ancestry_plink_dir, ancestry, "chr{}.dosage.psam".format(chrom))

            if os.path.exists(pvar):
                snp_data[chrom][ancestry] = count_pvar_snps(pvar)
            else:
                snp_data[chrom][ancestry] = 0
                print("WARNING: Missing pvar: {}".format(pvar), file=sys.stderr)

            if os.path.exists(psam):
                sample_data[chrom][ancestry] = count_psam_samples(psam)
            else:
                sample_data[chrom][ancestry] = 0
                print("WARNING: Missing psam: {}".format(psam), file=sys.stderr)

    # --- Write wide-format output ---
    # Columns: chr | original_snps | <ANC>_snps ... | original_samples | <ANC>_samples ...
    out_dir = os.path.dirname(args.output_table)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    col_sources = ["original"] + ancestries  # order: original first, then ancestries sorted

    snp_cols    = ["{}_snps".format(s)    for s in col_sources]
    sample_cols = ["{}_samples".format(s) for s in col_sources]
    header = ["chr"] + snp_cols + sample_cols

    with open(args.output_table, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for chrom in chromosomes:
            snp_vals    = [str(snp_data[chrom].get(s, 0))    for s in col_sources]
            sample_vals = [str(sample_data[chrom].get(s, 0)) for s in col_sources]
            fh.write("\t".join([str(chrom)] + snp_vals + sample_vals) + "\n")

    print("Written {} rows to {}".format(len(rows), args.output_table), file=sys.stderr)


if __name__ == "__main__":
    main()

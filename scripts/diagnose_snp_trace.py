#!/usr/bin/env python3
import glob
import gzip
import os


def find_variant_in_tsv(path: str, chrom: str, pos: str, rsid: str) -> str:
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        _ = handle.readline()
        for raw in handle:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0].replace("chr", "") == chrom and parts[1] == pos and parts[2] == rsid:
                return raw.strip()
    return ""


def find_variant_in_vcfgz(path: str, chrom: str, pos: str) -> str:
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0].replace("chr", "") == chrom and parts[1] == pos:
                return "\t".join(parts[:5])
    return ""


def find_variant_in_gwas(path: str, chrom: str, pos: str, rsid: str) -> str:
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        required = {"CHROM", "POS", "ID"}
        if not required.issubset(set(header)):
            return ""
        i_chr = header.index("CHROM")
        i_pos = header.index("POS")
        i_id = header.index("ID")
        for raw in handle:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) > max(i_chr, i_pos, i_id):
                if parts[i_chr].replace("chr", "") == chrom and parts[i_pos] == pos and parts[i_id] == rsid:
                    return raw.strip()
    return ""


def main():
    chrom = "11"
    pos = "128605169"
    rsid = "rs1236176"

    result_dirs = sorted([d for d in os.listdir(".") if d.startswith("results_") and os.path.isdir(d)])
    print("RESULT_DIR_CANDIDATES", result_dirs)

    for result_dir in result_dirs:
        print(f"\n=== Checking {result_dir} ===")

        maf_path = os.path.join(result_dir, "maf_summary", "merged.annotated.maf.tsv")
        if os.path.exists(maf_path):
            hit = find_variant_in_tsv(maf_path, chrom, pos, rsid)
            if hit:
                print("FOUND in merged.annotated.maf.tsv:", hit)
            else:
                print("NOT found in merged.annotated.maf.tsv")
        else:
            print("missing", maf_path)

        subset_vcf = os.path.join(result_dir, "cohort_region_subset", f"chr{chrom}.cohort_regions.vcf.gz")
        if os.path.exists(subset_vcf):
            hit = find_variant_in_vcfgz(subset_vcf, chrom, pos)
            if hit:
                print("FOUND in cohort_region_subset:", hit)
            else:
                print("NOT found in cohort_region_subset chr11 vcf")
        else:
            print("missing", subset_vcf)

        phased_vcf = os.path.join(result_dir, "cohort_phasing", f"chr{chrom}.phased.vcf.gz")
        if os.path.exists(phased_vcf):
            hit = find_variant_in_vcfgz(phased_vcf, chrom, pos)
            if hit:
                print("FOUND in cohort_phasing:", hit)
            else:
                print("NOT found in cohort_phasing chr11 phased vcf")
        else:
            print("missing", phased_vcf)

        gwas_files = sorted(glob.glob(os.path.join(result_dir, "tractor_gwas", "merged.model_*.gwas.tsv")))
        if gwas_files:
            print("GWAS files:", [os.path.basename(x) for x in gwas_files])
            any_hit = False
            for path in gwas_files:
                hit = find_variant_in_gwas(path, chrom, pos, rsid)
                if hit:
                    print("FOUND in", os.path.basename(path) + ":", hit)
                    any_hit = True
            if not any_hit:
                print("NOT found in any merged.model_*.gwas.tsv")
        else:
            print("No GWAS pair files found")

        final_path = os.path.join(result_dir, "maf_gwas_summary", "merged.maf_gwas_summary.tsv")
        if os.path.exists(final_path):
            hit = find_variant_in_tsv(final_path, chrom, pos, rsid)
            if hit:
                print("FOUND in merged.maf_gwas_summary.tsv:", hit)
            else:
                print("NOT found in merged.maf_gwas_summary.tsv")
        else:
            print("missing", final_path)


if __name__ == "__main__":
    main()

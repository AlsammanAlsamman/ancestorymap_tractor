#!/usr/bin/env python3
import csv
import glob
import gzip
import os
import subprocess
from pathlib import Path

BASE = Path('.')
RESULTS_DIR = BASE / 'results_hispanic_mod'
RSID_FILE = BASE / 'logs/diagnostics/target_missing_rsids.txt'
OUT_TSV = BASE / 'logs/diagnostics/missing_rsids_trace.tsv'
OUT_TXT = BASE / 'logs/diagnostics/missing_rsids_summary.txt'


def read_rsid_list(path: Path):
    with path.open('r', encoding='utf-8') as f:
        return [line.strip() for line in f if line.strip()]


def load_bim(rsids):
    rs_set = set(rsids)
    out = {}
    bim_path = BASE / 'inputs/Hispanic.bim'
    if not bim_path.exists():
        return out
    with bim_path.open('r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            p = line.strip().split()
            if len(p) < 6:
                continue
            rsid = p[1]
            if rsid in rs_set:
                out[rsid] = {
                    'chr': p[0].replace('chr', ''),
                    'pos': p[3],
                    'a1': p[4].upper(),
                    'a2': p[5].upper(),
                }
    return out


def scan_vcf_by_id(vcf_glob, rsids):
    rs_set = set(rsids)
    found = {}
    for fp in sorted(glob.glob(vcf_glob)):
        with gzip.open(fp, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if not line or line[0] == '#':
                    continue
                p = line.rstrip('\n').split('\t')
                if len(p) < 5:
                    continue
                rsid = p[2]
                if rsid in rs_set and rsid not in found:
                    found[rsid] = {
                        'chr': p[0].replace('chr', ''),
                        'pos': p[1],
                        'ref': p[3].upper(),
                        'alt': p[4].upper(),
                    }
    return found


def scan_tsv_first_by_rsid(path, rsids, rsid_col='rsid'):
    out = {}
    if not os.path.exists(path):
        return out
    rs_set = set(rsids)
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        header = f.readline().rstrip('\n').split('\t')
        if rsid_col not in header:
            return out
        i = header.index(rsid_col)
        for line in f:
            p = line.rstrip('\n').split('\t')
            if i < len(p):
                r = p[i]
                if r in rs_set and r not in out:
                    out[r] = p
    return out


def scan_gwas_pairs(rsids):
    rs_set = set(rsids)
    out = {r: [] for r in rsids}
    for fp in sorted(glob.glob(str(RESULTS_DIR / 'tractor_gwas' / 'merged.model_*.gwas.tsv'))):
        pair = os.path.basename(fp).replace('merged.model_', '').replace('.gwas.tsv', '')
        with open(fp, 'r', encoding='utf-8', errors='ignore') as f:
            header = f.readline().rstrip('\n').split('\t')
            if 'ID' not in header:
                continue
            i = header.index('ID')
            for line in f:
                p = line.rstrip('\n').split('\t')
                if i < len(p):
                    r = p[i]
                    if r in rs_set:
                        out[r].append(pair)
    return out


def bcftools_ref_at_pos(chrom, pos):
    if not chrom or not pos:
        return ('', '')
    cand = sorted(glob.glob(str(RESULTS_DIR / 'reference_panel_subset' / f'chr{chrom}.*.vcf.gz')))
    if not cand:
        return ('', '')
    fp = cand[0]
    cmd = ['bcftools', 'query', '-f', '%REF\t%ALT\n', '-r', f'{chrom}:{pos}-{pos}', fp]
    p = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if p.returncode != 0 or not p.stdout.strip():
        return ('', '')
    first = p.stdout.splitlines()[0].split('\t')
    if len(first) < 2:
        return ('', '')
    return (first[0].upper(), first[1].upper())


def infer_reason(row):
    if row['in_bim'] == '0':
        return 'absent_in_inputs_Hispanic.bim'
    if row['in_step3_raw'] == '0':
        return 'filtered_at_step3_loci_extraction_or_id_changed'
    if row['in_step3c_harmonized'] == '0':
        if row['refpanel_ref_alt']:
            return 'filtered_at_step3c_allele_mismatch_or_not_in_ref'
        return 'filtered_at_step3c_not_in_reference_subset'
    if row['in_step4_phased'] == '0':
        return 'filtered_at_step4_phasing_not_in_reference_or_qc'
    if row['in_merged_annotated_maf'] == '0':
        return 'filtered_before_or_during_maf_merge'
    if row['n_gwas_pairs'] == '0':
        return 'absent_in_tractor_gwas_outputs'
    if row['in_merged_maf_gwas_summary'] == '0':
        return 'filtered_at_step9_requires_complete_metrics'
    if row['in_step10_filtered'] == '0':
        return 'filtered_at_step10_thresholds'
    return 'present_in_final_filtered_set'


def main():
    rsids = read_rsid_list(RSID_FILE)

    bim = load_bim(rsids)
    step3 = scan_vcf_by_id(str(RESULTS_DIR / 'cohort_region_subset' / 'chr*.loci.vcf.gz'), rsids)
    step3c = scan_vcf_by_id(str(RESULTS_DIR / 'cohort_region_subset_harmonized' / 'chr*.loci.harmonized.vcf.gz'), rsids)
    step4 = scan_vcf_by_id(str(RESULTS_DIR / 'cohort_phasing' / 'chr*.phased.vcf.gz'), rsids)

    cohort_maf = scan_tsv_first_by_rsid(str(RESULTS_DIR / 'maf_summary' / 'cohort.maf.tsv'), rsids, 'rsid')
    merged_maf = scan_tsv_first_by_rsid(str(RESULTS_DIR / 'maf_summary' / 'merged.annotated.maf.tsv'), rsids, 'rsid')
    merged_summary = scan_tsv_first_by_rsid(str(RESULTS_DIR / 'maf_gwas_summary' / 'merged.maf_gwas_summary.tsv'), rsids, 'rsid')
    step10_filtered = scan_tsv_first_by_rsid(str(RESULTS_DIR / 'locus_ancestry_report' / 'filtered_snps_for_locus_report.tsv'), rsids, 'rsid')

    gwas_pairs = scan_gwas_pairs(rsids)

    rows = []
    for rsid in rsids:
        b = bim.get(rsid, {})
        r3 = step3.get(rsid, {})
        r3c = step3c.get(rsid, {})
        r4 = step4.get(rsid, {})

        chrom = r3.get('chr') or r3c.get('chr') or r4.get('chr') or b.get('chr', '')
        pos = r3.get('pos') or r3c.get('pos') or r4.get('pos') or b.get('pos', '')
        ref_ref, ref_alt = bcftools_ref_at_pos(chrom, pos)

        row = {
            'rsid': rsid,
            'chr': chrom,
            'pos': pos,
            'in_bim': '1' if rsid in bim else '0',
            'in_step3_raw': '1' if rsid in step3 else '0',
            'in_step3c_harmonized': '1' if rsid in step3c else '0',
            'in_step4_phased': '1' if rsid in step4 else '0',
            'in_cohort_maf': '1' if rsid in cohort_maf else '0',
            'in_merged_annotated_maf': '1' if rsid in merged_maf else '0',
            'n_gwas_pairs': str(len(set(gwas_pairs.get(rsid, [])))),
            'gwas_pairs': ','.join(sorted(set(gwas_pairs.get(rsid, [])))),
            'in_merged_maf_gwas_summary': '1' if rsid in merged_summary else '0',
            'in_step10_filtered': '1' if rsid in step10_filtered else '0',
            'raw_ref_alt': (r3.get('ref', '') + '/' + r3.get('alt', '')).strip('/'),
            'harm_ref_alt': (r3c.get('ref', '') + '/' + r3c.get('alt', '')).strip('/'),
            'phased_ref_alt': (r4.get('ref', '') + '/' + r4.get('alt', '')).strip('/'),
            'refpanel_ref_alt': (ref_ref + '/' + ref_alt).strip('/'),
            'bim_a1_a2': ((b.get('a1', '') + '/' + b.get('a2', '')).strip('/')),
        }
        row['likely_reason'] = infer_reason(row)
        rows.append(row)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        'rsid', 'chr', 'pos', 'in_bim', 'in_step3_raw', 'in_step3c_harmonized', 'in_step4_phased',
        'in_cohort_maf', 'in_merged_annotated_maf', 'n_gwas_pairs', 'gwas_pairs',
        'in_merged_maf_gwas_summary', 'in_step10_filtered', 'raw_ref_alt', 'harm_ref_alt',
        'phased_ref_alt', 'refpanel_ref_alt', 'bim_a1_a2', 'likely_reason'
    ]
    with OUT_TSV.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    reason_counts = {}
    for r in rows:
        reason_counts[r['likely_reason']] = reason_counts.get(r['likely_reason'], 0) + 1

    with OUT_TXT.open('w', encoding='utf-8') as f:
        f.write(f'total_rsids={len(rows)}\n')
        for k in sorted(reason_counts):
            f.write(f'{k}={reason_counts[k]}\n')

    print(f'Wrote {OUT_TSV}')
    print(f'Wrote {OUT_TXT}')


if __name__ == '__main__':
    main()

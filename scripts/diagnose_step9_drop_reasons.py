#!/usr/bin/env python3
import csv
from pathlib import Path

BASE = Path('.')
RESULTS = BASE / 'results_hispanic_mod'
TRACE = BASE / 'logs/diagnostics/missing_rsids_trace.tsv'
OUT = BASE / 'logs/diagnostics/missing_rsids_step9_detail.tsv'


def read_trace_ids():
    ids = []
    with TRACE.open('r', encoding='utf-8') as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            if row.get('likely_reason') == 'filtered_at_step9_requires_complete_metrics':
                ids.append(row['rsid'])
    return ids


def load_maf_rows(target_ids):
    target = set(target_ids)
    maf = {}
    path = RESULTS / 'maf_summary' / 'merged.annotated.maf.tsv'
    with path.open('r', encoding='utf-8', errors='ignore') as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            rsid = row.get('rsid', '')
            if rsid in target:
                maf[rsid] = {
                    'chr': str(row.get('chr', '')).replace('chr', ''),
                    'pos': str(row.get('pos', '')),
                    'ref': str(row.get('ref', '')).upper(),
                    'alt': str(row.get('alt', '')).upper(),
                }
    return maf


def inspect_pair(pair_file, pair_name, target_map):
    out = {k: {'present': 0, 'same': 0, 'flip': 0, 'mismatch': 0} for k in target_map}
    with pair_file.open('r', encoding='utf-8', errors='ignore') as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            rsid = row.get('ID', '')
            if rsid not in target_map:
                continue
            out[rsid]['present'] = 1
            maf = target_map[rsid]
            g_ref = str(row.get('REF', '')).upper()
            g_alt = str(row.get('ALT', '')).upper()
            if g_ref == maf['ref'] and g_alt == maf['alt']:
                out[rsid]['same'] = 1
            elif g_ref == maf['alt'] and g_alt == maf['ref']:
                out[rsid]['flip'] = 1
            else:
                out[rsid]['mismatch'] = 1
    return out


def main():
    target_ids = read_trace_ids()
    maf = load_maf_rows(target_ids)

    pair_files = sorted((RESULTS / 'tractor_gwas').glob('merged.model_*.gwas.tsv'))
    pair_names = [p.name.replace('merged.model_', '').replace('.gwas.tsv', '') for p in pair_files]

    summary = {rsid: {'pairs_total': len(pair_files), 'pairs_present': 0, 'pairs_oriented': 0, 'pairs_mismatch': 0} for rsid in maf}

    for pf, pn in zip(pair_files, pair_names):
        stats = inspect_pair(pf, pn, maf)
        for rsid, s in stats.items():
            if s['present']:
                summary[rsid]['pairs_present'] += 1
            if s['same'] or s['flip']:
                summary[rsid]['pairs_oriented'] += 1
            if s['mismatch']:
                summary[rsid]['pairs_mismatch'] += 1

    with OUT.open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['rsid', 'pairs_total', 'pairs_present', 'pairs_oriented', 'pairs_mismatch', 'likely_step9_drop_reason'])
        for rsid in sorted(summary):
            s = summary[rsid]
            if s['pairs_oriented'] < s['pairs_total']:
                reason = 'missing_metrics_due_to_orientation_or_absent_pair'
            else:
                reason = 'likely_missing_other_required_metrics'
            w.writerow([rsid, s['pairs_total'], s['pairs_present'], s['pairs_oriented'], s['pairs_mismatch'], reason])

    print(f'Wrote {OUT}')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
import csv
from pathlib import Path

BASE = Path('.')
RESULTS = BASE / 'results_hispanic_mod'
TRACE = BASE / 'logs/diagnostics/missing_rsids_trace.tsv'
OUT_TSV = BASE / 'logs/diagnostics/step9_13_rootcause.tsv'
OUT_SUMMARY = BASE / 'logs/diagnostics/step9_13_rootcause_summary.txt'

PAIR_NAMES = ['AFR_AMR', 'AFR_EAS', 'AFR_EUR', 'AMR_EAS', 'AMR_EUR', 'EUR_EAS']
ANCESTRY_LABELS = ['AFR', 'AMR', 'EUR', 'EAS']


def read_targets():
    targets = []
    with TRACE.open('r', encoding='utf-8') as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            if row.get('likely_reason') == 'filtered_at_step9_requires_complete_metrics':
                targets.append((row['rsid'], row['chr']))
    return targets


def dosage_stats_for_rsid(chr_name, anc_idx, rsid):
    path = RESULTS / 'tractor' / f'chr{chr_name}.phased.anc{anc_idx}.dosage.txt'
    if not path.exists():
        return None

    with path.open('r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            p = line.strip().split()
            if len(p) < 6:
                continue
            if p[2] != rsid:
                continue

            vals = []
            for x in p[5:]:
                try:
                    vals.append(float(x))
                except ValueError:
                    continue

            n = len(vals)
            if n == 0:
                return {
                    'n': 0,
                    'min': '',
                    'max': '',
                    'mean': '',
                    'var': '',
                    'nonzero_count': 0,
                    'all_zero': 1,
                }

            s = sum(vals)
            ss = sum(v * v for v in vals)
            mean = s / n
            var = (ss - (s * s) / n) / (n - 1) if n > 1 else 0.0
            nonzero = sum(1 for v in vals if v != 0.0)
            all_zero = 1 if nonzero == 0 else 0

            return {
                'n': n,
                'min': min(vals),
                'max': max(vals),
                'mean': mean,
                'var': var,
                'nonzero_count': nonzero,
                'all_zero': all_zero,
            }
    return None


def pair_metric_status(chr_name, pair_name, rsid):
    path = RESULTS / 'tractor_gwas' / f'chr{chr_name}.model_{pair_name}.gwas.tsv'
    if not path.exists():
        return {'found': 0, 'all_na_metrics': 1, 'na_count': -1, 'metric_count': -1}

    with path.open('r', encoding='utf-8', errors='ignore') as f:
        r = csv.DictReader(f, delimiter='\t')
        metric_cols = [
            'beta_local_ancestry', 'beta_dosage_i', 'beta_dosage_j',
            'se_local_ancestry', 'se_dosage_i', 'se_dosage_j',
            'p_local_ancestry', 'p_dosage_i', 'p_dosage_j'
        ]
        for row in r:
            if row.get('ID') != rsid:
                continue
            vals = [row.get(c, '') for c in metric_cols]
            na_count = sum(1 for v in vals if str(v).strip().upper() in {'', 'NA', 'NAN'})
            all_na = 1 if na_count == len(metric_cols) else 0
            return {'found': 1, 'all_na_metrics': all_na, 'na_count': na_count, 'metric_count': len(metric_cols)}

    return {'found': 0, 'all_na_metrics': 1, 'na_count': -1, 'metric_count': -1}


def classify(zero_var_ancestries, pair_all_na_count, pair_found_count):
    if pair_found_count == 0:
        return 'pair_rows_not_found'
    if pair_all_na_count == 6 and len(zero_var_ancestries) > 0:
        return 'all_pairs_na_with_zero_variance_dosage_predictors'
    if pair_all_na_count == 6 and len(zero_var_ancestries) == 0:
        return 'all_pairs_na_without_zero_variance_predictors'
    if 0 < pair_all_na_count < 6:
        return 'mixed_pairs_some_na_some_estimated'
    return 'not_step9_na_pattern'


def main():
    targets = read_targets()
    rows = []

    for rsid, chr_name in targets:
        ancestry_stats = {}
        zero_var_ancestries = []

        for idx, label in enumerate(ANCESTRY_LABELS):
            st = dosage_stats_for_rsid(chr_name, idx, rsid)
            ancestry_stats[label] = st
            if st is not None:
                var = st['var']
                if st['all_zero'] == 1 or (var != '' and abs(float(var)) < 1e-12):
                    zero_var_ancestries.append(label)

        pair_states = {}
        pair_all_na_count = 0
        pair_found_count = 0
        for pair in PAIR_NAMES:
            ps = pair_metric_status(chr_name, pair, rsid)
            pair_states[pair] = ps
            if ps['found'] == 1:
                pair_found_count += 1
            if ps['all_na_metrics'] == 1:
                pair_all_na_count += 1

        cause = classify(zero_var_ancestries, pair_all_na_count, pair_found_count)

        row = {
            'rsid': rsid,
            'chr': chr_name,
            'zero_var_ancestries': ','.join(zero_var_ancestries),
            'n_zero_var_ancestries': len(zero_var_ancestries),
            'pair_rows_found': pair_found_count,
            'pairs_all_na_metrics': pair_all_na_count,
            'root_cause': cause,
        }

        for label in ANCESTRY_LABELS:
            st = ancestry_stats[label]
            row[f'{label}_var'] = '' if st is None else st['var']
            row[f'{label}_nonzero_count'] = '' if st is None else st['nonzero_count']

        for pair in PAIR_NAMES:
            row[f'{pair}_all_na'] = pair_states[pair]['all_na_metrics']

        rows.append(row)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        'rsid', 'chr', 'zero_var_ancestries', 'n_zero_var_ancestries',
        'pair_rows_found', 'pairs_all_na_metrics', 'root_cause',
        'AFR_var', 'AMR_var', 'EUR_var', 'EAS_var',
        'AFR_nonzero_count', 'AMR_nonzero_count', 'EUR_nonzero_count', 'EAS_nonzero_count',
        'AFR_AMR_all_na', 'AFR_EAS_all_na', 'AFR_EUR_all_na',
        'AMR_EAS_all_na', 'AMR_EUR_all_na', 'EUR_EAS_all_na',
    ]

    with OUT_TSV.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    counts = {}
    for r in rows:
        counts[r['root_cause']] = counts.get(r['root_cause'], 0) + 1

    with OUT_SUMMARY.open('w', encoding='utf-8') as f:
        f.write(f'total_step9_rsids={len(rows)}\n')
        for k in sorted(counts):
            f.write(f'{k}={counts[k]}\n')

    print(f'Wrote {OUT_TSV}')
    print(f'Wrote {OUT_SUMMARY}')


if __name__ == '__main__':
    main()

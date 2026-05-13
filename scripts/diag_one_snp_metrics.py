#!/usr/bin/env python3
import csv
import glob

rsid = 'rs11221450'
cols = [
    'CHROM', 'POS', 'ID', 'REF', 'ALT',
    'p_local_ancestry', 'p_dosage_i', 'p_dosage_j',
    'beta_dosage_i', 'beta_dosage_j', 'ancestry_i', 'ancestry_j'
]

for fp in sorted(glob.glob('results_hispanic_mod/tractor_gwas/merged.model_*.gwas.tsv')):
    pair = fp.split('merged.model_')[1].split('.gwas.tsv')[0]
    found = False
    with open(fp, 'r', encoding='utf-8', errors='ignore') as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            if row.get('ID') == rsid:
                found = True
                vals = [row.get(k, '') for k in cols]
                print(pair + '\t' + '\t'.join(vals))
                break
    if not found:
        print(pair + '\tNOT_FOUND')

import sys
from pathlib import Path
import numpy as np, pandas as pd
B=Path('/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2/validation/2026-07-30_0.5.0-dev_search_redesign')
sys.path.insert(0,str(B/'code/mhc_bench')); import metrics as M
A=Path('/tmp/claude-847728/-mnt-vast-standard-home-a-hajialiasgarynaj01-u14286/d416d4b8-f030-4400-8270-87c634df447a/scratchpad/assign')
man=pd.read_csv(B/'runs/mhc_bench_sep_kmer_anchor/pool_manifest.csv').set_index('pool')

rows=[]
for f in sorted(A.glob('*.tsv')):
    pool=f.stem; d=pd.read_csv(f,sep='\t')
    al=d.allele.to_numpy(); cl=d.cluster_id.to_numpy()
    cur=M.evaluate(al,cl)
    # ORACLE MERGE: relabel each fragment by its majority allele, merge same-label fragments
    maj=d.groupby('cluster_id').allele.agg(lambda s:s.value_counts().idxmax())
    orc=M.evaluate(al,d.cluster_id.map(maj).to_numpy())
    rows.append(dict(pool=pool,alleles=man.loc[pool,'allele_count'],peptides=len(d),
        cur_ami=cur['ami'],cur_pur=cur['adjusted_purity_macro'],cur_rec=cur['bcubed_recall_macro'],
        cur_f1=cur['bcubed_f1_macro'],cur_k=cur['clusters'],
        orc_ami=orc['ami'],orc_pur=orc['adjusted_purity_macro'],orc_rec=orc['bcubed_recall_macro'],
        orc_f1=orc['bcubed_f1_macro'],orc_k=orc['clusters']))
r=pd.DataFrame(rows)
r['band']=pd.cut(r.alleles,bins=[1,6,12,30],labels=['2-6','7-12','13-30'])
r.to_csv('ceiling.csv',index=False)

print('=== PepCluster2 as-is  vs  ORACLE-MERGED (upper bound on any merge step) ===\n')
cols=['cur_k','orc_k','cur_ami','orc_ami','cur_pur','orc_pur','cur_rec','orc_rec','cur_f1','orc_f1']
print('OVERALL (48 test pools)');print(r[cols].mean().round(4).to_string());print()
print('BY BAND');print(r.groupby('band',observed=True)[cols].mean().round(4).to_string());print()
print('alleles per band:');print(r.groupby('band',observed=True).alleles.mean().round(1).to_string())
print()
print('MixMHCp reference (test): default AMI .392 / oracle_k AMI .491')
print('by band AMI  default: 2-6 .573  7-12 .441  13-30 .253')
print('by band AMI oracle_k: 2-6 .712  7-12 .535  13-30 .335')

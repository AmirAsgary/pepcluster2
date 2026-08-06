import pandas as pd, numpy as np
B='/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2/validation/2026-07-30_0.5.0-dev_search_redesign'
M=['ami','adjusted_purity_macro','bcubed_precision_macro','bcubed_recall_macro','bcubed_f1_macro','clusters']

def ours(tag,name):
    d=pd.read_csv(f'{B}/runs/mhc_bench_{tag}/tables/test_selected_runs.csv')
    for meth in d['method'].unique():
        s=d[d.method==meth].copy(); s['tool']=f'{name} [{meth}]'
        yield s

frames=[]
for tag,name in [('sep_kmer_anchor','PC2 k-mer'),('sep_aln_anchor','PC2 aln')]:
    frames+=list(ours(tag,name))

mm=pd.read_csv(f'{B}/benchmark/results/immuneapp/raw/mixmhcp.csv')
mm=mm[(mm.split=='test')&(mm.status=='ok')].copy()
mm['tool']='MixMHCp ('+mm.setting+')'
frames.append(mm)

d=pd.concat(frames,ignore_index=True)
print('=== TEST SPLIT, all metrics ===')
t=d.groupby('tool')[M].mean().round(4)
print(t.to_string())
print()
d['band']=pd.cut(d.allele_count,bins=[1,6,12,30],labels=['2-6','7-12','13-30'])
for m in ['ami','adjusted_purity_macro','bcubed_recall_macro','bcubed_f1_macro','clusters']:
    print(f'--- {m} by band ---')
    print(d.pivot_table(index='band',columns='tool',values=m,observed=True).round(4).to_string())
    print()

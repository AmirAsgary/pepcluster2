import pandas as pd, numpy as np
SP='/tmp/claude-847728/-mnt-vast-standard-home-a-hajialiasgarynaj01-u14286/d416d4b8-f030-4400-8270-87c634df447a/scratchpad/'
d=pd.read_csv(SP+'full.csv')
d['band']=pd.cut(d.alleles,bins=[1,6,12,30],labels=['2-6','7-12','13-30'])
d['cfg']=np.where(d['mode']=='bhc','BHC a='+d.alpha.map(lambda x:f'{x:.0e}'),
                  'BF a0='+d.a0.astype(str)+' t='+d.t.astype(int).astype(str))
print('=== BASELINES (48 test pools) ===')
b=d.drop_duplicates('pool')
print(f"PepCluster2 as-is   AMI {b.base_ami.mean():.4f}  F1 {b.base_f1.mean():.4f}  k {b.base_k.mean():.1f}")
print(f"Oracle merge        AMI {b.orc_ami.mean():.4f}  F1 {b.orc_f1.mean():.4f}  k {b.orc_k.mean():.1f}")
print( "MixMHCp default     AMI 0.3918  F1 0.3918  k 4.1")
print( "MixMHCp forced k    AMI 0.4915  F1 0.4728  k 12.4")
print('\n=== ALL CONFIGS, overall ===')
g=d.groupby('cfg').agg(mg_k=('mg_k','mean'),mg_ami=('mg_ami','mean'),mg_f1=('mg_f1','mean'),
   em_k=('em_k','mean'),em_ami=('em_ami','mean'),em_f1=('em_f1','mean'),
   em_pur=('em_pur','mean'),em_rec=('em_rec','mean')).round(4)
print(g.sort_values('em_ami',ascending=False).to_string())
best=g.em_ami.idxmax()
print(f'\nBEST BY EM AMI: {best}')
print('\n=== BEST CONFIG BY BAND ===')
s=d[d.cfg==best]
print(s.groupby('band',observed=True).agg(pools=('pool','size'),alleles=('alleles','mean'),
  base_ami=('base_ami','mean'),mg_k=('mg_k','mean'),mg_ami=('mg_ami','mean'),
  em_k=('em_k','mean'),em_ami=('em_ami','mean'),em_f1=('em_f1','mean'),
  orc_ami=('orc_ami','mean')).round(3).to_string())
print('\n=== EM vs ORACLE, best config ===')
print(f"pools where EM beats merge-only oracle: {(s.em_ami>s.orc_ami).sum()}/{len(s)}")
print('\n=== BHC (negative result) ===')
print(d[d['mode']=='bhc'].groupby('cfg')[['mg_k','mg_ami','em_k','em_ami']].mean().round(4).to_string())

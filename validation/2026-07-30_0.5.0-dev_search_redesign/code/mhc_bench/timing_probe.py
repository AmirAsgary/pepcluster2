import sys, time
from pathlib import Path
sys.path.insert(0,'code/mhc_bench')
import run_grid as R, pandas as pd
root=Path('runs/mhc_bench'); binary=Path('../../target/release/pepcluster2').resolve()
m=pd.read_csv(root/'pool_manifest.csv'); m=m[m.split=='inner']
# one pool near each size decade
picks=[]
for lo,hi in ((1000,1200),(3000,4000),(9000,12000),(28000,40000),(90000,100000)):
    d=m[(m.peptides>=lo)&(m.peptides<=hi)]
    if len(d): picks.append(d.iloc[0])
print('%-28s %8s %10s %10s'%('pool','peptides','graph_s','lazy_s'))
tot={}
for p in picks:
    row={'pool':p.pool,'split':'inner','outer_fold':p.outer_fold,'allele_count':p.allele_count,
         'peptides':p.peptides,'representative_order':'coverage',
         'alignment_threshold':0.50,'anchor_threshold':0.60}
    times={}
    for meth in ('graph','greedy_lazy'):
        t=time.time(); r=R.run_job({**row,'method':meth}, binary, root/'pools', 8)
        times[meth]=time.time()-t
        if r['status']!='ok': print('FAIL',meth,r.get('error'))
    print('%-28s %8d %10.1f %10.1f'%(p.pool,p.peptides,times['graph'],times['greedy_lazy']))
    tot[p.peptides]=times

import sys, time
from pathlib import Path
sys.path.insert(0,'code/mhc_bench')
import run_grid as R, pandas as pd
root=Path('runs/mhc_bench'); binary=Path('../../target/release/pepcluster2').resolve()
tmp=Path('/cbscratch/amirasgary2/pc2_mhc_tmp/probe'); tmp.mkdir(parents=True,exist_ok=True)
m=pd.read_csv(root/'pool_manifest.csv'); m=m[m.split=='inner']
print('%-9s %-12s %7s %8s %9s %11s'%('peptides','method','threads','wall_s','core_s','efficiency'))
for lo,hi in ((2000,6000),(28000,40000)):
    p=m[(m.peptides>=lo)&(m.peptides<=hi)].iloc[0]
    for meth in ('graph','greedy_lazy'):
        base=None
        for th in (1,2,4,8,16):
            row={'pool':p.pool,'split':'inner','outer_fold':p.outer_fold,'allele_count':p.allele_count,
                 'peptides':p.peptides,'method':meth,'representative_order':'coverage',
                 'alignment_threshold':0.45,'anchor_threshold':0.60}
            t=time.time(); r=R.run_job(row,binary,root/'pools',th,tmp); w=time.time()-t
            if base is None: base=w
            print('%-9d %-12s %7d %8.2f %9.1f %10.0f%%'%(p.peptides,meth,th,w,w*th,100*base/(w*th)))

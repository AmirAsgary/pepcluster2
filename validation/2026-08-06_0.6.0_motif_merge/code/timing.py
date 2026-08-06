import time, numpy as np, pandas as pd
from pathlib import Path
from scipy.special import gammaln
exec(open('merge.py').read().split("rows=[]")[0].replace("import sys, time","import sys, time"))
A=Path('.')/'assign'
# largest pool
sizes={f:sum(1 for _ in open(f)) for f in A.glob('*.tsv')}
f=max(sizes,key=sizes.get); print(f.stem,'peptides',sizes[f]-1)
d=pd.read_csv(f,sep='\t')
t=time.time(); ids,C=profiles(d); t_prof=time.time()-t
print(f'clusters={len(ids)}  profile build={t_prof:.3f}s')
bg=C.sum((0,1)); bg/=bg.sum(); lml=make_lml(10.0*bg*20+1e-9)
K=len(C)
t=time.time()
single=lml(C); BF=np.full((K,K),-np.inf)
for i in range(K):
    j=np.arange(i+1,K)
    if len(j): BF[i,j]=lml(C[i][None]+C[j])-single[i]-single[j]
print(f'initial K^2 BF matrix={time.time()-t:.3f}s  ({K*(K-1)//2} pairs)')
t=time.time(); order=agglomerate(C,lml); print(f'full agglomeration to k=1={time.time()-t:.3f}s')

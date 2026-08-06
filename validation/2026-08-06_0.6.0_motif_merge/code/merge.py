"""Prototype of the proposed Bayesian profile-merge stage.

Frame: 9 columns = first 4 + last 5 residues (standard MHC-I frame, handles
length variation, and matches where PepCluster2 already puts its weight).
Criterion: Dirichlet-multinomial marginal likelihood ratio, one shared profile
vs two separate profiles, summed over columns (positions independent).
Agglomerate greedily by log Bayes factor; record metrics at every k.
"""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.special import gammaln
B=Path('/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2/validation/2026-07-30_0.5.0-dev_search_redesign')
sys.path.insert(0,str(B/'code/mhc_bench')); import metrics as M
A=Path('/tmp/claude-847728/-mnt-vast-standard-home-a-hajialiasgarynaj01-u14286/d416d4b8-f030-4400-8270-87c634df447a/scratchpad/assign')
man=pd.read_csv(B/'runs/mhc_bench_sep_kmer_anchor/pool_manifest.csv').set_index('pool')
AA='ACDEFGHIKLMNPQRSTVWY'; IDX={a:i for i,a in enumerate(AA)}; NP=9

def frame(p):
    n=len(p)
    if n<NP: return None
    return [p[0],p[1],p[2],p[3],p[n-5],p[n-4],p[n-3],p[n-2],p[n-1]]

def profiles(d):
    ids=sorted(d.cluster_id.unique()); pos={c:i for i,c in enumerate(ids)}
    C=np.zeros((len(ids),NP,20))
    for pep,cid in zip(d.peptide.values,d.cluster_id.values):
        f=frame(pep)
        if f is None: continue
        r=pos[cid]
        for j,a in enumerate(f):
            if a in IDX: C[r,j,IDX[a]]+=1
    return ids,C

def make_lml(alpha):
    A0=alpha.sum(); ga=gammaln(alpha).sum()
    def lml(x):                       # x: (...,NP,20)
        n=x.sum(-1)
        return (gammaln(A0)-gammaln(A0+n)+gammaln(x+alpha).sum(-1)-ga).sum(-1)
    return lml

def agglomerate(C,lml):
    K=len(C); act=np.arange(K); cur=C.copy(); single=lml(cur)
    lab=np.arange(K); order=[]
    # pairwise logBF
    BF=np.full((K,K),-np.inf)
    for i in range(K):
        j=np.arange(i+1,K)
        if len(j): BF[i,j]=lml(cur[i][None]+cur[j])-single[i]-single[j]
    alive=np.ones(K,bool)
    while alive.sum()>1:
        m=BF.copy(); m[~alive,:]=-np.inf; m[:,~alive]=-np.inf
        i,j=np.unravel_index(np.argmax(m),m.shape); v=m[i,j]
        order.append((i,j,v))
        cur[i]=cur[i]+cur[j]; single[i]=lml(cur[i]); alive[j]=False
        BF[j,:]=-np.inf; BF[:,j]=-np.inf
        oth=np.where(alive)[0]; oth=oth[oth!=i]
        if len(oth):
            b=lml(cur[i][None]+cur[oth])-single[i]-single[oth]
            BF[np.minimum(i,oth),np.maximum(i,oth)]=b
    return order

def evaluate_path(d,ids,order,al):
    pos={c:i for i,c in enumerate(ids)}
    lab=np.array([pos[c] for c in d.cluster_id.values])
    parent=np.arange(len(ids))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    out=[]
    K=len(ids)
    for step,(i,j,v) in enumerate(order):
        parent[find(j)]=find(i)
        K-=1
        cl=np.array([find(x) for x in lab])
        out.append((K,v,M.evaluate(al,cl)))
    return out

rows=[]; t0=time.time()
for f in sorted(A.glob('*.tsv')):
    pool=f.stem; d=pd.read_csv(f,sep='\t'); al=d.allele.to_numpy()
    ids,C=profiles(d)
    bg=C.sum((0,1)); bg=bg/bg.sum()
    for a0 in (1.0,10.0):
        lml=make_lml(a0*bg*20+1e-9)
        order=agglomerate(C,lml)
        path=evaluate_path(d,ids,order,al)
        # automatic stop: last merge with logBF>0
        auto=[(K,s) for K,v,s in path if v>0]
        aK,aS=(auto[-1] if auto else (len(ids),M.evaluate(al,d.cluster_id.to_numpy())))
        best=max(path,key=lambda t:t[2]['ami'])
        rows.append(dict(pool=pool,alleles=man.loc[pool,'allele_count'],a0=a0,
            auto_k=aK,auto_ami=aS['ami'],auto_f1=aS['bcubed_f1_macro'],auto_pur=aS['adjusted_purity_macro'],
            best_k=best[0],best_ami=best[2]['ami'],best_f1=best[2]['bcubed_f1_macro']))
    print('.',end='',flush=True)
print(f'\n{time.time()-t0:.0f}s')
r=pd.DataFrame(rows); r.to_csv('merge_proto.csv',index=False)
r['band']=pd.cut(r.alleles,bins=[1,6,12,30],labels=['2-6','7-12','13-30'])
for a0 in (1.0,10.0):
    s=r[r.a0==a0]
    print(f'\n===== prior concentration a0={a0} =====')
    print('OVERALL'); print(s[['auto_k','auto_ami','auto_f1','auto_pur','best_k','best_ami','best_f1']].mean().round(4).to_string())
    print('BY BAND'); print(s.groupby('band',observed=True)[['auto_k','auto_ami','auto_f1','best_k','best_ami','best_f1']].mean().round(4).to_string())

"""Control experiment: does the PepCluster2 seed do any work?

Identical EM (same frame, same pseudocounts, same convergence rule) from four inits:
  A  merge  : PepCluster2 -> Bayesian merge (a0=3, t=0) -> EM        [our pipeline]
  B  rand@A : random init at the SAME k that A's merge produced, R restarts
  C  frag   : PepCluster2's ~175 fragments -> EM directly (no merge stage)
  D  rand@K : random init at the TRUE allele count, R restarts   [parallel to
              MixMHCp forced-k; uses the label, reference only]

Random init seeds each component from a small random subset of peptides, which is
the favourable/fair random baseline (uniform-random responsibilities would make all
components start at background and is a strawman).
"""
import sys,time
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse
sys.argv=['x']
exec(open('bhc.py').read().split("ap=argparse")[0])

A0_MERGE=3.0; T=0.0; A0_EM=3.0; RESTARTS=10; SEED_PEPS=20

def em_counts(X,C,w,alpha_pc,n_iter=60,tol=1e-6):
    N=len(X); K=len(C)
    oh=[sparse.csr_matrix((np.ones(N),(np.arange(N),X[:,j])),shape=(N,21)) for j in range(NP)]
    prev=-np.inf
    for _ in range(n_iter):
        th=C+alpha_pc; th=th/th.sum(-1,keepdims=True)
        lt=np.zeros((K,NP,21)); lt[:,:,:20]=np.log(th)
        ll=np.zeros((N,K))
        for j in range(NP): ll+=lt[:,j,:][:,X[:,j]].T
        ll+=np.log(np.maximum(w,1e-300))
        mx=ll.max(1,keepdims=True); e=np.exp(ll-mx); s=e.sum(1,keepdims=True)
        obj=float((np.log(s)+mx).sum()); R=e/s
        w=R.sum(0); w=np.maximum(w,1e-12); w/=w.sum()
        C=np.zeros((K,NP,20))
        for j in range(NP): C[:,j,:]=(R.T@oh[j])[:,:20]
        if abs(obj-prev)<tol*abs(obj): break
        prev=obj
    return R.argmax(1),obj

def init_labels(X,lab,K):
    return count_matrix(X,lab,K), np.maximum(np.bincount(lab,minlength=K),1e-9)/len(lab)

def init_random(X,K,rng):
    N=len(X); m=min(SEED_PEPS,max(2,N//(2*K)))
    C=np.zeros((K,NP,20))
    for k in range(K):
        idx=rng.choice(N,size=m,replace=False)
        C[k]=count_matrix(X[idx],np.zeros(m,dtype=np.int64),1)[0]
    return C, np.full(K,1.0/K)

rows=[]; t0=time.time()
for f in sorted(A.glob('*.tsv')):
    pool=f.stem; d=pd.read_csv(f,sep='\t')
    al=d.allele.to_numpy(); X=encode(d.peptide.values); N=len(X)
    ids,inv=np.unique(d.cluster_id.values,return_inverse=True); K0=len(ids)
    C0=count_matrix(X,inv,K0); sz=np.bincount(inv,minlength=K0)
    bg=C0.sum((0,1)); bg/=bg.sum(); pc=A0_EM*bg*20
    Ktrue=int(man.loc[pool,'alleles'].count(';')+1) if False else int(man.loc[pool,'allele_count'])
    rec=dict(pool=pool,alleles=Ktrue,peptides=N,K0=K0,
             base_ami=M.evaluate(al,inv)['ami'])
    # A: merge -> EM
    order=agglomerate(C0,sz,make_lml(A0_MERGE*bg*20+1e-9),'bf',1.0)
    sc=np.array([s for _,_,s in order]); pos=np.where(sc>T)[0]
    lab=labels_at(order,K0,int(pos[-1]+1) if len(pos) else 0)[inv]
    kA=int(lab.max())+1
    Ci,wi=init_labels(X,lab,kA); la,_=em_counts(X,Ci,wi,pc)
    a=M.evaluate(al,la); rec.update(A_k_seed=kA,A_k=a['clusters'],A_ami=a['ami'],A_f1=a['bcubed_f1_macro'])
    # C: fragments -> EM, no merge
    Ci,wi=init_labels(X,inv,K0); lc,_=em_counts(X,Ci,wi,pc)
    c=M.evaluate(al,lc); rec.update(C_k=c['clusters'],C_ami=c['ami'],C_f1=c['bcubed_f1_macro'])
    # B and D: random restarts
    for tag,K in (('B',kA),('D',Ktrue)):
        am=[];fm=[];ob=[]
        for r in range(RESTARTS):
            rng=np.random.default_rng(hash((pool,tag,r))%(2**31))
            Ci,wi=init_random(X,K,rng); lr,o=em_counts(X,Ci,wi,pc)
            e=M.evaluate(al,lr); am.append(e['ami']); fm.append(e['bcubed_f1_macro']); ob.append(o)
        am=np.array(am); fm=np.array(fm); ob=np.array(ob)
        rec.update({f'{tag}_k':K,f'{tag}_ami_mean':am.mean(),f'{tag}_ami_std':am.std(),
                    f'{tag}_ami_best':am.max(),f'{tag}_ami_worst':am.min(),
                    f'{tag}_f1_mean':fm.mean(),
                    f'{tag}_ami_maxlik':float(am[ob.argmax()])})
    rows.append(rec); print('.',end='',flush=True)
print(f'\n{time.time()-t0:.0f}s')
pd.DataFrame(rows).to_csv(SP/'controls.csv',index=False); print('wrote controls.csv',len(rows))

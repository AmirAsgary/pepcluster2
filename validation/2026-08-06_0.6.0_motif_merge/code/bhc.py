"""Merge stage v3: proper Bayesian Hierarchical Clustering (Heller & Ghahramani 2005)
with a Dirichlet-multinomial likelihood, followed by EM refinement.

For a candidate merge of subtrees T_i, T_j into T_k (n_k = n_i + n_j):
    d_k    = alpha*Gamma(n_k) + d_i*d_j          (leaves: d = alpha)
    pi_k   = alpha*Gamma(n_k) / d_k
    P(D_k|T_k) = pi_k*P(D_k|H1) + (1-pi_k)*P(D_i|T_i)*P(D_j|T_j)
    r_k    = pi_k*P(D_k|H1) / P(D_k|T_k)         -> merge while r_k > 0.5

The split hypothesis uses the SUBTREE marginal P(D_i|T_i), not the single-cluster
marginal. That is the term my first attempt omitted, and why it collapsed to k=1.
Everything is in logs.  Merge score = log-odds; r_k > 0.5 <=> score > 0.
"""
import sys, time, argparse
from pathlib import Path
import numpy as np, pandas as pd
from scipy.special import gammaln
from scipy import sparse
from scipy.special import logsumexp

B=Path('/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2/validation/2026-07-30_0.5.0-dev_search_redesign')
sys.path.insert(0,str(B/'code/mhc_bench')); import metrics as M
SP=Path('/tmp/claude-847728/-mnt-vast-standard-home-a-hajialiasgarynaj01-u14286/d416d4b8-f030-4400-8270-87c634df447a/scratchpad')
A=SP/'assign'
man=pd.read_csv(B/'runs/mhc_bench_sep_kmer_anchor/pool_manifest.csv').set_index('pool')
AA='ACDEFGHIKLMNPQRSTVWY'; IDX={a:i for i,a in enumerate(AA)}; NP=9; MISS=20

def encode(peps):
    """9 columns. L>=9: pos 1-4 -> cols 1-4, pos L-4..L -> cols 5-9.
       L==8:  pos 1-4 -> cols 1-4, GAP at col 5, pos 5-8 -> cols 6-9."""
    X=np.full((len(peps),NP),MISS,dtype=np.int64)
    for i,p in enumerate(peps):
        L=len(p)
        if   L>=NP: src=[0,1,2,3,L-5,L-4,L-3,L-2,L-1]; cols=[0,1,2,3,4,5,6,7,8]
        elif L==8:  src=[0,1,2,3,4,5,6,7];             cols=[0,1,2,3,5,6,7,8]
        else: continue
        for c,s in zip(cols,src): X[i,c]=IDX.get(p[s],MISS)
    return X

def count_matrix(X,lab,K):
    C=np.zeros((K,NP,20))
    for j in range(NP):
        ok=X[:,j]<20
        if ok.any(): np.add.at(C[:,j,:],(lab[ok],X[ok,j]),1.0)
    return C

def make_lml(alpha):
    A0=alpha.sum(); ga=gammaln(alpha).sum()
    def lml(x):
        n=x.sum(-1)
        return (gammaln(A0)-gammaln(A0+n)+gammaln(x+alpha).sum(-1)-ga).sum(-1)
    return lml

def agglomerate(C0,sz,lml,mode,crp_alpha):
    """mode 'bf': score = logBF.   mode 'bhc': Heller-Ghahramani log-odds."""
    K=len(C0); C=C0.copy(); n=sz.astype(float).copy()
    H1=lml(C)                       # log P(D|H1) for each current cluster
    la=np.log(crp_alpha)
    logd=np.full(K,la)              # leaves: d = alpha
    logPT=H1.copy()                 # leaves: P(D|T) = P(D|H1)
    def score(i,oth):
        h1=lml(C[i][None]+C[oth])
        if mode=='bf': return h1-H1[i]-H1[oth]
        return (la+gammaln(n[i]+n[oth])-logd[i]-logd[oth]) + h1 - logPT[i] - logPT[oth]
    S=np.full((K,K),-np.inf)
    for i in range(K):
        j=np.arange(i+1,K)
        if len(j): S[i,j]=score(i,j)
    alive=np.ones(K,bool); order=[]
    while alive.sum()>1:
        m=S.copy(); m[~alive,:]=-np.inf; m[:,~alive]=-np.inf
        i,j=np.unravel_index(np.argmax(m),m.shape); v=float(m[i,j])
        order.append((int(i),int(j),v))
        nk=n[i]+n[j]; Ck=C[i]+C[j]; h1k=lml(Ck)
        if mode=='bhc':
            lag=la+gammaln(nk); ldd=logd[i]+logd[j]
            ldk=np.logaddexp(lag,ldd)
            lpi=lag-ldk; l1mpi=ldd-ldk
            logPT[i]=np.logaddexp(lpi+h1k, l1mpi+logPT[i]+logPT[j]); logd[i]=ldk
        C[i]=Ck; n[i]=nk; H1[i]=h1k; alive[j]=False
        S[j,:]=-np.inf; S[:,j]=-np.inf
        oth=np.where(alive)[0]; oth=oth[oth!=i]
        if len(oth): S[np.minimum(i,oth),np.maximum(i,oth)]=score(i,oth)
    return order

def labels_at(order,K0,steps):
    par=np.arange(K0)
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for i,j,_ in order[:steps]: par[find(j)]=find(i)
    root=np.array([find(x) for x in range(K0)])
    _,comp=np.unique(root,return_inverse=True); return comp

def em(X,lab0,alpha_pc,n_iter=60,tol=1e-6):
    N=len(X); K=int(lab0.max())+1
    oh=[sparse.csr_matrix((np.ones(N),(np.arange(N),X[:,j])),shape=(N,21)) for j in range(NP)]
    C=count_matrix(X,lab0,K)
    w=np.bincount(lab0,minlength=K).astype(float); w=np.maximum(w,1e-9); w/=w.sum()
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
    return R.argmax(1)

ap=argparse.ArgumentParser()
ap.add_argument('--pools',type=int,default=0)
ap.add_argument('--a0',type=float,nargs='+',default=[1.0,10.0])
ap.add_argument('--alpha',type=float,nargs='+',default=[1e-2,1.0,1e2,1e4,1e6])
ap.add_argument('--out',default='bhc.csv')
args=ap.parse_args()

files=sorted(A.glob('*.tsv'))
if args.pools: files=files[:args.pools]
rows=[]; t0=time.time()
for f in files:
    pool=f.stem; d=pd.read_csv(f,sep='\t')
    al=d.allele.to_numpy(); X=encode(d.peptide.values)
    ids,inv=np.unique(d.cluster_id.values,return_inverse=True); K0=len(ids)
    C0=count_matrix(X,inv,K0); sz=np.bincount(inv,minlength=K0)
    bg=C0.sum((0,1)); bg=bg/bg.sum()
    base=M.evaluate(al,inv)
    # oracle ceiling for this pool
    maj=d.groupby('cluster_id').allele.agg(lambda s:s.value_counts().idxmax())
    orc=M.evaluate(al,d.cluster_id.map(maj).to_numpy())
    for a0 in args.a0:
        lml=make_lml(a0*bg*20+1e-9)
        cfgs=[('bf',1.0)]+[('bhc',a) for a in args.alpha]
        for mode,ca in cfgs:
            order=agglomerate(C0,sz,lml,mode,ca)
            sc=np.array([s for _,_,s in order])
            pos=np.where(sc>0)[0]; steps=int(pos[-1]+1) if len(pos) else 0
            comp=labels_at(order,K0,steps); lab=comp[inv]
            mg=M.evaluate(al,lab)
            rec=dict(pool=pool,alleles=man.loc[pool,'allele_count'],peptides=len(d),
                a0=a0,mode=mode,alpha=ca,
                base_ami=base['ami'],base_f1=base['bcubed_f1_macro'],base_k=base['clusters'],
                orc_ami=orc['ami'],orc_f1=orc['bcubed_f1_macro'],
                mg_k=mg['clusters'],mg_ami=mg['ami'],mg_f1=mg['bcubed_f1_macro'],
                mg_pur=mg['adjusted_purity_macro'],mg_rec=mg['bcubed_recall_macro'])
            if lab.max()>0:
                le=em(X,lab,a0*bg*20)
                e=M.evaluate(al,le)
                rec.update(em_k=e['clusters'],em_ami=e['ami'],em_f1=e['bcubed_f1_macro'],
                           em_pur=e['adjusted_purity_macro'],em_rec=e['bcubed_recall_macro'])
            else:
                rec.update(em_k=mg['clusters'],em_ami=mg['ami'],em_f1=mg['bcubed_f1_macro'],
                           em_pur=mg['adjusted_purity_macro'],em_rec=mg['bcubed_recall_macro'])
            rows.append(rec)
    print('.',end='',flush=True)
print(f'\n{time.time()-t0:.0f}s')
pd.DataFrame(rows).to_csv(SP/args.out,index=False); print('wrote',args.out,len(rows),'rows')

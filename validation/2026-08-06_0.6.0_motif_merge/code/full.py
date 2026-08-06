"""Full sweep. Two knobs, both of which are priors on granularity:
   a0 : Dirichlet prior concentration (per-position pseudocount total = a0*20)
   t  : merge threshold, stop when logBF <= t.  Equivalent to a partition prior
        P(partition) ~ exp(-t*k), i.e. a flat per-cluster penalty. Because t is
        constant across pairs it does NOT change the merge order, so one
        dendrogram per a0 serves the whole t sweep for free.
   BHC/CRP included at two alpha values to document the negative result.
"""
import sys,time
from pathlib import Path
import numpy as np, pandas as pd
sys.argv=['x']
exec(open('bhc.py').read().split("ap=argparse")[0])

A0S=[0.3,1.0,3.0]; TS=[0.0,20.0,60.0]; BHC_A=[1e2,1e6]
rows=[]; t0=time.time()
for f in sorted(A.glob('*.tsv')):
    pool=f.stem; d=pd.read_csv(f,sep='\t')
    al=d.allele.to_numpy(); X=encode(d.peptide.values)
    ids,inv=np.unique(d.cluster_id.values,return_inverse=True); K0=len(ids)
    C0=count_matrix(X,inv,K0); sz=np.bincount(inv,minlength=K0)
    bg=C0.sum((0,1)); bg/=bg.sum()
    base=M.evaluate(al,inv)
    maj=d.groupby('cluster_id').allele.agg(lambda s:s.value_counts().idxmax())
    orc=M.evaluate(al,d.cluster_id.map(maj).to_numpy())
    common=dict(pool=pool,alleles=man.loc[pool,'allele_count'],peptides=len(d),
        base_ami=base['ami'],base_f1=base['bcubed_f1_macro'],base_k=base['clusters'],
        orc_ami=orc['ami'],orc_f1=orc['bcubed_f1_macro'],orc_k=orc['clusters'])
    def emit(mode,a0,alpha,t,order,sc):
        pos=np.where(sc>t)[0]; steps=int(pos[-1]+1) if len(pos) else 0
        lab=labels_at(order,K0,steps)[inv]; mg=M.evaluate(al,lab)
        r=dict(common,mode=mode,a0=a0,alpha=alpha,t=t,
               mg_k=mg['clusters'],mg_ami=mg['ami'],mg_f1=mg['bcubed_f1_macro'],
               mg_pur=mg['adjusted_purity_macro'],mg_rec=mg['bcubed_recall_macro'])
        if lab.max()>0:
            e=M.evaluate(al,em(X,lab,a0*bg*20))
            r.update(em_k=e['clusters'],em_ami=e['ami'],em_f1=e['bcubed_f1_macro'],
                     em_pur=e['adjusted_purity_macro'],em_rec=e['bcubed_recall_macro'])
        else:
            r.update(em_k=mg['clusters'],em_ami=mg['ami'],em_f1=mg['bcubed_f1_macro'],
                     em_pur=mg['adjusted_purity_macro'],em_rec=mg['bcubed_recall_macro'])
        rows.append(r)
    for a0 in A0S:
        lml=make_lml(a0*bg*20+1e-9)
        order=agglomerate(C0,sz,lml,'bf',1.0); sc=np.array([s for _,_,s in order])
        for t in TS: emit('bf',a0,np.nan,t,order,sc)
        if a0==1.0:
            for ca in BHC_A:
                o2=agglomerate(C0,sz,lml,'bhc',ca); s2=np.array([s for _,_,s in o2])
                emit('bhc',a0,ca,0.0,o2,s2)
    print('.',end='',flush=True)
print(f'\n{time.time()-t0:.0f}s')
pd.DataFrame(rows).to_csv(SP/'full.csv',index=False); print('wrote full.csv',len(rows))

import json, subprocess, sys, tempfile, time, concurrent.futures
from pathlib import Path
import numpy as np, pandas as pd

B=Path('/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2/validation/2026-07-30_0.5.0-dev_search_redesign')
sys.path.insert(0,str(B/'code/mhc_bench')); import metrics as M
BIN=Path('/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2/target/release/pepcluster2')
POOLS=B/'runs/mhc_bench_sep_kmer_anchor/pools'
TMP=Path('/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pc2_bench_tmp'); TMP.mkdir(exist_ok=True)
OUT=Path('/tmp/claude-847728/-mnt-vast-standard-home-a-hajialiasgarynaj01-u14286/d416d4b8-f030-4400-8270-87c634df447a/scratchpad/assign'); OUT.mkdir(exist_ok=True)

man=pd.read_csv(B/'runs/mhc_bench_sep_kmer_anchor/pool_manifest.csv')
test=man[man.split=='test'].sort_values('peptides',ascending=False)

def run(pool):
    dst=OUT/f'{pool}.tsv'
    if dst.exists(): return pool,'cached'
    fasta=POOLS/f'{pool}.fasta'
    with tempfile.TemporaryDirectory(prefix='oracle_',dir=TMP) as tmp:
        o=Path(tmp)/'out'
        cmd=[str(BIN),'--input',str(fasta),'--output-dir',str(o),
             '--mode','separate_kmer_anchor','--kmer-similarity-threshold','0.25',
             '--anchor-combination-similarity-threshold','0.35',
             '--representative-order','coverage','--threads','2',
             '--candidate-buffer-mb','512','--compact-output','--tmp-dir',str(Path(tmp)/'tmp'),
             '--clustering-method','graph','--no-prefilter']
        r=subprocess.run(cmd,capture_output=True,text=True)
        if r.returncode: return pool,f'FAIL {(r.stderr or r.stdout)[-200:]}'
        a=pd.read_csv(o/'node_clusters.tsv',sep='\t')
        lab=pd.read_csv(POOLS/f'{pool}.labels.tsv',sep='\t')
        m=lab.merge(a[['sequence','cluster_id']],left_on='peptide',right_on='sequence',how='inner')
        if len(m)!=len(lab): return pool,f'FAIL merge {len(m)}/{len(lab)}'
        m[['peptide','allele','cluster_id']].to_csv(dst,sep='\t',index=False)
    return pool,'ok'

t0=time.time()
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
    for pool,st in ex.map(run,test['pool'].tolist()):
        if st!='ok' and st!='cached': print(pool,st,flush=True)
print(f'clustering done in {time.time()-t0:.0f}s',flush=True)

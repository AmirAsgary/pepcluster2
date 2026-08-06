# External tool benchmark

PepCluster2 against published motif-deconvolution tools, on the same pools and
the same metrics.

## Layout

```
tools/                  installed third-party tools
  env/                  conda env (perl + List::MoreUtils) MixMHCp needs
  MixMHCp/              installed and working
  gibbscluster/         not installed - see below
code/
  run_mixmhcp.py        run + score MixMHCp
  run_gibbscluster.py   run + score GibbsCluster (parser untested, see below)
  compare_tools.py      cross-tool tables, figure and REPORT.md
results/
  immuneapp/            one directory per dataset
    raw/                each tool's own per-pool output
    tables/             per_pool.csv, summary.csv
    REPORT.md
    tool_comparison.png
```

A second dataset goes in `results/<its-name>/` and needs no code change: the
runners take `--pools`, `--manifest` and `--out`, and `compare_tools.py` takes
`--results` and `--dataset`.

## Datasets

`immuneapp` - the pools this study already used, built from the ImmuneApp
supplementary tables: 193 pools (120 inner, 25 outer, 48 test), 972-24,888
peptides, 2-30 alleles, balanced to within 5%, no two alleles above 97% similar.
Unchanged from the tuning analysis, so the PepCluster2 numbers here are the same
ones its own report quotes.

## Reproducing

```bash
export PATH="$PWD/tools/env/bin:$PATH"
RUNS=../runs/mhc_bench_sep_kmer_anchor

python3 code/run_mixmhcp.py \
  --pools $RUNS/pools --manifest $RUNS/pool_manifest.csv \
  --binary tools/MixMHCp/MixMHCp \
  --out results/immuneapp/raw/mixmhcp.csv \
  --cores 48 --tmp-root /cbscratch/amirasgary2/pc2_bench_tmp

python3 code/compare_tools.py --runs ../runs \
  --results results/immuneapp --dataset immuneapp
```

## Settings

Each tool is run twice.

`default` is the tool as documented, letting it choose its own number of motifs
or groups. This is the fair headline number.

`forced k` fixes that number to the pool's true allele count. No user could do
this without already knowing the answer, so it is not a fair comparison — it is
reported to show how much of a tool's result comes from its model rather than
from its model selection.

## Installing MixMHCp

Already done. `git clone https://github.com/GfellerLab/MixMHCp`, set `lib_path`
in the `MixMHCp` script to the absolute `lib/` path, then
`g++ -O2 lib/MixMHCp.cc -o lib/MixMHCp.x`. It also needs perl with
`List::MoreUtils`, which is why `tools/env` exists.

MixMHCp exits non-zero when its optional R logo/length plots fail, long after
clustering has succeeded. The runner therefore ignores the exit code and gates on
the responsibility table being present and covering every input peptide.

## Installing GibbsCluster — needs you

GibbsCluster is not freely downloadable. It sits behind a DTU academic licence
form that asks for name, email, affiliation and position, and requires accepting
the licence. That is a legal agreement tied to a person, so it has to be you:

1. https://services.healthtech.dtu.dk/services/GibbsCluster-2.0/ → "Downloads"
2. Choose version 2.0f, platform Linux, fill the form and accept the licence
3. Unpack the tarball into `tools/gibbscluster/`
4. Edit `gibbscluster` and set `setenv GIBBS` to that directory, per its readme
5. Confirm the parser before a full sweep:

```bash
python3 code/run_gibbscluster.py --selftest \
  --pools $RUNS/pools --manifest $RUNS/pool_manifest.csv \
  --binary tools/gibbscluster/gibbscluster \
  --out /dev/null --tmp-root /cbscratch/amirasgary2/pc2_bench_tmp
```

The self-test runs the smallest pool and prints the parsed cluster count and
scores. `run_gibbscluster.py` locates the peptide and cluster columns by header
name rather than position, so a column-order change fails loudly instead of
being mis-parsed — but the parser has never been executed against the real tool,
so do not skip the self-test.

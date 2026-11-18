# Source data

`pep_hla_dataset.parquet` is not tracked. It is a mono-allelic MS
immunopeptidomics table supplied for this study: 412,644 peptides over 91 HLA
alleles and 6 studies, one row per peptide-per-cell-line observation.

Columns used by `code/build_pools.py`:

| column | use |
|---|---|
| `peptide` | the sequence; isoleucine is already collapsed to leucine upstream |
| `modifiedSequence` | checked only, not used: stripping its bracketed mass deltas returns `peptide` exactly for all 412,644 rows, so it cannot recover I/L |
| `cell_context` | `<study>-<instrument>-<allele>`; the allele is the trailing token and becomes the label |

Place the file here and run `python3 code/build_pools.py` to regenerate the
160 pools. The build is deterministic: seeds are `zlib.crc32` of each pool
name, recorded in `pool_manifest.csv`.

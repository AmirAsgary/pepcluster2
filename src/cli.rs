use crate::scoring::{ScoringMode, SimdMode};
use std::env;
use std::path::PathBuf;

pub const VERSION: &str = "0.4.3";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PrefilterChoice {
    Auto,
    Force,
    Disable,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ClusteringMethod {
    Graph,
    Greedy,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GreedySelection {
    KmerDegree,
    LazyExact,
}

impl GreedySelection {
    pub fn name(self) -> &'static str {
        match self {
            Self::KmerDegree => "kmer-degree",
            Self::LazyExact => "lazy-exact",
        }
    }
}

impl ClusteringMethod {
    pub fn name(self) -> &'static str {
        match self {
            Self::Graph => "graph",
            Self::Greedy => "greedy",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Config {
    pub input: PathBuf,
    pub output_dir: PathBuf,
    pub tmp_dir: PathBuf,
    pub mode: ScoringMode,
    pub clustering_method: ClusteringMethod,
    pub greedy_selection: GreedySelection,
    pub threshold: f64,
    pub threshold_explicit: bool,
    pub alignment_threshold: f64,
    pub anchor_threshold: f64,
    pub alignment_threshold_overridden: bool,
    pub anchor_threshold_overridden: bool,
    pub gap_open: f64,
    pub gap_extension: f64,
    pub terminal_gap_open: f64,
    pub terminal_gap_extension: f64,
    pub minimum_terminal_match_length: usize,
    pub kmer_seed_threshold: f64,
    pub kmer_table: PathBuf,
    pub threads: usize,
    pub iteration_cap: Option<usize>,
    pub merge_cap: Option<usize>,
    pub prefilter: PrefilterChoice,
    pub full_sensitive_after_prefilter: bool,
    pub min_cluster_size: usize,
    pub merge: bool,
    pub strict: bool,
    pub index_only: bool,
    pub keep_tmp: bool,
    pub write_edges: bool,
    pub candidate_buffer_mb: usize,
    pub max_memory_gb: f64,
    pub simd_mode: SimdMode,
    pub compact_output: bool,
    pub write_cluster_fastas: bool,
    pub write_scored_pairs: bool,
}

fn help() -> &'static str {
    "pepcluster2 0.4.3
Experimental shift-aware clustering for MHC-I peptides.

USAGE:
  pepcluster2 -i INPUT.fasta -o OUTPUT_DIR [OPTIONS]

REQUIRED:
  -i, --input PATH              Input FASTA (plain text)
  -o, --output-dir PATH         Output directory

SCORING:
      --mode NAME               combined_kmer_anchor|combined_full_anchor|separate_aln_anchor
                                [default: separate_aln_anchor]
  -t, --threshold FLOAT         Combined-mode threshold; when explicitly supplied in
                                separate mode, sets both component thresholds
                                [default combined threshold: 0.60]
      --alignment-similarity-threshold FLOAT
                                Override the alignment threshold in separate mode
                                [default: 0.50]
      --anchor-combination-similarity-threshold FLOAT
                                Override the anchor threshold in separate mode
                                [default: 0.60]
      --kmer-seed-threshold F   Similar 2-mer retrieval threshold [default: 0.50]

ALIGNMENT:
      --gap-open FLOAT          Internal affine gap-open penalty [default: -4]
      --gap-extension FLOAT     Internal gap-extension penalty [default: -1]
      --terminal-overhang-gap-open FLOAT
                                Terminal-overhang gap-open penalty [default: -2]
      --terminal-overhang-gap-extension FLOAT
                                Terminal-overhang extension penalty [default: -1]
      --minimum-terminal-match-length INT
                                Required matched columns at each terminus [default: 2]

CLUSTERING:
      --clustering-method NAME  graph|greedy [default: graph]
      --greedy-selection NAME   kmer-degree|lazy-exact [default: kmer-degree]
      --iteration-cap INT       Maximum iterative passes [default: unset]
      --merge-cap INT           Early merge-rejection sample size [default: all]
      --min-cluster-size INT    Minimum size for per-cluster FASTA [default: 2]
      --no-merge                Disable strict representative-covering merge

GRAPH PREFILTER:
      --force-prefilter         Always run graph high-confidence prefilter
      --no-prefilter            Disable graph prefilter
      --full-sensitive-after-prefilter
                                Score all sensitive graph candidates after prefilter

PERFORMANCE AND STORAGE:
  -p, --threads INT             Worker threads; 0 uses all available CPUs [default: 0]
      --simd auto|on|off        Reserved scoring selector [default: auto]
      --tmp-dir PATH            Temporary directory [default: OUTPUT_DIR/tmp]
      --kmer-table PATH         Cached binary 2-mer table [default: TMP_DIR/kmer2_similarity_q.bin]
      --candidate-buffer-mb INT Candidate spill buffer [default: 256]
      --max-memory-gb FLOAT     Maximum estimated graph memory [default: 8]
      --keep-tmp                Keep temporary graph files
      --write-edges             Write graph edges.tsv
      --compact-output          Write compact node assignments
      --write-cluster-fastas    Write per-cluster FASTA files
      --write-scored-pairs      Validation only: write unique exactly-scored pairs

VALIDATION:
      --index-only              Count k-mer candidates without exact scoring
      --strict                  Fail instead of skipping invalid records
  -h, --help                    Print help
  -V, --version                 Print version
"
}

fn next_value(args: &[String], i: &mut usize, flag: &str) -> Result<String, String> {
    *i += 1;
    args.get(*i)
        .cloned()
        .ok_or_else(|| format!("missing value after {flag}"))
}

fn parse_number<T: std::str::FromStr>(value: String, flag: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("invalid value for {flag}: {value}"))
}

pub fn parse() -> Result<Option<Config>, String> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() {
        println!("{}", help());
        return Ok(None);
    }
    let mut input = None;
    let mut output_dir = None;
    let mut tmp_dir = None;
    let mut mode = ScoringMode::SeparateAlnAnchor;
    let mut clustering_method = ClusteringMethod::Graph;
    let mut greedy_selection = GreedySelection::KmerDegree;
    let mut threshold = 0.60;
    let mut threshold_explicit = false;
    let mut alignment_threshold = None;
    let mut anchor_threshold = None;
    let mut gap_open: f64 = -4.0;
    let mut gap_extension: f64 = -1.0;
    let mut terminal_gap_open: f64 = -2.0;
    let mut terminal_gap_extension: f64 = -1.0;
    let mut minimum_terminal_match_length = 2usize;
    let mut kmer_seed_threshold = 0.50;
    let mut kmer_table = None;
    let mut threads = 0usize;
    let mut iteration_cap = None;
    let mut merge_cap = None;
    let mut prefilter = PrefilterChoice::Auto;
    let mut full_sensitive_after_prefilter = false;
    let mut min_cluster_size = 2usize;
    let mut merge = true;
    let mut strict = false;
    let mut index_only = false;
    let mut keep_tmp = false;
    let mut write_edges = false;
    let mut candidate_buffer_mb = 256usize;
    let mut max_memory_gb = 8.0;
    let mut simd_mode = SimdMode::Auto;
    let mut compact_output = false;
    let mut write_cluster_fastas = false;
    let mut write_scored_pairs = false;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                println!("{}", help());
                return Ok(None);
            }
            "-V" | "--version" => {
                println!("pepcluster2 {VERSION}");
                return Ok(None);
            }
            "-i" | "--input" => input = Some(PathBuf::from(next_value(&args, &mut i, "--input")?)),
            "-o" | "--output-dir" => {
                output_dir = Some(PathBuf::from(next_value(&args, &mut i, "--output-dir")?))
            }
            "--tmp-dir" => tmp_dir = Some(PathBuf::from(next_value(&args, &mut i, "--tmp-dir")?)),
            "--mode" => {
                mode = ScoringMode::parse(&next_value(&args, &mut i, "--mode")?)?;
            }
            "--clustering-method" => {
                clustering_method = match next_value(&args, &mut i, "--clustering-method")?.as_str()
                {
                    "graph" => ClusteringMethod::Graph,
                    "greedy" => ClusteringMethod::Greedy,
                    value => return Err(format!("invalid --clustering-method: {value}")),
                };
            }
            "--greedy-selection" => {
                greedy_selection = match next_value(&args, &mut i, "--greedy-selection")?.as_str() {
                    "kmer-degree" => GreedySelection::KmerDegree,
                    "lazy-exact" => GreedySelection::LazyExact,
                    value => return Err(format!("invalid --greedy-selection: {value}")),
                };
            }
            "-t" | "--threshold" => {
                threshold = parse_number(next_value(&args, &mut i, "--threshold")?, "--threshold")?;
                threshold_explicit = true;
            }
            "--alignment-similarity-threshold" => {
                alignment_threshold = Some(parse_number(
                    next_value(&args, &mut i, "--alignment-similarity-threshold")?,
                    "--alignment-similarity-threshold",
                )?)
            }
            "--anchor-combination-similarity-threshold" => {
                anchor_threshold = Some(parse_number(
                    next_value(&args, &mut i, "--anchor-combination-similarity-threshold")?,
                    "--anchor-combination-similarity-threshold",
                )?)
            }
            "--gap-open" => {
                gap_open = parse_number(next_value(&args, &mut i, "--gap-open")?, "--gap-open")?
            }
            "--gap-extension" => {
                gap_extension = parse_number(
                    next_value(&args, &mut i, "--gap-extension")?,
                    "--gap-extension",
                )?
            }
            "--terminal-overhang-gap-open" => {
                terminal_gap_open = parse_number(
                    next_value(&args, &mut i, "--terminal-overhang-gap-open")?,
                    "--terminal-overhang-gap-open",
                )?
            }
            "--terminal-overhang-gap-extension" => {
                terminal_gap_extension = parse_number(
                    next_value(&args, &mut i, "--terminal-overhang-gap-extension")?,
                    "--terminal-overhang-gap-extension",
                )?
            }
            "--minimum-terminal-match-length" => {
                minimum_terminal_match_length = parse_number(
                    next_value(&args, &mut i, "--minimum-terminal-match-length")?,
                    "--minimum-terminal-match-length",
                )?
            }
            "--kmer-seed-threshold" => {
                kmer_seed_threshold = parse_number(
                    next_value(&args, &mut i, "--kmer-seed-threshold")?,
                    "--kmer-seed-threshold",
                )?
            }
            "--kmer-table" => {
                kmer_table = Some(PathBuf::from(next_value(&args, &mut i, "--kmer-table")?))
            }
            "-p" | "--threads" => {
                threads = parse_number(next_value(&args, &mut i, "--threads")?, "--threads")?
            }
            "--iteration-cap" | "--iterations" => {
                iteration_cap = Some(parse_number(
                    next_value(&args, &mut i, "--iteration-cap")?,
                    "--iteration-cap",
                )?)
            }
            "--merge-cap" => {
                merge_cap = Some(parse_number(
                    next_value(&args, &mut i, "--merge-cap")?,
                    "--merge-cap",
                )?)
            }
            "--min-cluster-size" => {
                min_cluster_size = parse_number(
                    next_value(&args, &mut i, "--min-cluster-size")?,
                    "--min-cluster-size",
                )?
            }
            "--candidate-buffer-mb" | "--edge-buffer-mb" => {
                candidate_buffer_mb = parse_number(
                    next_value(&args, &mut i, "--candidate-buffer-mb")?,
                    "--candidate-buffer-mb",
                )?
            }
            "--max-memory-gb" => {
                max_memory_gb = parse_number(
                    next_value(&args, &mut i, "--max-memory-gb")?,
                    "--max-memory-gb",
                )?
            }
            "--simd" => {
                simd_mode = match next_value(&args, &mut i, "--simd")?.as_str() {
                    "auto" => SimdMode::Auto,
                    "on" => SimdMode::On,
                    "off" => SimdMode::Off,
                    value => return Err(format!("invalid --simd mode: {value}")),
                }
            }
            "--no-merge" => merge = false,
            "--force-prefilter" => {
                if prefilter == PrefilterChoice::Disable {
                    return Err("--force-prefilter conflicts with --no-prefilter".into());
                }
                prefilter = PrefilterChoice::Force;
            }
            "--no-prefilter" => {
                if prefilter == PrefilterChoice::Force {
                    return Err("--no-prefilter conflicts with --force-prefilter".into());
                }
                prefilter = PrefilterChoice::Disable;
            }
            "--full-sensitive-after-prefilter" => full_sensitive_after_prefilter = true,
            "--strict" => strict = true,
            "--index-only" => index_only = true,
            "--keep-tmp" => keep_tmp = true,
            "--write-edges" => write_edges = true,
            "--compact-output" => compact_output = true,
            "--write-cluster-fastas" => write_cluster_fastas = true,
            "--write-scored-pairs" => write_scored_pairs = true,
            value => return Err(format!("unknown argument: {value}\n\n{}", help())),
        }
        i += 1;
    }

    let input = input.ok_or("--input is required")?;
    let output_dir = output_dir.ok_or("--output-dir is required")?;
    let tmp_dir = tmp_dir.unwrap_or_else(|| output_dir.join("tmp"));
    let kmer_table = kmer_table.unwrap_or_else(|| tmp_dir.join("kmer2_similarity_q.bin"));
    let alignment_threshold_overridden = alignment_threshold.is_some();
    let anchor_threshold_overridden = anchor_threshold.is_some();
    let alignment_threshold = alignment_threshold.unwrap_or_else(|| {
        if mode == ScoringMode::SeparateAlnAnchor && !threshold_explicit {
            0.50
        } else {
            threshold
        }
    });
    let anchor_threshold = anchor_threshold.unwrap_or(threshold);

    for (flag, value) in [
        ("--threshold", threshold),
        ("--alignment-similarity-threshold", alignment_threshold),
        (
            "--anchor-combination-similarity-threshold",
            anchor_threshold,
        ),
        ("--kmer-seed-threshold", kmer_seed_threshold),
    ] {
        if !(0.0..=1.0).contains(&value) {
            return Err(format!("{flag} must be between 0 and 1"));
        }
    }
    for (flag, value) in [
        ("--gap-open", gap_open),
        ("--gap-extension", gap_extension),
        ("--terminal-overhang-gap-open", terminal_gap_open),
        ("--terminal-overhang-gap-extension", terminal_gap_extension),
    ] {
        if !value.is_finite() || value > 0.0 {
            return Err(format!("{flag} must be a finite non-positive number"));
        }
    }
    if minimum_terminal_match_length == 0 || minimum_terminal_match_length > 3 {
        return Err("--minimum-terminal-match-length must be between 1 and 3".into());
    }
    if iteration_cap == Some(0) {
        return Err("--iteration-cap must be at least 1".into());
    }
    if merge_cap == Some(0) {
        return Err("--merge-cap must be at least 1".into());
    }
    if candidate_buffer_mb == 0 {
        return Err("--candidate-buffer-mb must be at least 1".into());
    }
    if min_cluster_size == 0 {
        return Err("--min-cluster-size must be at least 1".into());
    }
    if max_memory_gb <= 0.0 {
        return Err("--max-memory-gb must be positive".into());
    }
    if clustering_method == ClusteringMethod::Greedy {
        if prefilter == PrefilterChoice::Force {
            return Err(
                "--force-prefilter is available only with --clustering-method graph".into(),
            );
        }
        if full_sensitive_after_prefilter {
            return Err(
                "--full-sensitive-after-prefilter is available only with --clustering-method graph"
                    .into(),
            );
        }
        if write_edges {
            return Err("--write-edges requires --clustering-method graph".into());
        }
    } else if greedy_selection != GreedySelection::KmerDegree {
        return Err("--greedy-selection is available only with --clustering-method greedy".into());
    }

    Ok(Some(Config {
        input,
        output_dir,
        tmp_dir,
        mode,
        clustering_method,
        greedy_selection,
        threshold,
        threshold_explicit,
        alignment_threshold,
        anchor_threshold,
        alignment_threshold_overridden,
        anchor_threshold_overridden,
        gap_open,
        gap_extension,
        terminal_gap_open,
        terminal_gap_extension,
        minimum_terminal_match_length,
        kmer_seed_threshold,
        kmer_table,
        threads,
        iteration_cap,
        merge_cap,
        prefilter,
        full_sensitive_after_prefilter,
        min_cluster_size,
        merge,
        strict,
        index_only,
        keep_tmp,
        write_edges,
        candidate_buffer_mb,
        max_memory_gb,
        simd_mode,
        compact_output,
        write_cluster_fastas,
        write_scored_pairs,
    }))
}

use crate::graph::RepresentativeOrder;
use crate::index::TerminalSeed;
use crate::motif::MotifParams;
use crate::scoring::{ScoringMode, SimdMode};
use std::env;
use std::path::PathBuf;

pub const VERSION: &str = "0.7.0";

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
    /// Terminal k-mer similarity threshold, used by `separate_kmer_anchor`.
    pub kmer_threshold: f64,
    pub alignment_threshold_overridden: bool,
    pub anchor_threshold_overridden: bool,
    pub kmer_threshold_overridden: bool,
    pub gap_open: f64,
    pub gap_extension: f64,
    pub terminal_gap_open: f64,
    pub terminal_gap_extension: f64,
    pub minimum_terminal_match_length: usize,
    pub kmer_seed_threshold: f64,
    pub terminal_seed: TerminalSeed,
    pub representative_order: RepresentativeOrder,
    /// Hysteresis on synchronous reassignment, in similarity units.
    pub reassignment_margin: f64,
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
    /// Build the optional motif layer above the similarity clusters.
    pub merge_motifs: bool,
    pub motif: MotifParams,
}

fn help() -> &'static str {
    "pepcluster2 0.7.0
Experimental shift-aware clustering for MHC-I peptides.

USAGE:
  pepcluster2 -i INPUT.fasta -o OUTPUT_DIR [OPTIONS]

REQUIRED:
  -i, --input PATH              Input FASTA (plain text)
  -o, --output-dir PATH         Output directory

SCORING:
      --mode NAME               combined_kmer_anchor|combined_full_anchor|
                                separate_aln_anchor|separate_kmer_anchor
                                [default: separate_aln_anchor]
  -t, --threshold FLOAT         Combined-mode threshold; when explicitly supplied in
                                separate mode, sets both component thresholds
                                [default combined threshold: 0.60]
      --alignment-similarity-threshold FLOAT
                                Override the alignment threshold in separate mode
                                [default: 0.50]
      --anchor-combination-similarity-threshold FLOAT
                                Override the anchor threshold in either separate mode
                                [default: 0.60]
      --kmer-similarity-threshold FLOAT
                                Terminal k-mer similarity threshold, used by
                                separate_kmer_anchor. This is the similarity of the
                                first three and last three residues, compared
                                position by position with no alignment; the core
                                contributes nothing. Distinct from
                                --kmer-seed-threshold, which only retrieves
                                candidates and never accepts a relationship
                                [default: 0.60]
      --kmer-seed-threshold F   Similar 2-mer retrieval threshold [default: 0.40]
      --terminal-seed NAME      Indexed terminal column pairs:
                                all-column-pairs indexes (1,2), (1,3) and (2,3);
                                contiguous reproduces 0.4.3 and is insensitive to
                                shifted terminal columns
                                [default: all-column-pairs]

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
      --representative-order NAME
                                coverage selects representatives by dynamic set
                                cover; intrinsic visits peptides in an order that
                                depends only on the peptide itself, which makes
                                subset partitions nested [default: coverage]
      --reassignment-margin F   A peptide leaves its representative only when
                                another beats it by more than this margin. Zero
                                reproduces 0.4.3, where exact ties moved and made
                                clusters sensitive to dataset composition
                                [default: 0.01]
      --iteration-cap INT       Maximum iterative passes [default: unset]
      --merge-cap INT           Early merge-rejection sample size [default: all]
      --min-cluster-size INT    Minimum size for per-cluster FASTA [default: 2]
      --no-merge                Disable strict representative-covering merge

MOTIF LAYER (optional, off by default):
      --merge-motifs            Merge similarity clusters into motif-level groups
                                and write them as a separate output layer. A
                                similarity cluster is a ball around a
                                representative; a binding motif is a product of
                                per-position preferences, which one ball cannot
                                cover, so one motif fragments into many clusters.
                                This stage compares clusters as amino-acid
                                profiles and merges those a Dirichlet-multinomial
                                model says came from one profile. The motif
                                partition does NOT satisfy the
                                representative-to-member invariant and never
                                replaces the similarity clusters
      --motif-prior-concentration FLOAT
                                Dirichlet pseudocounts per motif column, divided
                                over the background residue frequencies. Larger
                                values smooth harder and merge more readily
                                [default: 10]
      --no-motif-merge          Skip the agglomerative merge and give EM one
                                component per similarity cluster, letting it find
                                its own count. Agreement is statistically
                                indistinguishable from merging first, but the
                                partition is finer: higher precision, lower
                                recall. A supported variant, not a diagnostic
      --motif-count INT         Seed EM with this many components and skip the
                                merge. Seeds are the similarity clusters that are
                                FARTHEST APART in profile space, not the largest:
                                the largest are drawn from far fewer distinct
                                motifs than their number suggests, so seeding on
                                size hands EM duplicates that it then merges.
                                Exactly this many motifs are returned. EM on its
                                own merges components the data does not separate,
                                so any it empties reclaim the peptide that fits
                                them best; motifs recovered that way can resemble
                                each other, which is the price of a strict count.
                                The alleles of a sample are usually known from
                                typing, so supplying this is ordinary use rather
                                than an oracle
      --motif-merge-threshold FLOAT
                                Merge while the best log Bayes factor exceeds
                                this. Equivalent to a prior over partitions
                                proportional to exp(-t * clusters), so larger
                                values keep more motifs [default: 25]
      --no-motif-em             DIAGNOSTIC ONLY. Skip EM and keep the merged
                                partition. Not a supported variant: merging alone
                                reaches AMI 0.43 against 0.60 with EM. Retained so
                                the published hyperparameter grid stays
                                reproducible
      --motif-em-prior-concentration FLOAT
                                Dirichlet pseudocounts smoothing the EM profiles.
                                This is the parameter that matters most: over the
                                swept range it moves AMI by 0.34, far more than
                                any merge setting [default: 3]
      --motif-em-max-iterations INT
                                EM iteration cap [default: 200]
      --motif-em-tolerance FLOAT
                                Relative log-likelihood change at which EM stops
                                [default: 1e-6]

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
      --max-memory-gb FLOAT     Memory ceiling. Caps the estimated graph memory,
                                and with --greedy-selection lazy-exact also caps
                                the cache of scored candidate lists. Past the cap
                                a list is not retained and the next pop recomputes
                                it, so a small value costs time, never accuracy.
                                Peak cache use is reported in run_stats.json
                                [default: 8]
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

impl Config {
    /// Reassignment margin quantized to the same thousandths as the ranking
    /// weight the refinement compares.
    pub fn reassignment_margin_q(&self) -> u16 {
        (self.reassignment_margin * 1000.0).round().clamp(0.0, 1000.0) as u16
    }
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
    let mut kmer_threshold = None;
    let mut gap_open: f64 = -4.0;
    let mut gap_extension: f64 = -1.0;
    let mut terminal_gap_open: f64 = -2.0;
    let mut terminal_gap_extension: f64 = -1.0;
    let mut minimum_terminal_match_length = 2usize;
    let mut kmer_seed_threshold = 0.40;
    let mut terminal_seed = TerminalSeed::AllColumnPairs;
    let mut representative_order = RepresentativeOrder::Coverage;
    let mut reassignment_margin = 0.01;
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
    let mut merge_motifs = false;
    // Selected by nested cross-validation on the inner folds; see
    // validation/2026-08-06_0.6.0_motif_merge/REPORT.md.
    let mut motif_prior_concentration = 10.0f64;
    let mut motif_merge_threshold = 25.0f64;
    let mut motif_em = true;
    let mut motif_merge = true;
    let mut motif_em_prior_concentration = 3.0f64;
    let mut motif_em_max_iterations = 200usize;
    let mut motif_em_tolerance = 1e-6f64;
    let mut motif_count: Option<usize> = None;

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
            "--kmer-similarity-threshold" | "--kmer_similarity_threshold" => {
                kmer_threshold = Some(parse_number(
                    next_value(&args, &mut i, "--kmer-similarity-threshold")?,
                    "--kmer-similarity-threshold",
                )?)
            }
            "--kmer-seed-threshold" => {
                kmer_seed_threshold = parse_number(
                    next_value(&args, &mut i, "--kmer-seed-threshold")?,
                    "--kmer-seed-threshold",
                )?
            }
            "--terminal-seed" => {
                terminal_seed = TerminalSeed::parse(&next_value(&args, &mut i, "--terminal-seed")?)?
            }
            "--representative-order" => {
                representative_order = RepresentativeOrder::parse(&next_value(
                    &args,
                    &mut i,
                    "--representative-order",
                )?)?
            }
            "--reassignment-margin" => {
                reassignment_margin = parse_number(
                    next_value(&args, &mut i, "--reassignment-margin")?,
                    "--reassignment-margin",
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
            "--merge-motifs" => merge_motifs = true,
            "--no-motif-em" => motif_em = false,
            "--no-motif-merge" => motif_merge = false,
            "--motif-prior-concentration" => {
                motif_prior_concentration = parse_number(
                    next_value(&args, &mut i, "--motif-prior-concentration")?,
                    "--motif-prior-concentration",
                )?
            }
            "--motif-count" => {
                motif_count = Some(parse_number(
                    next_value(&args, &mut i, "--motif-count")?,
                    "--motif-count",
                )?)
            }
            "--motif-merge-threshold" => {
                motif_merge_threshold = parse_number(
                    next_value(&args, &mut i, "--motif-merge-threshold")?,
                    "--motif-merge-threshold",
                )?
            }
            "--motif-em-prior-concentration" => {
                motif_em_prior_concentration = parse_number(
                    next_value(&args, &mut i, "--motif-em-prior-concentration")?,
                    "--motif-em-prior-concentration",
                )?
            }
            "--motif-em-max-iterations" => {
                motif_em_max_iterations = parse_number(
                    next_value(&args, &mut i, "--motif-em-max-iterations")?,
                    "--motif-em-max-iterations",
                )?
            }
            "--motif-em-tolerance" => {
                motif_em_tolerance = parse_number(
                    next_value(&args, &mut i, "--motif-em-tolerance")?,
                    "--motif-em-tolerance",
                )?
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
    let kmer_threshold_overridden = kmer_threshold.is_some();
    let alignment_threshold = alignment_threshold.unwrap_or_else(|| {
        if mode == ScoringMode::SeparateAlnAnchor && !threshold_explicit {
            0.50
        } else {
            threshold
        }
    });
    let anchor_threshold = anchor_threshold.unwrap_or(threshold);
    // With neither component threshold supplied, --threshold governs both.
    let kmer_threshold = kmer_threshold.unwrap_or(threshold);

    for (flag, value) in [
        ("--threshold", threshold),
        ("--alignment-similarity-threshold", alignment_threshold),
        (
            "--anchor-combination-similarity-threshold",
            anchor_threshold,
        ),
        ("--kmer-seed-threshold", kmer_seed_threshold),
        ("--kmer-similarity-threshold", kmer_threshold),
        ("--reassignment-margin", reassignment_margin),
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
    for (flag, value) in [
        ("--motif-prior-concentration", motif_prior_concentration),
        (
            "--motif-em-prior-concentration",
            motif_em_prior_concentration,
        ),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(format!("{flag} must be a finite positive number"));
        }
    }
    if !motif_merge_threshold.is_finite() {
        return Err("--motif-merge-threshold must be finite".into());
    }
    if !motif_em_tolerance.is_finite() || motif_em_tolerance < 0.0 {
        return Err("--motif-em-tolerance must be a finite non-negative number".into());
    }
    if motif_em_max_iterations == 0 {
        return Err("--motif-em-max-iterations must be at least 1".into());
    }
    if motif_count == Some(0) {
        return Err("--motif-count must be at least 1".into());
    }
    if motif_count.is_some() && !motif_merge {
        return Err(
            "--motif-count already bypasses the merge; --no-motif-merge is redundant"
                .into(),
        );
    }
    if !motif_merge && !motif_em {
        return Err("--no-motif-merge with --no-motif-em would do nothing".into());
    }
    if motif_count.is_some() && !motif_em {
        // The requested count seeds components from the largest clusters only;
        // without EM nothing would place the remaining peptides.
        return Err("--motif-count requires the EM stage; remove --no-motif-em".into());
    }
    if !merge_motifs {
        // Fail loudly rather than silently ignoring motif settings: a sweep that
        // forgets --merge-motifs would otherwise report identical results for
        // every configuration and look like the parameters do nothing.
        for flag in [
            "--motif-prior-concentration",
            "--motif-merge-threshold",
            "--motif-em-prior-concentration",
            "--motif-em-max-iterations",
            "--motif-em-tolerance",
            "--motif-count",
            "--no-motif-em",
            "--no-motif-merge",
        ] {
            if args.iter().any(|a| a == flag) {
                return Err(format!("{flag} requires --merge-motifs"));
            }
        }
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
    if representative_order == RepresentativeOrder::Intrinsic
        && clustering_method == ClusteringMethod::Greedy
        && greedy_selection == GreedySelection::LazyExact
    {
        return Err(
            "--greedy-selection lazy-exact is a dynamic set-cover rule and cannot be combined with --representative-order intrinsic; use --greedy-selection kmer-degree"
                .into(),
        );
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
        kmer_threshold,
        alignment_threshold_overridden,
        anchor_threshold_overridden,
        kmer_threshold_overridden,
        gap_open,
        gap_extension,
        terminal_gap_open,
        terminal_gap_extension,
        minimum_terminal_match_length,
        kmer_seed_threshold,
        terminal_seed,
        representative_order,
        reassignment_margin,
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
        merge_motifs,
        motif: MotifParams {
            merge: motif_merge,
            target_count: motif_count,
            prior_concentration: motif_prior_concentration,
            merge_threshold: motif_merge_threshold,
            em: motif_em,
            em_prior_concentration: motif_em_prior_concentration,
            em_max_iterations: motif_em_max_iterations,
            em_tolerance: motif_em_tolerance,
        },
    }))
}

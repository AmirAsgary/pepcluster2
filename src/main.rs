mod cli;
mod edge_store;
mod fasta;
mod graph;
mod greedy;
mod index;
mod kmer;
mod model;
mod motif;
mod output;
mod pair_trace;
mod scoring;

use cli::{ClusteringMethod, Config, GreedySelection, PrefilterChoice, VERSION};
use edge_store::{
    generate_edges, generate_prefilter_edges, merge_edge_files, CandidateScope, EdgeBuildStats,
    EdgeMode,
};
use fasta::{load_nodes, DynError, FastaStats};
use graph::{
    canonicalize, greedy_set_cover, iterate_to_convergence, load_graph, prefilter_masks,
    representative_coverage, Clustering, IterationStats,
};
use greedy::GreedyRunStats;
use index::{
    build_exact_index, build_similar_key_relations, retrieve_candidates, IndexStats, TerminalSeed,
};
use kmer::KmerSimilarityTable;
use model::Node;
use motif::MotifResult;
use output::{
    write_cluster_outputs, write_compact_outputs, write_edges, write_motif_outputs,
    write_provisional_clusters,
};
use rayon::prelude::*;
use rayon::ThreadPool;
use scoring::Scorer;
use std::env;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

fn json_escape(value: &str) -> String {
    format!(
        "\"{}\"",
        value
            .replace('\\', "\\\\")
            .replace('"', "\\\"")
            .replace('\n', "\\n")
    )
}

fn build_pool(requested: usize) -> Result<ThreadPool, DynError> {
    let builder = rayon::ThreadPoolBuilder::new();
    Ok(if requested == 0 {
        builder
    } else {
        builder.num_threads(requested)
    }
    .build()?)
}

fn work_directory(tmp_dir: &Path) -> Result<PathBuf, DynError> {
    let stamp = SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs();
    let path = tmp_dir.join(format!("pepcluster2-{}-{stamp}", std::process::id()));
    fs::create_dir_all(&path)?;
    Ok(path)
}

fn available_disk_bytes(path: &Path) -> Option<u128> {
    let output = Command::new("df").arg("-Pk").arg(path).output().ok()?;
    if !output.status.success() {
        return None;
    }
    let text = String::from_utf8(output.stdout).ok()?;
    let fields: Vec<&str> = text.lines().last()?.split_whitespace().collect();
    fields.get(3)?.parse::<u128>().ok().map(|kb| kb * 1024)
}

/// Sampled fraction of index hits that survive the sound anchor upper bound.
/// The non-prefilter path spills only bound-passing pairs, so the raw
/// occurrence upper bound overestimates temporary disk by roughly this factor.
fn sampled_bound_retention(
    nodes: &[Node],
    buckets: &[Vec<u32>],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    seed: TerminalSeed,
) -> f64 {
    const SAMPLE: usize = 2048;
    if nodes.len() < 2 {
        return 1.0;
    }
    let step = (nodes.len() / SAMPLE).max(1);
    let (retrieved, kept) = (0..nodes.len())
        .step_by(step)
        .collect::<Vec<usize>>()
        .par_iter()
        .map(|&id| {
            let candidates =
                retrieve_candidates(Some(id as u32), &nodes[id], buckets, table, seed);
            let passed = candidates
                .iter()
                .filter(|other| {
                    scorer.anchor_bound_passes(&nodes[id], &nodes[**other as usize])
                })
                .count() as u64;
            (candidates.len() as u64, passed)
        })
        .reduce(|| (0u64, 0u64), |a, b| (a.0 + b.0, a.1 + b.1));
    if retrieved == 0 {
        1.0
    } else {
        kept as f64 / retrieved as f64
    }
}

fn estimated_non_prefilter_disk_bytes(index: IndexStats, bound_retention: f64) -> u128 {
    let raw = index.candidate_occurrence_upper_bound.saturating_mul(16);
    ((raw as f64) * bound_retention.clamp(0.0, 1.0)) as u128
}

fn choose_prefilter(
    config: &Config,
    index: IndexStats,
    bound_retention: f64,
) -> (bool, u128, Option<u128>, String) {
    let estimated = estimated_non_prefilter_disk_bytes(index, bound_retention);
    let available = available_disk_bytes(&config.tmp_dir);
    if config.clustering_method == ClusteringMethod::Greedy {
        return (
            false,
            estimated,
            available,
            "not applicable to greedy clustering".into(),
        );
    }
    match config.prefilter {
        PrefilterChoice::Force => (true, estimated, available, "forced".into()),
        PrefilterChoice::Disable => (false, estimated, available, "disabled".into()),
        PrefilterChoice::Auto => {
            let active = available.is_some_and(|free| estimated > free.saturating_mul(80) / 100);
            let reason = available.map_or_else(
                || "automatic: available disk could not be measured; prefilter disabled".into(),
                |free| format!("automatic: estimated {estimated} bytes (sampled anchor-bound retention {bound_retention:.4}); available {free} bytes; 80% safety limit"),
            );
            (active, estimated, available, reason)
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn write_config(
    config: &Config,
    scorer: Option<&Scorer>,
    threads: usize,
    prefilter_active: bool,
    estimated_disk: u128,
    available_disk: Option<u128>,
    reason: &str,
) -> Result<(), DynError> {
    let mut w = BufWriter::new(File::create(config.output_dir.join("run_config.txt"))?);
    writeln!(w, "pepcluster2_version={VERSION}")?;
    writeln!(w, "input={}", config.input.display())?;
    writeln!(w, "output_dir={}", config.output_dir.display())?;
    writeln!(w, "tmp_dir={}", config.tmp_dir.display())?;
    writeln!(w, "mode={}", config.mode.name())?;
    writeln!(w, "clustering_method={}", config.clustering_method.name())?;
    writeln!(w, "greedy_selection={}", config.greedy_selection.name())?;
    writeln!(w, "threshold={}", config.threshold)?;
    writeln!(w, "threshold_explicit={}", config.threshold_explicit)?;
    writeln!(
        w,
        "alignment_similarity_threshold={}",
        config.alignment_threshold
    )?;
    writeln!(
        w,
        "anchor_combination_similarity_threshold={}",
        config.anchor_threshold
    )?;
    writeln!(w, "kmer_similarity_threshold={}", config.kmer_threshold)?;
    writeln!(
        w,
        "kmer_threshold_overridden={}",
        config.kmer_threshold_overridden
    )?;
    writeln!(
        w,
        "alignment_threshold_overridden={}",
        config.alignment_threshold_overridden
    )?;
    writeln!(
        w,
        "anchor_threshold_overridden={}",
        config.anchor_threshold_overridden
    )?;
    writeln!(w, "gap_open={}", config.gap_open)?;
    writeln!(w, "gap_extension={}", config.gap_extension)?;
    writeln!(w, "terminal_overhang_gap_open={}", config.terminal_gap_open)?;
    writeln!(
        w,
        "terminal_overhang_gap_extension={}",
        config.terminal_gap_extension
    )?;
    writeln!(
        w,
        "minimum_terminal_match_length={}",
        config.minimum_terminal_match_length
    )?;
    writeln!(w, "kmer_seed_threshold={}", config.kmer_seed_threshold)?;
    writeln!(w, "terminal_seed={}", config.terminal_seed.name())?;
    writeln!(w, "reassignment_margin={}", config.reassignment_margin)?;
    writeln!(
        w,
        "representative_order={}",
        config.representative_order.name()
    )?;
    writeln!(w, "prefilter_choice={:?}", config.prefilter)?;
    writeln!(w, "prefilter_active={prefilter_active}")?;
    writeln!(
        w,
        "full_sensitive_after_prefilter={}",
        config.full_sensitive_after_prefilter
    )?;
    writeln!(w, "prefilter_decision={reason}")?;
    writeln!(w, "estimated_non_prefilter_disk_bytes={estimated_disk}")?;
    writeln!(
        w,
        "available_tmp_disk_bytes={}",
        available_disk
            .map(|x| x.to_string())
            .unwrap_or_else(|| "unknown".into())
    )?;
    writeln!(w, "kmer_table={}", config.kmer_table.display())?;
    writeln!(w, "threads={threads}")?;
    writeln!(
        w,
        "scorer={}",
        scorer.map(Scorer::simd_name).unwrap_or("not_used")
    )?;
    writeln!(
        w,
        "iteration_cap={}",
        config
            .iteration_cap
            .map(|x| x.to_string())
            .unwrap_or_else(|| "unset".into())
    )?;
    writeln!(
        w,
        "merge_cap={}",
        config
            .merge_cap
            .map(|x| x.to_string())
            .unwrap_or_else(|| "all".into())
    )?;
    writeln!(w, "strict_merge={}", config.merge)?;
    writeln!(w, "strict_input={}", config.strict)?;
    writeln!(w, "candidate_buffer_mb={}", config.candidate_buffer_mb)?;
    writeln!(w, "max_memory_gb={}", config.max_memory_gb)?;
    writeln!(w, "merge_motifs={}", config.merge_motifs)?;
    writeln!(
        w,
        "motif_prior_concentration={}",
        config.motif.prior_concentration
    )?;
    writeln!(w, "motif_merge_threshold={}", config.motif.merge_threshold)?;
    writeln!(w, "motif_em={}", config.motif.em)?;
    writeln!(
        w,
        "motif_em_prior_concentration={}",
        config.motif.em_prior_concentration
    )?;
    writeln!(
        w,
        "motif_em_max_iterations={}",
        config.motif.em_max_iterations
    )?;
    writeln!(w, "motif_em_tolerance={}", config.motif.em_tolerance)?;
    writeln!(w, "compact_output={}", config.compact_output)?;
    writeln!(w, "write_cluster_fastas={}", config.write_cluster_fastas)?;
    writeln!(w, "write_scored_pairs={}", config.write_scored_pairs)?;
    w.flush()?;
    fs::write(
        config.output_dir.join("command.txt"),
        format!("{}\n", env::args().collect::<Vec<_>>().join(" ")),
    )?;
    Ok(())
}

#[derive(Default)]
struct EdgeReports {
    prefilter: Option<EdgeBuildStats>,
    sensitive: Option<EdgeBuildStats>,
    final_edges: u64,
    provisional_clusters: usize,
}

fn cluster_sizes(nodes: &[Node], clustering: &Clustering) -> (Vec<u64>, Vec<u64>) {
    let mut unique = vec![0u64; clustering.representatives.len()];
    let mut peptides = vec![0u64; clustering.representatives.len()];
    for (node, &cluster) in clustering.cluster_of.iter().enumerate() {
        unique[cluster as usize] += 1;
        peptides[cluster as usize] += nodes[node].frequency;
    }
    (unique, peptides)
}

#[allow(clippy::too_many_arguments)]
fn write_run_reports(
    config: &Config,
    fasta: FastaStats,
    nodes: &[Node],
    index: IndexStats,
    prefilter_active: bool,
    decision: &str,
    edges: &EdgeReports,
    greedy_stats: GreedyRunStats,
    clustering: &Clustering,
    iteration: IterationStats,
    motifs: Option<&MotifResult>,
    elapsed: f64,
    stages: &[(&str, f64)],
) -> Result<(), DynError> {
    let (unique_sizes, peptide_sizes) = cluster_sizes(nodes, clustering);
    let mut sorted = peptide_sizes.clone();
    sorted.sort_unstable();
    let singleton = peptide_sizes.iter().filter(|x| **x == 1).count();
    let median = if sorted.is_empty() {
        0.0
    } else if sorted.len() % 2 == 1 {
        sorted[sorted.len() / 2] as f64
    } else {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) as f64 / 2.0
    };
    let sensitive = edges.sensitive.as_ref();
    let pre = edges.prefilter.as_ref();
    let graph_computed = sensitive.map(|x| x.unique_candidate_pairs).unwrap_or(0)
        + pre.map(|x| x.unique_candidate_pairs).unwrap_or(0);
    let computed = if config.clustering_method == ClusteringMethod::Graph {
        graph_computed
    } else {
        greedy_stats.candidate_pairs_scored
    };
    let possible = nodes.len() as u128 * nodes.len().saturating_sub(1) as u128 / 2;
    let fraction = if possible == 0 {
        0.0
    } else {
        computed as f64 / possible as f64
    };
    // Cost decomposition. Candidate volume, cheap-bound rejections and
    // constrained-alignment evaluations are different costs and a single
    // "candidate pairs scored" figure hides which one dominates.
    let index_hits = if config.clustering_method == ClusteringMethod::Graph {
        sensitive.map(|x| x.seed_candidate_occurrences).unwrap_or(0)
            + pre.map(|x| x.seed_candidate_occurrences).unwrap_or(0)
    } else {
        greedy_stats.index_candidate_occurrences
    };
    let bound_rejected = if config.clustering_method == ClusteringMethod::Graph {
        sensitive.map(|x| x.anchor_bound_rejected).unwrap_or(0)
            + pre.map(|x| x.anchor_bound_rejected).unwrap_or(0)
    } else {
        greedy_stats.anchor_bound_rejected
    };
    // Counted globally, so it covers candidate scoring, representative updates,
    // merge validation and reassignment on every clustering path.
    let alignment_evaluations = scoring::alignment_evaluations();

    let mut summary = String::new();
    summary.push_str("PEPCLUSTER2 RUN SUMMARY\n=======================\n");
    summary.push_str(&format!(
        "Version: {VERSION}\nScoring mode: {}\nClustering method: {}\nGreedy selection: {}\nTerminal seed: {}\nRepresentative order: {}\n",
        config.mode.name(),
        config.clustering_method.name(),
        config.greedy_selection.name(),
        config.terminal_seed.name(),
        config.representative_order.name()
    ));
    summary.push_str(&format!("Input records: {}\nAccepted peptides: {}\nExcluded peptides: {}\nUnique peptide sequences: {}\n", fasta.records, fasta.accepted, fasta.skipped, nodes.len()));
    summary.push_str(&format!(
        "Occupied composite keys: {}\nLargest index bucket: {}\n",
        index.occupied_keys, index.largest_bucket
    ));
    summary.push_str(&format!(
        "Prefilter active: {prefilter_active}\nPrefilter decision: {decision}\n"
    ));
    if prefilter_active && !config.full_sensitive_after_prefilter {
        summary.push_str("WARNING: scoped graph completion is approximate and may differ from the non-prefilter graph.\n");
    }
    summary.push_str(&format!("Index candidate hits: {index_hits}\nRejected by anchor upper bound: {bound_rejected}\n"));
    summary.push_str(&format!("Candidate pairs scored: {computed}\nFraction of all unique-sequence pairs scored: {fraction:.8}\nConstrained-alignment evaluations: {alignment_evaluations}\n"));
    summary.push_str(&format!("Eligible graph edges: {}\nGreedy eligible assignments/retrievals: {}\nRepresentative-update pair scores: {}\nMerge pair scores: {}\n", edges.final_edges, greedy_stats.eligible_pairs, greedy_stats.representative_pair_scores, greedy_stats.merge_pair_scores));
    summary.push_str(&format!("Iterations: {}\nConverged: {}\nReassignment moves: {}\nRepresentative changes: {}\nValidated merges: {}\nValidation failures: {}\n", iteration.iterations, iteration.converged, iteration.reassignment_moves, iteration.representative_changes, iteration.merges, iteration.validation_failures));
    summary.push_str(&format!("Clusters: {}\nSingleton clusters: {}\nMean peptide cluster size: {:.3}\nMedian peptide cluster size: {:.3}\nLargest peptide cluster: {}\nLargest unique-sequence cluster: {}\nElapsed seconds: {:.6}\n", clustering.representatives.len(), singleton, fasta.accepted as f64 / clustering.representatives.len().max(1) as f64, median, peptide_sizes.iter().max().unwrap_or(&0), unique_sizes.iter().max().unwrap_or(&0), elapsed));
    if let Some(m) = motifs {
        summary.push_str(&format!(
            "\nMOTIF LAYER\nMerged motif groups: {}\nOccupied motifs: {}\nAccepted merges: {}\nEM iterations: {}\nEM converged: {}\nMotif prior concentration: {}\nMotif merge threshold: {}\nEM prior concentration: {}\n",
            m.merged_count,
            m.motif_count,
            m.merges,
            m.em_iterations,
            m.em_converged,
            config.motif.prior_concentration,
            config.motif.merge_threshold,
            config.motif.em_prior_concentration
        ));
        summary.push_str("The motif partition does not satisfy the representative-to-member invariant and is reported separately from the similarity clusters.\n");
    }
    summary.push_str("\nSTAGE TIMINGS\n");
    for (name, seconds) in stages {
        summary.push_str(&format!("{name}: {seconds:.6} s\n"));
    }
    fs::write(config.output_dir.join("run_summary.txt"), summary)?;

    let stage_json = stages
        .iter()
        .map(|(n, s)| format!("    {}: {:.6}", json_escape(n), s))
        .collect::<Vec<_>>()
        .join(",\n");
    let motif_json = motifs.map_or_else(
        || "  \"motif_layer\": false,\n".to_string(),
        |m| {
            format!(
                concat!(
                    "  \"motif_layer\": true,\n  \"motif_merged_groups\": {},\n  \"motif_occupied\": {},\n",
                    "  \"motif_merges\": {},\n  \"motif_em_iterations\": {},\n  \"motif_em_converged\": {},\n",
                    "  \"motif_prior_concentration\": {},\n  \"motif_merge_threshold\": {},\n  \"motif_em_prior_concentration\": {},\n"
                ),
                m.merged_count,
                m.motif_count,
                m.merges,
                m.em_iterations,
                m.em_converged,
                config.motif.prior_concentration,
                config.motif.merge_threshold,
                config.motif.em_prior_concentration
            )
        },
    );
    let json = format!(concat!(
        "{{\n  \"run_type\": \"cluster\",\n  \"pepcluster2_version\": {},\n  \"scoring_mode\": {},\n  \"clustering_method\": {},\n",
        "  \"greedy_selection\": {},\n  \"terminal_seed\": {},\n  \"representative_order\": {},\n  \"kmer_seed_threshold\": {},\n  \"reassignment_margin\": {},\n",
        "  \"input_records\": {},\n  \"accepted_records\": {},\n  \"skipped_records\": {},\n  \"unique_sequences\": {},\n",
        "  \"prefilter_active\": {},\n  \"prefilter_candidate_pairs\": {},\n  \"valid_prefilter_edges\": {},\n  \"sensitive_candidate_pairs\": {},\n",
        "  \"index_candidate_hits\": {},\n  \"anchor_bound_rejected\": {},\n  \"alignment_evaluations\": {},\n",
        "  \"candidate_pairs_computed\": {},\n  \"fraction_all_pairs_computed\": {:.12},\n  \"graph_edge_count\": {},\n",
        "  \"greedy_eligible_pairs\": {},\n  \"representative_pair_scores\": {},\n  \"merge_pair_scores\": {},\n",
        "  \"final_clusters\": {},\n  \"singleton_clusters\": {},\n  \"iterations\": {},\n  \"converged\": {},\n  \"reassignment_moves\": {},\n  \"representative_changes\": {},\n  \"strict_merges\": {},\n  \"validation_failures\": {},\n",
        "  \"largest_cluster_peptides\": {},\n{}  \"elapsed_seconds\": {:.6},\n  \"stage_seconds\": {{\n{}\n  }}\n}}\n"),
        json_escape(VERSION), json_escape(config.mode.name()), json_escape(config.clustering_method.name()), json_escape(config.greedy_selection.name()),
        json_escape(config.terminal_seed.name()), json_escape(config.representative_order.name()), config.kmer_seed_threshold, config.reassignment_margin,
        fasta.records, fasta.accepted, fasta.skipped, nodes.len(), prefilter_active,
        pre.map(|x| x.unique_candidate_pairs).unwrap_or(0), pre.map(|x| x.unique_valid_edges).unwrap_or(0), sensitive.map(|x| x.unique_candidate_pairs).unwrap_or(0),
        index_hits, bound_rejected, alignment_evaluations,
        computed, fraction, edges.final_edges, greedy_stats.eligible_pairs, greedy_stats.representative_pair_scores, greedy_stats.merge_pair_scores,
        clustering.representatives.len(), singleton, iteration.iterations, iteration.converged, iteration.reassignment_moves, iteration.representative_changes, iteration.merges, iteration.validation_failures,
        peptide_sizes.iter().max().unwrap_or(&0), motif_json, elapsed, stage_json);
    fs::write(config.output_dir.join("run_stats.json"), json)?;
    Ok(())
}

fn make_scorer(config: &Config) -> Result<Scorer, DynError> {
    Scorer::new(
        config.mode,
        config.threshold,
        config.alignment_threshold,
        config.anchor_threshold,
        config.kmer_threshold,
        config.gap_open,
        config.gap_extension,
        config.terminal_gap_open,
        config.terminal_gap_extension,
        config.minimum_terminal_match_length,
        config.simd_mode,
    )
    .map_err(|error| error.into())
}

fn run(config: Config) -> Result<(), DynError> {
    let total = Instant::now();
    fs::create_dir_all(&config.output_dir)?;
    fs::create_dir_all(&config.tmp_dir)?;
    let pool = build_pool(config.threads)?;
    let threads = pool.current_num_threads();
    if config.write_scored_pairs {
        pair_trace::enable(threads);
    }
    eprintln!(
        "PepCluster2 {VERSION}; mode={}; method={}; threads={threads}",
        config.mode.name(),
        config.clustering_method.name()
    );
    let mut stages = Vec::<(&str, f64)>::new();

    let started = Instant::now();
    let (nodes, fasta_stats) = load_nodes(&config.input, config.strict)?;
    if nodes.len() > u32::MAX as usize {
        return Err("more than 2^32-1 unique peptide sequences are unsupported".into());
    }
    stages.push(("load_fasta", started.elapsed().as_secs_f64()));
    eprintln!(
        "[1/6] {} accepted records, {} unique sequences, {} excluded",
        fasta_stats.accepted,
        nodes.len(),
        fasta_stats.skipped
    );

    let started = Instant::now();
    let table =
        KmerSimilarityTable::open_or_create(&config.kmer_table, config.kmer_seed_threshold)?;
    let buckets = build_exact_index(&nodes, config.terminal_seed);
    let (relations, index_stats) = build_similar_key_relations(&buckets, &table);
    stages.push(("build_kmer_index", started.elapsed().as_secs_f64()));
    let scorer = make_scorer(&config)?;
    let bound_retention = pool.install(|| {
        sampled_bound_retention(&nodes, &buckets, &table, &scorer, config.terminal_seed)
    });
    let (prefilter_active, estimated_disk, available_disk, decision) =
        choose_prefilter(&config, index_stats, bound_retention);
    eprintln!(
        "[2/6] {} occupied keys; anchor-bound retention {bound_retention:.4}; prefilter={} ({})",
        index_stats.occupied_keys, prefilter_active, decision
    );

    if config.index_only {
        write_config(
            &config,
            None,
            threads,
            prefilter_active,
            estimated_disk,
            available_disk,
            &decision,
        )?;
        let json = format!("{{\n  \"run_type\": \"index_only\",\n  \"pepcluster2_version\": {},\n  \"unique_sequences\": {},\n  \"occupied_composite_keys\": {},\n  \"candidate_occurrence_upper_bound\": {},\n  \"sampled_anchor_bound_retention\": {:.6},\n  \"estimated_non_prefilter_disk_bytes\": {},\n  \"prefilter_active\": {}\n}}\n", json_escape(VERSION), nodes.len(), index_stats.occupied_keys, index_stats.candidate_occurrence_upper_bound, bound_retention, estimated_disk, prefilter_active);
        fs::write(config.output_dir.join("run_stats.json"), json)?;
        return Ok(());
    }

    write_config(
        &config,
        Some(&scorer),
        threads,
        prefilter_active,
        estimated_disk,
        available_disk,
        &decision,
    )?;
    let mut edge_reports = EdgeReports::default();
    let mut greedy_stats = GreedyRunStats::default();
    let mut edge_file: Option<PathBuf> = None;

    let started = Instant::now();
    let (clustering, iteration_stats) = match config.clustering_method {
        ClusteringMethod::Graph => {
            let work_dir = work_directory(&config.tmp_dir)?;
            let memory_limit = (config.max_memory_gb * 1024.0_f64.powi(3)) as u64;
            let buffer_bytes = config.candidate_buffer_mb * 1024 * 1024;
            let (final_file, final_count, initial) = if prefilter_active {
                let pre = pool.install(|| {
                    generate_prefilter_edges(
                        &nodes,
                        &table,
                        &scorer,
                        config.terminal_seed,
                        &work_dir,
                        buffer_bytes,
                        config.keep_tmp,
                        threads,
                    )
                })?;
                let pre_graph = load_graph(
                    &pre.edge_file,
                    nodes.len(),
                    pre.unique_valid_edges,
                    memory_limit,
                )?;
                let provisional =
                    greedy_set_cover(&nodes, &pre_graph, config.representative_order);
                if config.write_edges {
                    write_provisional_clusters(
                        &config.output_dir.join("prefilter_provisional_clusters.tsv"),
                        &nodes,
                        &provisional,
                    )?;
                }
                edge_reports.provisional_clusters = provisional.representatives.len();
                let (reps, unassigned) = prefilter_masks(nodes.len(), &provisional);
                let scope = if config.full_sensitive_after_prefilter {
                    CandidateScope::All
                } else {
                    CandidateScope::PrefilterCompletion {
                        representatives: &reps,
                        unassigned: &unassigned,
                    }
                };
                let sensitive = pool.install(|| {
                    generate_edges(
                        &relations,
                        &buckets,
                        &nodes,
                        &scorer,
                        &work_dir,
                        buffer_bytes,
                        config.keep_tmp,
                        threads,
                        EdgeMode::Sensitive,
                        scope,
                        "sensitive",
                    )
                })?;
                let merged = work_dir.join("final_edges.bin");
                let count = merge_edge_files(&[&pre.edge_file, &sensitive.edge_file], &merged)?;
                let graph = load_graph(&merged, nodes.len(), count, memory_limit)?;
                let initial = greedy_set_cover(&nodes, &graph, config.representative_order);
                edge_reports.prefilter = Some(pre);
                edge_reports.sensitive = Some(sensitive);
                (merged, count, initial)
            } else {
                let sensitive = pool.install(|| {
                    generate_edges(
                        &relations,
                        &buckets,
                        &nodes,
                        &scorer,
                        &work_dir,
                        buffer_bytes,
                        config.keep_tmp,
                        threads,
                        EdgeMode::Sensitive,
                        CandidateScope::All,
                        "sensitive",
                    )
                })?;
                let graph = load_graph(
                    &sensitive.edge_file,
                    nodes.len(),
                    sensitive.unique_valid_edges,
                    memory_limit,
                )?;
                let initial = greedy_set_cover(&nodes, &graph, config.representative_order);
                let file = sensitive.edge_file.clone();
                let count = sensitive.unique_valid_edges;
                edge_reports.sensitive = Some(sensitive);
                (file, count, initial)
            };
            edge_reports.final_edges = final_count;
            let graph = load_graph(&final_file, nodes.len(), final_count, memory_limit)?;
            let result = pool.install(|| {
                iterate_to_convergence(
                    &nodes,
                    &scorer,
                    &graph,
                    initial,
                    config.iteration_cap,
                    config.merge,
                    config.merge_cap,
                    config.reassignment_margin_q(),
                )
            });
            let clustering = canonicalize(result.0);
            let (covered, total_nodes) = representative_coverage(&graph, &clustering);
            if covered != total_nodes {
                return Err(format!(
                    "final graph representative validation failed for {} sequences",
                    total_nodes - covered
                )
                .into());
            }
            if config.write_edges {
                write_edges(
                    &final_file,
                    &config.output_dir.join("edges.tsv"),
                    &nodes,
                    &scorer,
                )?;
            }
            edge_file = Some(final_file);
            if !config.keep_tmp {
                fs::remove_dir_all(&work_dir)?;
            }
            (clustering, result.1)
        }
        ClusteringMethod::Greedy => {
            let (initial, initial_stats) = pool.install(|| match config.greedy_selection {
                GreedySelection::KmerDegree => greedy::initial_clustering(
                    &nodes,
                    &buckets,
                    &table,
                    &scorer,
                    config.terminal_seed,
                    config.representative_order,
                ),
                GreedySelection::LazyExact => greedy::initial_clustering_lazy_exact(
                    &nodes,
                    &buckets,
                    &table,
                    &scorer,
                    config.terminal_seed,
                ),
            });
            greedy_stats = initial_stats;
            let result = pool.install(|| {
                greedy::iterate_to_convergence(
                    &nodes,
                    &table,
                    &scorer,
                    initial,
                    config.iteration_cap,
                    config.merge,
                    config.merge_cap,
                    config.terminal_seed,
                    config.reassignment_margin_q(),
                    &mut greedy_stats,
                )
            });
            let (covered, total_nodes) =
                greedy::representative_coverage(&nodes, &scorer, &result.0);
            if covered != total_nodes {
                return Err(format!(
                    "final greedy representative validation failed for {} sequences",
                    total_nodes - covered
                )
                .into());
            }
            result
        }
    };
    stages.push((
        "candidate_scoring_and_clustering",
        started.elapsed().as_secs_f64(),
    ));
    if config.write_scored_pairs {
        pair_trace::write(&config.output_dir.join("scored_pairs.bin"))?;
    }
    eprintln!(
        "[3/6] {} clusters after {} iterations; converged={}",
        clustering.representatives.len(),
        iteration_stats.iterations,
        iteration_stats.converged
    );

    let motifs = if config.merge_motifs {
        let started = Instant::now();
        let result = pool.install(|| motif::build_motifs(&nodes, &clustering, &config.motif));
        stages.push(("motif_merge", started.elapsed().as_secs_f64()));
        eprintln!(
            "[4/6] motif layer: {} clusters -> {} merged -> {} motifs (EM {} iterations, converged={})",
            clustering.representatives.len(),
            result.merged_count,
            result.motif_count,
            result.em_iterations,
            result.em_converged
        );
        Some(result)
    } else {
        None
    };

    let started = Instant::now();
    if let Some(result) = motifs.as_ref() {
        write_motif_outputs(&config.output_dir, &nodes, &clustering, result)?;
    }
    if config.compact_output {
        write_compact_outputs(&config.output_dir, &nodes, &scorer, &clustering)?;
    } else {
        write_cluster_outputs(
            &config.input,
            &config.output_dir,
            config.strict,
            &nodes,
            &scorer,
            &clustering,
            config.min_cluster_size,
            config.write_cluster_fastas,
        )?;
    }
    stages.push(("write_outputs", started.elapsed().as_secs_f64()));
    write_run_reports(
        &config,
        fasta_stats,
        &nodes,
        index_stats,
        prefilter_active,
        &decision,
        &edge_reports,
        greedy_stats,
        &clustering,
        iteration_stats,
        motifs.as_ref(),
        total.elapsed().as_secs_f64(),
        &stages,
    )?;
    let _ = edge_file;
    eprintln!(
        "[6/6] results: {} ({:.2}s)",
        config.output_dir.display(),
        total.elapsed().as_secs_f64()
    );
    Ok(())
}

fn main() {
    let result = match cli::parse() {
        Ok(Some(config)) => run(config),
        Ok(None) => return,
        Err(error) => Err(error.into()),
    };
    if let Err(error) = result {
        eprintln!("error: {error}");
        std::process::exit(2);
    }
}

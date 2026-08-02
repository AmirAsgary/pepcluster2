//! Exhaustive reference partitions for the PepCluster2 validation.
//!
//! Every peptide pair is scored exactly, which removes candidate search from the
//! comparison entirely. The complete edge set then feeds the identical
//! clustering procedure the tool runs, so the reference is the partition
//! PepCluster2 would produce if its candidate search were perfect. Comparing a
//! run against it isolates candidate-search loss.
//!
//! * `pipeline_clusters.tsv` — the complete clustering procedure.
//! * `reassign_only_clusters.tsv` — the same with merging disabled, which
//!   isolates the contribution of the merge stage.
//!
//! Both are emitted for the `coverage` and `intrinsic` representative orders so
//! that a run is always compared against a reference using the same selection
//! rule.

#[path = "../../../../../src/edge_store.rs"]
mod edge_store;
#[path = "../../../../../src/fasta.rs"]
mod fasta;
#[path = "../../../../../src/graph.rs"]
mod graph;
#[path = "../../../../../src/index.rs"]
mod index;
#[path = "../../../../../src/kmer.rs"]
mod kmer;
#[path = "../../../../../src/model.rs"]
mod model;
#[path = "../../../../../src/pair_trace.rs"]
mod pair_trace;
#[path = "../../../../../src/scoring.rs"]
mod scoring;

use fasta::load_nodes;
use graph::{
    canonicalize, greedy_set_cover, iterate_to_convergence, Clustering, RepresentativeOrder,
};
use model::{Edge, Graph, Neighbor, Node};
use rayon::prelude::*;
use scoring::{Scorer, ScoringMode, SimdMode};
use std::env;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

fn build_graph(node_count: usize, edges: &[Edge]) -> Graph {
    let mut degrees = vec![0usize; node_count];
    for edge in edges {
        degrees[edge.u as usize] += 1;
        degrees[edge.v as usize] += 1;
    }
    let mut offsets = vec![0usize; node_count + 1];
    for i in 0..node_count {
        offsets[i + 1] = offsets[i] + degrees[i];
    }
    let mut cursor = offsets[..node_count].to_vec();
    let mut neighbors = vec![Neighbor { node: 0, weight: 0 }; offsets[node_count]];
    for edge in edges {
        for (from, to) in [(edge.u, edge.v), (edge.v, edge.u)] {
            let position = cursor[from as usize];
            neighbors[position] = Neighbor {
                node: to,
                weight: edge.weight,
            };
            cursor[from as usize] += 1;
        }
    }
    for i in 0..node_count {
        neighbors[offsets[i]..offsets[i + 1]].sort_unstable_by_key(|x| x.node);
    }
    Graph {
        offsets,
        neighbors,
        edge_count: edges.len() as u64,
    }
}

fn write_edges(path: &Path, edges: &[Edge]) -> std::io::Result<()> {
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, File::create(path)?);
    writer.write_all(b"PC2TRUE1")?;
    writer.write_all(&(edges.len() as u64).to_le_bytes())?;
    for edge in edges {
        let packed = ((edge.u as u64) << 32) | edge.v as u64;
        writer.write_all(&packed.to_le_bytes())?;
    }
    writer.flush()
}

fn write_partition(
    path: &Path,
    nodes: &[Node],
    clustering: &Clustering,
    prefix: &str,
) -> std::io::Result<usize> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "sequence\tcluster_id\trepresentative_sequence")?;
    for (node, &cluster) in clustering.cluster_of.iter().enumerate() {
        writeln!(
            writer,
            "{}\t{}_{:06}\t{}",
            String::from_utf8_lossy(&nodes[node].sequence),
            prefix,
            cluster + 1,
            String::from_utf8_lossy(
                &nodes[clustering.representatives[cluster as usize] as usize].sequence
            )
        )?;
    }
    writer.flush()?;
    let mut sizes = vec![0usize; clustering.representatives.len()];
    for cluster in &clustering.cluster_of {
        sizes[*cluster as usize] += 1;
    }
    Ok(sizes.into_iter().filter(|size| *size == 1).count())
}

struct ReferenceReport {
    clusters: usize,
    singletons: usize,
}

/// Two references from one edge set: the complete clustering procedure, and the
/// same procedure with merging disabled so the merge stage can be attributed.
fn emit_reference(
    output: &Path,
    nodes: &[Node],
    scorer: &Scorer,
    graph: &Graph,
    order: RepresentativeOrder,
    suffix: &str,
    reassignment_margin: u16,
) -> Result<[ReferenceReport; 2], Box<dyn std::error::Error + Send + Sync>> {
    let mut reports = Vec::with_capacity(2);
    for (merge, name, prefix) in [
        (false, "reassign_only_clusters", "RA"),
        (true, "pipeline_clusters", "PL"),
    ] {
        let start = greedy_set_cover(nodes, graph, order);
        let (partition, _) = iterate_to_convergence(
            nodes,
            scorer,
            graph,
            start,
            None,
            merge,
            None,
            reassignment_margin,
        );
        let partition = canonicalize(partition);
        let singletons = write_partition(
            &output.join(format!("{name}{suffix}.tsv")),
            nodes,
            &partition,
            prefix,
        )?;
        reports.push(ReferenceReport {
            clusters: partition.representatives.len(),
            singletons,
        });
    }
    Ok(reports.try_into().unwrap_or_else(|_| unreachable!()))
}

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let args: Vec<String> = env::args().collect();
    if !(4..=5).contains(&args.len()) {
        return Err(
            "usage: exhaustive-reference INPUT.fasta OUTPUT_DIR THREADS [REASSIGNMENT_MARGIN_Q]"
                .into(),
        );
    }
    let input = PathBuf::from(&args[1]);
    let output = PathBuf::from(&args[2]);
    let threads: usize = args[3].parse()?;
    // Must match the clustering default, or the reference is not the partition
    // the tool would produce from a perfect search.
    let reassignment_margin: u16 = match args.get(4) {
        Some(value) => value.parse()?,
        None => 10,
    };
    fs::create_dir_all(&output)?;
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;

    let total_started = Instant::now();
    let (nodes, fasta_stats) = load_nodes(&input, true)?;
    let scorer = Scorer::new(
        ScoringMode::SeparateAlnAnchor,
        0.60,
        0.50,
        0.60,
        0.60,
        -4.0,
        -1.0,
        -2.0,
        -1.0,
        2,
        SimdMode::Auto,
    )?;
    let possible_pairs = nodes.len() as u64 * nodes.len().saturating_sub(1) as u64 / 2;

    let scoring_started = Instant::now();
    let mut edges: Vec<Edge> = (0..nodes.len())
        .into_par_iter()
        .fold(Vec::new, |mut local, i| {
            for j in i + 1..nodes.len() {
                if let Some(weight) = scorer.score_scalar(&nodes[i], &nodes[j]) {
                    local.push(Edge {
                        u: i as u32,
                        v: j as u32,
                        weight,
                    });
                }
            }
            local
        })
        .reduce(Vec::new, |mut left, mut right| {
            left.append(&mut right);
            left
        });
    edges.par_sort_unstable_by_key(|edge| (edge.u, edge.v));
    let scoring_seconds = scoring_started.elapsed().as_secs_f64();

    let graph = build_graph(nodes.len(), &edges);
    write_edges(&output.join("true_pairs.bin"), &edges)?;

    let clustering_started = Instant::now();
    let [reassign, pipeline] = emit_reference(
        &output,
        &nodes,
        &scorer,
        &graph,
        RepresentativeOrder::Coverage,
        "",
        reassignment_margin,
    )?;
    let [reassign_intrinsic, pipeline_intrinsic] = emit_reference(
        &output,
        &nodes,
        &scorer,
        &graph,
        RepresentativeOrder::Intrinsic,
        "_intrinsic",
        reassignment_margin,
    )?;
    let clustering_seconds = clustering_started.elapsed().as_secs_f64();

    let stats = format!(
        concat!(
            "{{\n",
            "  \"input_records\": {},\n",
            "  \"unique_sequences\": {},\n",
            "  \"all_possible_pairs\": {},\n",
            "  \"true_eligible_pairs\": {},\n",
            "  \"reassign_only_clusters\": {},\n",
            "  \"reassign_only_singleton_clusters\": {},\n",
            "  \"pipeline_clusters\": {},\n",
            "  \"pipeline_singleton_clusters\": {},\n",
            "  \"reassign_only_clusters_intrinsic\": {},\n",
            "  \"reassign_only_singleton_clusters_intrinsic\": {},\n",
            "  \"pipeline_clusters_intrinsic\": {},\n",
            "  \"pipeline_singleton_clusters_intrinsic\": {},\n",
            "  \"scoring_seconds\": {:.6},\n",
            "  \"clustering_seconds\": {:.6},\n",
            "  \"total_seconds\": {:.6}\n",
            "}}\n"
        ),
        fasta_stats.accepted,
        nodes.len(),
        possible_pairs,
        edges.len(),
        reassign.clusters,
        reassign.singletons,
        pipeline.clusters,
        pipeline.singletons,
        reassign_intrinsic.clusters,
        reassign_intrinsic.singletons,
        pipeline_intrinsic.clusters,
        pipeline_intrinsic.singletons,
        scoring_seconds,
        clustering_seconds,
        total_started.elapsed().as_secs_f64(),
    );
    fs::write(output.join("run_stats.json"), stats)?;
    Ok(())
}

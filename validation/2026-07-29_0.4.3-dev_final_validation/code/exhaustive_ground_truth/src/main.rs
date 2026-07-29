#[path = "../../../../../src/fasta.rs"]
mod fasta;
#[path = "../../../../../src/model.rs"]
mod model;
#[path = "../../../../../src/pair_trace.rs"]
mod pair_trace;
#[path = "../../../../../src/scoring.rs"]
mod scoring;

use fasta::load_nodes;
use model::{Edge, Graph, Neighbor, Node};
use rayon::prelude::*;
use scoring::{Scorer, ScoringMode, SimdMode};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::env;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Clone, Copy, Eq, PartialEq)]
struct HeapEntry {
    coverage: u32,
    weight_sum: u64,
    frequency: u64,
    node: u32,
}

impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.coverage
            .cmp(&other.coverage)
            .then(self.weight_sum.cmp(&other.weight_sum))
            .then(self.frequency.cmp(&other.frequency))
            .then_with(|| other.node.cmp(&self.node))
    }
}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

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

fn greedy_set_cover(nodes: &[Node], graph: &Graph) -> (Vec<u32>, Vec<u32>) {
    let n = nodes.len();
    let mut assigned = vec![false; n];
    let mut cluster_of = vec![u32::MAX; n];
    let mut coverage: Vec<u32> = (0..n)
        .map(|i| 1 + graph.neighbors(i as u32).len() as u32)
        .collect();
    let mut weight_sum: Vec<u64> = (0..n)
        .map(|i| {
            graph
                .neighbors(i as u32)
                .iter()
                .map(|x| x.weight as u64)
                .sum()
        })
        .collect();
    let mut heap = BinaryHeap::with_capacity(n);
    for node in 0..n {
        heap.push(HeapEntry {
            coverage: coverage[node],
            weight_sum: weight_sum[node],
            frequency: nodes[node].frequency,
            node: node as u32,
        });
    }
    let mut representatives = Vec::new();
    while let Some(entry) = heap.pop() {
        let rep = entry.node as usize;
        if assigned[rep] || coverage[rep] != entry.coverage || weight_sum[rep] != entry.weight_sum {
            continue;
        }
        let cluster = representatives.len() as u32;
        representatives.push(entry.node);
        let mut newly_assigned = vec![entry.node];
        newly_assigned.extend(
            graph
                .neighbors(entry.node)
                .iter()
                .filter(|x| !assigned[x.node as usize])
                .map(|x| x.node),
        );
        newly_assigned.sort_unstable();
        newly_assigned.dedup();
        for &node in &newly_assigned {
            assigned[node as usize] = true;
            cluster_of[node as usize] = cluster;
        }
        for &removed in &newly_assigned {
            for neighbor in graph.neighbors(removed) {
                let candidate = neighbor.node as usize;
                if !assigned[candidate] {
                    coverage[candidate] -= 1;
                    weight_sum[candidate] -= neighbor.weight as u64;
                    heap.push(HeapEntry {
                        coverage: coverage[candidate],
                        weight_sum: weight_sum[candidate],
                        frequency: nodes[candidate].frequency,
                        node: candidate as u32,
                    });
                }
            }
        }
    }
    (cluster_of, representatives)
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

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        return Err("usage: exhaustive-ground-truth INPUT.fasta OUTPUT_DIR THREADS".into());
    }
    let input = PathBuf::from(&args[1]);
    let output = PathBuf::from(&args[2]);
    let threads: usize = args[3].parse()?;
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
    let clustering_started = Instant::now();
    let (cluster_of, representatives) = greedy_set_cover(&nodes, &graph);
    let clustering_seconds = clustering_started.elapsed().as_secs_f64();
    write_edges(&output.join("true_pairs.bin"), &edges)?;

    let mut assignments = BufWriter::new(File::create(output.join("ground_truth_clusters.tsv"))?);
    writeln!(assignments, "sequence\tcluster_id\trepresentative_sequence")?;
    for (node, &cluster) in cluster_of.iter().enumerate() {
        writeln!(
            assignments,
            "{}\tGT_{:06}\t{}",
            String::from_utf8_lossy(&nodes[node].sequence),
            cluster + 1,
            String::from_utf8_lossy(&nodes[representatives[cluster as usize] as usize].sequence)
        )?;
    }
    assignments.flush()?;

    let singleton_clusters = {
        let mut sizes = vec![0usize; representatives.len()];
        for cluster in &cluster_of {
            sizes[*cluster as usize] += 1;
        }
        sizes.into_iter().filter(|size| *size == 1).count()
    };
    let stats = format!(
        concat!(
            "{{\n",
            "  \"input_records\": {},\n",
            "  \"unique_sequences\": {},\n",
            "  \"all_possible_pairs\": {},\n",
            "  \"true_eligible_pairs\": {},\n",
            "  \"ground_truth_clusters\": {},\n",
            "  \"singleton_clusters\": {},\n",
            "  \"scoring_seconds\": {:.6},\n",
            "  \"clustering_seconds\": {:.6},\n",
            "  \"total_seconds\": {:.6}\n",
            "}}\n"
        ),
        fasta_stats.accepted,
        nodes.len(),
        possible_pairs,
        edges.len(),
        representatives.len(),
        singleton_clusters,
        scoring_seconds,
        clustering_seconds,
        total_started.elapsed().as_secs_f64(),
    );
    fs::write(output.join("run_stats.json"), stats)?;
    Ok(())
}

use crate::fasta::DynError;
use crate::index::KeyRelation;
use crate::kmer::{KmerSimilarityTable, N_DIMERS};
use crate::model::{Edge, Node};
use crate::scoring::Scorer;
use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

const PAIR_DISK_BYTES: usize = 8;
const RELATIONS_PER_TASK: usize = 8192;
const SCORE_BATCH_PAIRS: usize = 1_048_576;

#[derive(Default)]
struct CandidateTaskResult {
    seed_candidate_occurrences: u64,
    candidate_occurrences: u64,
    chunk_paths: Vec<PathBuf>,
}

#[derive(Clone, Debug)]
pub struct EdgeBuildStats {
    pub seed_candidate_occurrences: u64,
    pub candidate_occurrences: u64,
    pub unique_candidate_pairs: u64,
    pub unique_valid_edges: u64,
    pub candidate_chunk_count: usize,
    pub edge_file: PathBuf,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EdgeMode {
    Sensitive,
    Prefilter,
}

#[derive(Clone, Copy)]
pub enum CandidateScope<'a> {
    All,
    PrefilterCompletion {
        representatives: &'a [bool],
        unassigned: &'a [bool],
    },
}

impl CandidateScope<'_> {
    #[inline]
    fn permits(self, a: u32, b: u32) -> bool {
        match self {
            CandidateScope::All => true,
            CandidateScope::PrefilterCompletion {
                representatives,
                unassigned,
            } => {
                representatives[a as usize]
                    || representatives[b as usize]
                    // Unassigned peptides must be compared with every
                    // provisional member, not only provisional
                    // representatives. Otherwise a valid future
                    // representative can never recruit them.
                    || unassigned[a as usize]
                    || unassigned[b as usize]
            }
        }
    }
}

fn write_pair(writer: &mut impl Write, pair: (u32, u32)) -> io::Result<()> {
    writer.write_all(&pair.0.to_le_bytes())?;
    writer.write_all(&pair.1.to_le_bytes())
}

fn read_pair(reader: &mut impl Read) -> io::Result<Option<(u32, u32)>> {
    let mut first = [0u8; 4];
    let n = reader.read(&mut first)?;
    if n == 0 {
        return Ok(None);
    }
    if n != 4 {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "truncated pair record",
        ));
    }
    let mut second = [0u8; 4];
    reader.read_exact(&mut second)?;
    Ok(Some((
        u32::from_le_bytes(first),
        u32::from_le_bytes(second),
    )))
}

fn write_edge(writer: &mut impl Write, edge: Edge) -> io::Result<()> {
    writer.write_all(&edge.u.to_le_bytes())?;
    writer.write_all(&edge.v.to_le_bytes())?;
    writer.write_all(&edge.weight.to_le_bytes())
}

pub fn read_edge(reader: &mut impl Read) -> io::Result<Option<Edge>> {
    let Some((u, v)) = read_pair(reader)? else {
        return Ok(None);
    };
    let mut weight = [0u8; 2];
    reader.read_exact(&mut weight)?;
    Ok(Some(Edge {
        u,
        v,
        weight: u16::from_le_bytes(weight),
    }))
}

fn spill_pairs(
    pairs: &mut Vec<(u32, u32)>,
    work_dir: &Path,
    file_counter: &AtomicUsize,
    label: &str,
) -> Result<Option<PathBuf>, DynError> {
    if pairs.is_empty() {
        return Ok(None);
    }
    pairs.sort_unstable();
    pairs.dedup();
    let id = file_counter.fetch_add(1, Ordering::Relaxed);
    let path = work_dir.join(format!("{label}_candidate_chunk_{id:06}.bin"));
    let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(&path)?);
    for &pair in pairs.iter() {
        write_pair(&mut writer, pair)?;
    }
    writer.flush()?;
    pairs.clear();
    Ok(Some(path))
}

fn generate_candidate_task(
    relations: &[KeyRelation],
    buckets: &[Vec<u32>],
    work_dir: &Path,
    file_counter: &AtomicUsize,
    pair_buffer_limit: usize,
    nodes: &[Node],
    mode: EdgeMode,
    scope: CandidateScope<'_>,
    label: &str,
) -> Result<CandidateTaskResult, DynError> {
    let mut result = CandidateTaskResult::default();
    let mut pairs = Vec::<(u32, u32)>::with_capacity(pair_buffer_limit.min(2_000_000));
    for relation in relations {
        let first = &buckets[relation.first as usize];
        let second = &buckets[relation.second as usize];
        if relation.first == relation.second {
            for i in 0..first.len() {
                for j in (i + 1)..first.len() {
                    let pair = (first[i], first[j]);
                    result.seed_candidate_occurrences += 1;
                    if !scope.permits(pair.0, pair.1)
                        || (mode == EdgeMode::Prefilter
                            && Scorer::distinct_shared_anchor_types(
                                &nodes[pair.0 as usize],
                                &nodes[pair.1 as usize],
                            ) < 2)
                    {
                        continue;
                    }
                    pairs.push(pair);
                    result.candidate_occurrences += 1;
                    if pairs.len() >= pair_buffer_limit {
                        if let Some(path) = spill_pairs(&mut pairs, work_dir, file_counter, label)?
                        {
                            result.chunk_paths.push(path);
                        }
                    }
                }
            }
        } else {
            for &a in first {
                for &b in second {
                    if a == b {
                        continue;
                    }
                    let pair = (a.min(b), a.max(b));
                    result.seed_candidate_occurrences += 1;
                    if !scope.permits(pair.0, pair.1)
                        || (mode == EdgeMode::Prefilter
                            && Scorer::distinct_shared_anchor_types(
                                &nodes[pair.0 as usize],
                                &nodes[pair.1 as usize],
                            ) < 2)
                    {
                        continue;
                    }
                    pairs.push(pair);
                    result.candidate_occurrences += 1;
                    if pairs.len() >= pair_buffer_limit {
                        if let Some(path) = spill_pairs(&mut pairs, work_dir, file_counter, label)?
                        {
                            result.chunk_paths.push(path);
                        }
                    }
                }
            }
        }
    }
    if let Some(path) = spill_pairs(&mut pairs, work_dir, file_counter, label)? {
        result.chunk_paths.push(path);
    }
    Ok(result)
}

fn merge_candidate_chunks(
    chunk_paths: &[PathBuf],
    output_path: &Path,
    remove_chunks: bool,
) -> Result<u64, DynError> {
    let mut readers: Vec<BufReader<File>> = chunk_paths
        .iter()
        .map(|path| File::open(path).map(|file| BufReader::with_capacity(1024 * 1024, file)))
        .collect::<Result<_, _>>()?;
    let mut heap = BinaryHeap::<Reverse<(u32, u32, usize)>>::new();
    for (source, reader) in readers.iter_mut().enumerate() {
        if let Some((u, v)) = read_pair(reader)? {
            heap.push(Reverse((u, v, source)));
        }
    }
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, File::create(output_path)?);
    let mut previous: Option<(u32, u32)> = None;
    let mut unique = 0u64;
    while let Some(Reverse((u, v, source))) = heap.pop() {
        let pair = (u, v);
        if previous != Some(pair) {
            write_pair(&mut writer, pair)?;
            previous = Some(pair);
            unique += 1;
        }
        if let Some((next_u, next_v)) = read_pair(&mut readers[source])? {
            heap.push(Reverse((next_u, next_v, source)));
        }
    }
    writer.flush()?;
    drop(readers);
    if remove_chunks {
        for path in chunk_paths {
            fs::remove_file(path)?;
        }
    }
    Ok(unique)
}

fn score_unique_candidates(
    candidate_file: &Path,
    edge_file: &Path,
    nodes: &[Node],
    scorer: &Scorer,
    mode: EdgeMode,
) -> Result<u64, DynError> {
    let mut reader = BufReader::with_capacity(8 * 1024 * 1024, File::open(candidate_file)?);
    let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, File::create(edge_file)?);
    let mut pairs = Vec::<(u32, u32)>::with_capacity(SCORE_BATCH_PAIRS);
    let mut valid = 0u64;
    loop {
        pairs.clear();
        while pairs.len() < SCORE_BATCH_PAIRS {
            match read_pair(&mut reader)? {
                Some(pair) => pairs.push(pair),
                None => break,
            }
        }
        if pairs.is_empty() {
            break;
        }
        let scored: Vec<Vec<Edge>> = pairs
            .par_chunks(8192)
            .map(|chunk| {
                let mut edges = Vec::new();
                scorer.filter_pairs(nodes, chunk, mode == EdgeMode::Prefilter, &mut edges);
                edges
            })
            .collect();
        for edges in scored {
            valid += edges.len() as u64;
            for edge in edges {
                write_edge(&mut writer, edge)?;
            }
        }
        if pairs.len() < SCORE_BATCH_PAIRS {
            break;
        }
    }
    writer.flush()?;
    Ok(valid)
}

#[inline]
fn terminal_seed_hit(a: &Node, b: &Node, table: &KmerSimilarityTable) -> bool {
    let af = [
        a.codes[0] as usize * 20 + a.codes[1] as usize,
        a.codes[1] as usize * 20 + a.codes[2] as usize,
    ];
    let bf = [
        b.codes[0] as usize * 20 + b.codes[1] as usize,
        b.codes[1] as usize * 20 + b.codes[2] as usize,
    ];
    let ae = [
        a.codes[3] as usize * 20 + a.codes[4] as usize,
        a.codes[4] as usize * 20 + a.codes[5] as usize,
    ];
    let be = [
        b.codes[3] as usize * 20 + b.codes[4] as usize,
        b.codes[4] as usize * 20 + b.codes[5] as usize,
    ];
    af.iter()
        .any(|x| bf.iter().any(|y| table.are_neighbours(*x, *y)))
        && ae
            .iter()
            .any(|x| be.iter().any(|y| table.are_neighbours(*x, *y)))
}

fn generate_anchor_bucket_task(
    anchor_buckets: &[Vec<u32>],
    nodes: &[Node],
    table: &KmerSimilarityTable,
    work_dir: &Path,
    file_counter: &AtomicUsize,
    pair_buffer_limit: usize,
    label: &str,
) -> Result<CandidateTaskResult, DynError> {
    let mut result = CandidateTaskResult::default();
    let mut pairs = Vec::<(u32, u32)>::with_capacity(pair_buffer_limit.min(2_000_000));
    for bucket in anchor_buckets {
        for i in 0..bucket.len() {
            for j in (i + 1)..bucket.len() {
                result.seed_candidate_occurrences += 1;
                let pair = (bucket[i].min(bucket[j]), bucket[i].max(bucket[j]));
                if !terminal_seed_hit(&nodes[pair.0 as usize], &nodes[pair.1 as usize], table) {
                    continue;
                }
                result.candidate_occurrences += 1;
                pairs.push(pair);
                if pairs.len() >= pair_buffer_limit {
                    if let Some(path) = spill_pairs(&mut pairs, work_dir, file_counter, label)? {
                        result.chunk_paths.push(path);
                    }
                }
            }
        }
    }
    if let Some(path) = spill_pairs(&mut pairs, work_dir, file_counter, label)? {
        result.chunk_paths.push(path);
    }
    Ok(result)
}

/// High-confidence prefilter candidate generation starts from pairs of two
/// distinct exact anchor values. This avoids expanding every permissive k-mer
/// seed occurrence before applying the exact-anchor requirement.
pub fn generate_prefilter_edges(
    nodes: &[Node],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    work_dir: &Path,
    total_pair_buffer_bytes: usize,
    keep_tmp: bool,
    threads: usize,
) -> Result<EdgeBuildStats, DynError> {
    fs::create_dir_all(work_dir)?;
    let mut anchor_buckets = vec![Vec::<u32>::new(); N_DIMERS * N_DIMERS];
    for (node_id, node) in nodes.iter().enumerate() {
        let types = Scorer::distinct_anchor_types(node);
        for i in 0..types.len() {
            for j in (i + 1)..types.len() {
                anchor_buckets[types[i] as usize * N_DIMERS + types[j] as usize]
                    .push(node_id as u32);
            }
        }
    }
    let file_counter = AtomicUsize::new(0);
    let per_worker_bytes = total_pair_buffer_bytes / threads.max(1);
    let pair_buffer_limit = (per_worker_bytes / PAIR_DISK_BYTES).max(4096);
    let task_results: Vec<Result<CandidateTaskResult, DynError>> = anchor_buckets
        .par_chunks(512)
        .map(|chunk| {
            generate_anchor_bucket_task(
                chunk,
                nodes,
                table,
                work_dir,
                &file_counter,
                pair_buffer_limit,
                "prefilter",
            )
        })
        .collect();
    let mut seed_candidate_occurrences = 0u64;
    let mut candidate_occurrences = 0u64;
    let mut chunk_paths = Vec::new();
    for task in task_results {
        let task = task?;
        seed_candidate_occurrences += task.seed_candidate_occurrences;
        candidate_occurrences += task.candidate_occurrences;
        chunk_paths.extend(task.chunk_paths);
    }
    chunk_paths.sort();
    let candidate_file = work_dir.join("prefilter_candidate_pairs.bin");
    let unique_candidate_pairs = merge_candidate_chunks(&chunk_paths, &candidate_file, !keep_tmp)?;
    let edge_file = work_dir.join("prefilter_edges.bin");
    let unique_valid_edges = score_unique_candidates(
        &candidate_file,
        &edge_file,
        nodes,
        scorer,
        EdgeMode::Prefilter,
    )?;
    if !keep_tmp {
        fs::remove_file(&candidate_file)?;
    }
    Ok(EdgeBuildStats {
        seed_candidate_occurrences,
        candidate_occurrences,
        unique_candidate_pairs,
        unique_valid_edges,
        candidate_chunk_count: chunk_paths.len(),
        edge_file,
    })
}

pub fn generate_edges(
    relations: &[KeyRelation],
    buckets: &[Vec<u32>],
    nodes: &[Node],
    scorer: &Scorer,
    work_dir: &Path,
    total_pair_buffer_bytes: usize,
    keep_tmp: bool,
    threads: usize,
    mode: EdgeMode,
    scope: CandidateScope<'_>,
    label: &str,
) -> Result<EdgeBuildStats, DynError> {
    fs::create_dir_all(work_dir)?;
    let file_counter = AtomicUsize::new(0);
    let per_worker_bytes = total_pair_buffer_bytes / threads.max(1);
    let pair_buffer_limit = (per_worker_bytes / PAIR_DISK_BYTES).max(4096);
    let task_results: Vec<Result<CandidateTaskResult, DynError>> = relations
        .par_chunks(RELATIONS_PER_TASK)
        .map(|chunk| {
            generate_candidate_task(
                chunk,
                buckets,
                work_dir,
                &file_counter,
                pair_buffer_limit,
                nodes,
                mode,
                scope,
                label,
            )
        })
        .collect();
    let mut seed_candidate_occurrences = 0u64;
    let mut candidate_occurrences = 0u64;
    let mut chunk_paths = Vec::new();
    for task in task_results {
        let task = task?;
        seed_candidate_occurrences += task.seed_candidate_occurrences;
        candidate_occurrences += task.candidate_occurrences;
        chunk_paths.extend(task.chunk_paths);
    }
    chunk_paths.sort();
    let candidate_file = work_dir.join(format!("{label}_candidate_pairs.bin"));
    let unique_candidate_pairs = merge_candidate_chunks(&chunk_paths, &candidate_file, !keep_tmp)?;
    let edge_file = work_dir.join(format!("{label}_edges.bin"));
    let unique_valid_edges =
        score_unique_candidates(&candidate_file, &edge_file, nodes, scorer, mode)?;
    if !keep_tmp {
        fs::remove_file(&candidate_file)?;
    }
    Ok(EdgeBuildStats {
        seed_candidate_occurrences,
        candidate_occurrences,
        unique_candidate_pairs,
        unique_valid_edges,
        candidate_chunk_count: chunk_paths.len(),
        edge_file,
    })
}

/// Merge sorted edge streams and retain the largest weight for a duplicate pair.
pub fn merge_edge_files(inputs: &[&Path], output: &Path) -> Result<u64, DynError> {
    let mut readers: Vec<BufReader<File>> = inputs
        .iter()
        .map(|path| File::open(path).map(|f| BufReader::with_capacity(4 * 1024 * 1024, f)))
        .collect::<Result<_, _>>()?;
    let mut heap = BinaryHeap::<Reverse<(u32, u32, u16, usize)>>::new();
    for (source, reader) in readers.iter_mut().enumerate() {
        if let Some(edge) = read_edge(reader)? {
            heap.push(Reverse((edge.u, edge.v, edge.weight, source)));
        }
    }
    let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, File::create(output)?);
    let mut pending: Option<Edge> = None;
    let mut count = 0u64;
    while let Some(Reverse((u, v, weight, source))) = heap.pop() {
        if let Some(mut edge) = pending {
            if edge.u == u && edge.v == v {
                edge.weight = edge.weight.max(weight);
                pending = Some(edge);
            } else {
                write_edge(&mut writer, edge)?;
                count += 1;
                pending = Some(Edge { u, v, weight });
            }
        } else {
            pending = Some(Edge { u, v, weight });
        }
        if let Some(edge) = read_edge(&mut readers[source])? {
            heap.push(Reverse((edge.u, edge.v, edge.weight, source)));
        }
    }
    if let Some(edge) = pending {
        write_edge(&mut writer, edge)?;
        count += 1;
    }
    writer.flush()?;
    Ok(count)
}

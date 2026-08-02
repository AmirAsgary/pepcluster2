use crate::fasta::DynError;
use crate::index::{terminal_dimers, KeyRelation, TerminalSeed};
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
    anchor_bound_rejected: u64,
    candidate_occurrences: u64,
    chunk_paths: Vec<PathBuf>,
}

#[derive(Clone, Debug)]
pub struct EdgeBuildStats {
    /// Index hits before any rejection, counted with multiplicity.
    pub seed_candidate_occurrences: u64,
    /// Index hits discarded by the sound anchor upper bound.
    pub anchor_bound_rejected: u64,
    /// Index hits retained for spilling, counted with multiplicity.
    pub candidate_occurrences: u64,
    /// Distinct candidate pairs written to the candidate file and exactly scored.
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

/// Fill `buffer` completely, or report a clean end of stream.
///
/// `Read::read` may return fewer bytes than requested even when the stream has
/// more, and `BufReader` does exactly that whenever a record straddles its
/// internal buffer boundary. Records must therefore be assembled in a loop; a
/// single `read` is only correct for streams that never return a short read,
/// which files on a networked filesystem are not.
fn fill_record(reader: &mut impl Read, buffer: &mut [u8]) -> io::Result<bool> {
    let mut filled = 0usize;
    while filled < buffer.len() {
        match reader.read(&mut buffer[filled..]) {
            Ok(0) => {
                if filled == 0 {
                    return Ok(false);
                }
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated record",
                ));
            }
            Ok(read) => filled += read,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    Ok(true)
}

fn read_pair(reader: &mut impl Read) -> io::Result<Option<(u32, u32)>> {
    let mut record = [0u8; PAIR_DISK_BYTES];
    if !fill_record(reader, &mut record)? {
        return Ok(None);
    }
    Ok(Some((
        u32::from_le_bytes(record[..4].try_into().expect("four bytes")),
        u32::from_le_bytes(record[4..8].try_into().expect("four bytes")),
    )))
}

fn write_edge(writer: &mut impl Write, edge: Edge) -> io::Result<()> {
    writer.write_all(&edge.u.to_le_bytes())?;
    writer.write_all(&edge.v.to_le_bytes())?;
    writer.write_all(&edge.weight.to_le_bytes())
}

pub fn read_edge(reader: &mut impl Read) -> io::Result<Option<Edge>> {
    let mut record = [0u8; PAIR_DISK_BYTES + 2];
    if !fill_record(reader, &mut record)? {
        return Ok(None);
    }
    Ok(Some(Edge {
        u: u32::from_le_bytes(record[..4].try_into().expect("four bytes")),
        v: u32::from_le_bytes(record[4..8].try_into().expect("four bytes")),
        weight: u16::from_le_bytes(record[8..10].try_into().expect("two bytes")),
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

#[allow(clippy::too_many_arguments)]
fn generate_candidate_task(
    relations: &[KeyRelation],
    buckets: &[Vec<u32>],
    work_dir: &Path,
    file_counter: &AtomicUsize,
    pair_buffer_limit: usize,
    nodes: &[Node],
    scorer: &Scorer,
    mode: EdgeMode,
    scope: CandidateScope<'_>,
    label: &str,
) -> Result<CandidateTaskResult, DynError> {
    let mut result = CandidateTaskResult::default();
    let mut pairs = Vec::<(u32, u32)>::with_capacity(pair_buffer_limit.min(2_000_000));
    for relation in relations {
        let first = &buckets[relation.first as usize];
        let second = &buckets[relation.second as usize];
        let same_key = relation.first == relation.second;
        for (position, &a) in first.iter().enumerate() {
            let partners = if same_key {
                &second[position + 1..]
            } else {
                &second[..]
            };
            for &b in partners {
                if a == b {
                    continue;
                }
                let pair = (a.min(b), a.max(b));
                result.seed_candidate_occurrences += 1;
                if !scope.permits(pair.0, pair.1) {
                    continue;
                }
                let left = &nodes[pair.0 as usize];
                let right = &nodes[pair.1 as usize];
                if mode == EdgeMode::Prefilter
                    && Scorer::distinct_shared_anchor_types(left, right) < 2
                {
                    continue;
                }
                // Lossless rejection: the relaxed anchor assignment is an upper
                // bound on the exact anchor-combination similarity, so a pair
                // failing it cannot be accepted and never reaches disk.
                if !scorer.anchor_bound_passes(left, right) {
                    result.anchor_bound_rejected += 1;
                    continue;
                }
                pairs.push(pair);
                result.candidate_occurrences += 1;
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

/// Direct form of the index rule: at least one neighbouring front column pair
/// and at least one neighbouring end column pair.
#[inline]
fn terminal_seed_hit(
    a: &Node,
    b: &Node,
    table: &KmerSimilarityTable,
    seed: TerminalSeed,
) -> bool {
    let mut matched = [false; 2];
    for (slot, offset) in [0usize, 3usize].into_iter().enumerate() {
        let (left, n_left) = terminal_dimers(a, offset, seed);
        let (right, n_right) = terminal_dimers(b, offset, seed);
        matched[slot] = left[..n_left].iter().any(|x| {
            right[..n_right]
                .iter()
                .any(|y| table.are_neighbours(*x as usize, *y as usize))
        });
        if !matched[slot] {
            return false;
        }
    }
    matched[0] && matched[1]
}

#[allow(clippy::too_many_arguments)]
fn generate_anchor_bucket_task(
    anchor_buckets: &[Vec<u32>],
    nodes: &[Node],
    table: &KmerSimilarityTable,
    seed: TerminalSeed,
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
                if !terminal_seed_hit(
                    &nodes[pair.0 as usize],
                    &nodes[pair.1 as usize],
                    table,
                    seed,
                ) {
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
#[allow(clippy::too_many_arguments)]
pub fn generate_prefilter_edges(
    nodes: &[Node],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    seed: TerminalSeed,
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
                seed,
                work_dir,
                &file_counter,
                pair_buffer_limit,
                "prefilter",
            )
        })
        .collect();
    let mut seed_candidate_occurrences = 0u64;
    let mut anchor_bound_rejected = 0u64;
    let mut candidate_occurrences = 0u64;
    let mut chunk_paths = Vec::new();
    for task in task_results {
        let task = task?;
        seed_candidate_occurrences += task.seed_candidate_occurrences;
        anchor_bound_rejected += task.anchor_bound_rejected;
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
        anchor_bound_rejected,
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
                scorer,
                mode,
                scope,
                label,
            )
        })
        .collect();
    let mut seed_candidate_occurrences = 0u64;
    let mut anchor_bound_rejected = 0u64;
    let mut candidate_occurrences = 0u64;
    let mut chunk_paths = Vec::new();
    for task in task_results {
        let task = task?;
        seed_candidate_occurrences += task.seed_candidate_occurrences;
        anchor_bound_rejected += task.anchor_bound_rejected;
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
        anchor_bound_rejected,
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

#[cfg(test)]
mod tests {
    use super::*;

    /// A reader that never returns more than `chunk` bytes per call, which is
    /// what `BufReader` does when a record straddles its buffer boundary and
    /// what a networked filesystem does on large files.
    struct ShortReader<'a> {
        data: &'a [u8],
        chunk: usize,
    }

    impl Read for ShortReader<'_> {
        fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
            let take = self.data.len().min(out.len()).min(self.chunk);
            out[..take].copy_from_slice(&self.data[..take]);
            self.data = &self.data[take..];
            Ok(take)
        }
    }

    #[test]
    fn records_survive_short_reads() {
        let pairs = [(1u32, 2u32), (3, 4), (5, 6), (7, 8)];
        let mut bytes = Vec::new();
        for pair in pairs {
            write_pair(&mut bytes, pair).unwrap();
        }
        for chunk in 1..=PAIR_DISK_BYTES + 1 {
            let mut reader = ShortReader {
                data: &bytes,
                chunk,
            };
            let mut seen = Vec::new();
            while let Some(pair) = read_pair(&mut reader).unwrap() {
                seen.push(pair);
            }
            assert_eq!(seen, pairs, "pair stream corrupted at chunk size {chunk}");
        }

        let edges = [
            Edge {
                u: 1,
                v: 2,
                weight: 900,
            },
            Edge {
                u: 3,
                v: 4,
                weight: 10,
            },
        ];
        let mut bytes = Vec::new();
        for edge in edges {
            write_edge(&mut bytes, edge).unwrap();
        }
        for chunk in 1..=PAIR_DISK_BYTES + 3 {
            let mut reader = ShortReader {
                data: &bytes,
                chunk,
            };
            let mut seen = Vec::new();
            while let Some(edge) = read_edge(&mut reader).unwrap() {
                seen.push(edge);
            }
            assert_eq!(seen, edges, "edge stream corrupted at chunk size {chunk}");
        }
    }

    #[test]
    fn a_partial_trailing_record_is_an_error() {
        let mut bytes = Vec::new();
        write_pair(&mut bytes, (1, 2)).unwrap();
        bytes.push(0);
        let mut reader = ShortReader {
            data: &bytes,
            chunk: 3,
        };
        assert_eq!(read_pair(&mut reader).unwrap(), Some((1, 2)));
        assert!(read_pair(&mut reader).is_err());
    }
}

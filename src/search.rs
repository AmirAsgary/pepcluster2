//! Easy search: score every query peptide against a database of peptides and
//! report the distribution of its terminal k-mer similarity, rather than
//! clustering either set.
//!
//! # Why this can be exact and still cheap
//!
//! Terminal k-mer similarity (Section 6 of `ALGORITHM.md`) reads only the first
//! three and last three residues of each peptide. Two database peptides sharing
//! a terminal hexamer are therefore *indistinguishable* to this quantity: every
//! query scores identically against both. The database is collapsed onto its
//! distinct terminal hexamers, carrying a count per hexamer, and every statistic
//! is recovered from the weighted distribution. Nothing is approximated - the
//! reported numbers are exactly those of a full query-by-database comparison.
//!
//! It decomposes further. The similarity is
//!
//! ```text
//! front = clamp(max(0, sum of three normalised scores) / 3, 0, 1000)
//! end   = same over the last three residues
//! sim   = (front + end + 1) / 2
//! ```
//!
//! `front` depends only on the two peptides' first three residues and `end` only
//! on their last three, so for one query it is enough to score its front 3-mer
//! against each *distinct* database front 3-mer, and likewise at the end. Each
//! distinct hexamer is then two table lookups and an addition. There are at most
//! 20^3 = 8000 distinct 3-mers per terminus, so the per-query cost is bounded by
//! the number of distinct hexamers rather than by the number of peptides.
//!
//! Because similarity is quantised to thousandths, the distribution over the
//! database is a histogram of 1001 bins and quantiles are read straight off it.
//!
//! # Two populations, reported side by side
//!
//! The k-mer seed of Section 5 decides which pairs the clustering would ever
//! retrieve: at least one neighbouring front column pair and at least one
//! neighbouring end column pair. Restricting the statistics to seed hits is not
//! the same question as computing them over the whole database, and the answers
//! differ sharply - a median over the whole database is usually zero and says
//! little, while a median over hits describes the neighbourhood a query actually
//! sits in. Both are reported, along with `n_hits`, so neither has to be guessed
//! from the other.
//!
//! The seed condition factorises exactly. Retrieval enumerates the full cross
//! product of front and end keys (`composite_keys` in `index.rs`), so a database
//! peptide is retrieved precisely when one of its front keys neighbours one of
//! the query's and one of its end keys neighbours one of the query's. That test
//! also reads only the terminal residues, so it too is answered from the
//! collapsed database.

use crate::fasta::{load_nodes, DynError};
use crate::index::{build_exact_index, retrieve_candidates, terminal_dimers, TerminalSeed};
use crate::kmer::KmerSimilarityTable;
use crate::model::Node;
use crate::scoring::{normalized_residue_scores, Scorer};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// Similarity is quantised to thousandths, so the distribution needs 1001 bins.
const BINS: usize = 1001;
#[inline]
fn trimer_index(codes: &[u8]) -> usize {
    codes[0] as usize * 400 + codes[1] as usize * 20 + codes[2] as usize
}

/// One terminus, byte-for-byte the rule in `Scorer::aligned_three_mer_mean`.
///
/// The `.max(0)` applies to the SUM, per terminus, before dividing. A terminus
/// of three mismatching residues therefore contributes zero rather than a
/// negative amount and cannot drag the other terminus down, which is why the
/// whole quantity is not the mean of the six positions.
#[inline]
fn terminus_score(residue: &[i32; 400], left: usize, right: usize) -> i32 {
    let (a0, a1, a2) = (left / 400, (left / 20) % 20, left % 20);
    let (b0, b1, b2) = (right / 400, (right / 20) % 20, right % 20);
    let sum = residue[a0 * 20 + b0] + residue[a1 * 20 + b1] + residue[a2 * 20 + b2];
    sum.max(0).div_euclid(3).clamp(0, 1000)
}

/// The database collapsed onto its distinct terminal hexamers.
struct Database {
    /// Distinct front 3-mers, as indices into the 8000-way 3-mer space.
    front_trimers: Vec<usize>,
    /// Distinct end 3-mers.
    end_trimers: Vec<usize>,
    /// One entry per distinct hexamer: (index into `front_trimers`,
    /// index into `end_trimers`, how many database peptides carry it).
    hexamers: Vec<(u32, u32, u64)>,
    /// Front seed keys of each distinct front 3-mer, for the seed test.
    front_keys: Vec<([u32; 3], usize)>,
    end_keys: Vec<([u32; 3], usize)>,
    /// Database peptides represented, counting duplicates.
    total: u64,
    /// Members of each distinct hexamer as database node indices, in compressed
    /// row form: `members[offsets[h]..offsets[h + 1]]`. Only built for
    /// `--store-search`, which is the one thing that needs to name the target
    /// rather than only count it. Two u32 arrays, so 12M peptides cost ~96 MB.
    hexamer_offsets: Vec<u32>,
    hexamer_members: Vec<u32>,
}

fn build_database(nodes: &[Node], seed: TerminalSeed, with_members: bool) -> Database {
    let mut front_map: HashMap<usize, u32> = HashMap::new();
    let mut end_map: HashMap<usize, u32> = HashMap::new();
    let mut front_trimers: Vec<usize> = Vec::new();
    let mut end_trimers: Vec<usize> = Vec::new();
    let mut hexamer_map: HashMap<(u32, u32), usize> = HashMap::new();
    let mut hexamers: Vec<(u32, u32, u64)> = Vec::new();
    let mut total: u64 = 0;
    let mut node_hexamer: Vec<u32> = if with_members {
        Vec::with_capacity(nodes.len())
    } else {
        Vec::new()
    };

    for node in nodes {
        let front = trimer_index(&node.codes[0..3]);
        let end = trimer_index(&node.codes[3..6]);
        let front_id = *front_map.entry(front).or_insert_with(|| {
            front_trimers.push(front);
            (front_trimers.len() - 1) as u32
        });
        let end_id = *end_map.entry(end).or_insert_with(|| {
            end_trimers.push(end);
            (end_trimers.len() - 1) as u32
        });
        // A node already carries the frequency of the duplicate sequences it
        // collapsed, so the count is a peptide count, not a node count.
        let weight = node.frequency;
        total += weight;
        let hexamer_id = match hexamer_map.entry((front_id, end_id)) {
            std::collections::hash_map::Entry::Occupied(entry) => {
                let id = *entry.get();
                hexamers[id].2 += weight;
                id
            }
            std::collections::hash_map::Entry::Vacant(entry) => {
                let id = hexamers.len();
                entry.insert(id);
                hexamers.push((front_id, end_id, weight));
                id
            }
        };
        if with_members {
            node_hexamer.push(hexamer_id as u32);
        }
    }

    // Seed keys depend only on a terminus, so they are derived once per distinct
    // 3-mer rather than once per peptide.
    let keys_of = |trimers: &[usize], offset: usize| -> Vec<([u32; 3], usize)> {
        trimers
            .iter()
            .map(|&trimer| {
                let mut codes = [0u8; 6];
                codes[offset] = (trimer / 400) as u8;
                codes[offset + 1] = ((trimer / 20) % 20) as u8;
                codes[offset + 2] = (trimer % 20) as u8;
                let probe = Node {
                    sequence: Vec::new(),
                    sequence_codes: Vec::new(),
                    anchor: [0; 6],
                    codes,
                    combination_mask: 0,
                    frequency: 0,
                };
                terminal_dimers(&probe, offset, seed)
            })
            .collect()
    };
    let front_keys = keys_of(&front_trimers, 0);
    let end_keys = keys_of(&end_trimers, 3);

    // Compressed row form, built by counting then placing, so no per-hexamer
    // Vec is ever allocated.
    let (hexamer_offsets, hexamer_members) = if with_members {
        let mut offsets = vec![0u32; hexamers.len() + 1];
        for &hexamer in &node_hexamer {
            offsets[hexamer as usize + 1] += 1;
        }
        for index in 1..offsets.len() {
            offsets[index] += offsets[index - 1];
        }
        let mut cursor = offsets.clone();
        let mut members = vec![0u32; node_hexamer.len()];
        for (node_index, &hexamer) in node_hexamer.iter().enumerate() {
            let slot = &mut cursor[hexamer as usize];
            members[*slot as usize] = node_index as u32;
            *slot += 1;
        }
        (offsets, members)
    } else {
        (Vec::new(), Vec::new())
    };

    Database {
        front_trimers,
        end_trimers,
        hexamers,
        front_keys,
        end_keys,
        total,
        hexamer_offsets,
        hexamer_members,
    }
}

/// Value at rank `k` of the distribution described by `histogram`.
fn value_at_rank(histogram: &[u64; BINS], mut k: u64) -> f64 {
    for (value, &count) in histogram.iter().enumerate() {
        if count > k {
            return value as f64;
        }
        k -= count;
    }
    (BINS - 1) as f64
}

/// Quantile of a weighted histogram, matching numpy's default linear
/// interpolation so the output can be checked against the Python reference.
fn quantile(histogram: &[u64; BINS], total: u64, q: f64) -> f64 {
    if total == 0 {
        return f64::NAN;
    }
    let position = (total - 1) as f64 * q;
    let lower = position.floor();
    let upper = position.ceil();
    let low = value_at_rank(histogram, lower as u64);
    if (upper - lower).abs() < f64::EPSILON {
        return low;
    }
    let high = value_at_rank(histogram, upper as u64);
    low + (position - lower) * (high - low)
}

fn maximum(histogram: &[u64; BINS]) -> f64 {
    for value in (0..BINS).rev() {
        if histogram[value] > 0 {
            return value as f64;
        }
    }
    f64::NAN
}

struct Row {
    max_all: f64,
    q99_all: f64,
    q98_all: f64,
    q97_all: f64,
    median_all: f64,
    n_hits: u64,
    max_hit: f64,
    q99_hit: f64,
    q98_hit: f64,
    q97_hit: f64,
    median_hit: f64,
}

/// Everything about one query that the collapsed pass needs: its terminus
/// scores against every distinct database terminus, and its seed reachability.
/// Computed once per query and reused for every hexamer.
type QueryTables = (Vec<i32>, Vec<i32>, Vec<bool>, Vec<bool>);

fn query_tables(
    query: &Node,
    database: &Database,
    residue: &[i32; 400],
    table: &KmerSimilarityTable,
    seed: TerminalSeed,
) -> QueryTables {
    let query_front = trimer_index(&query.codes[0..3]);
    let query_end = trimer_index(&query.codes[3..6]);
    let front_scores: Vec<i32> = database
        .front_trimers
        .iter()
        .map(|&trimer| terminus_score(residue, query_front, trimer))
        .collect();
    let end_scores: Vec<i32> = database
        .end_trimers
        .iter()
        .map(|&trimer| terminus_score(residue, query_end, trimer))
        .collect();

    // Seed reachability, one boolean per distinct terminus. A dimer of the
    // database terminus must neighbour a dimer of the query's.
    let (query_front_keys, n_query_front) = terminal_dimers(query, 0, seed);
    let (query_end_keys, n_query_end) = terminal_dimers(query, 3, seed);
    let reachable = |keys: &([u32; 3], usize), query_keys: &[u32; 3], n: usize| -> bool {
        for &key in &keys.0[..keys.1] {
            for &other in &query_keys[..n] {
                if table.are_neighbours(key as usize, other as usize) {
                    return true;
                }
            }
        }
        false
    };
    let front_reachable: Vec<bool> = database
        .front_keys
        .iter()
        .map(|keys| reachable(keys, &query_front_keys, n_query_front))
        .collect();
    let end_reachable: Vec<bool> = database
        .end_keys
        .iter()
        .map(|keys| reachable(keys, &query_end_keys, n_query_end))
        .collect();
    (front_scores, end_scores, front_reachable, end_reachable)
}

fn search_one(
    query: &Node,
    database: &Database,
    residue: &[i32; 400],
    table: &KmerSimilarityTable,
    seed: TerminalSeed,
) -> Row {
    let (front_scores, end_scores, front_reachable, end_reachable) =
        query_tables(query, database, residue, table, seed);

    let mut all = [0u64; BINS];
    let mut hits = [0u64; BINS];
    let mut n_hits: u64 = 0;
    for &(front_id, end_id, count) in &database.hexamers {
        let front = front_scores[front_id as usize];
        let end = end_scores[end_id as usize];
        let similarity = ((front + end + 1) / 2) as usize;
        all[similarity] += count;
        if front_reachable[front_id as usize] && end_reachable[end_id as usize] {
            hits[similarity] += count;
            n_hits += count;
        }
    }

    let scale = 1000.0;
    Row {
        max_all: maximum(&all) / scale,
        q99_all: quantile(&all, database.total, 0.99) / scale,
        q98_all: quantile(&all, database.total, 0.98) / scale,
        q97_all: quantile(&all, database.total, 0.97) / scale,
        median_all: quantile(&all, database.total, 0.5) / scale,
        n_hits,
        max_hit: maximum(&hits) / scale,
        q99_hit: quantile(&hits, n_hits, 0.99) / scale,
        q98_hit: quantile(&hits, n_hits, 0.98) / scale,
        q97_hit: quantile(&hits, n_hits, 0.97) / scale,
        median_hit: quantile(&hits, n_hits, 0.5) / scale,
    }
}

/// One stored hit: which database peptide, how similar, and whether the k-mer
/// seed would have retrieved it.
pub struct Hit {
    pub target: u32,
    pub similarity: u16,
    pub seed_hit: bool,
}

/// Highest bin such that at least `wanted` targets score at or above it.
///
/// Similarity is quantised to 1001 values, so the top-N cut is read straight off
/// a counting histogram instead of sorting. Sorting the ~10 M distinct hexamers
/// of a 12 M database, per query, would dominate everything else in the run.
fn cutoff_for(counts: &[u64; BINS], wanted: u64, floor: u16) -> u16 {
    let mut seen = 0u64;
    for value in (floor as usize..BINS).rev() {
        seen += counts[value];
        if seen >= wanted {
            return value as u16;
        }
    }
    floor
}

/// Top hits of one query, via the collapsed database.
///
/// Two passes over the distinct hexamers rather than one pass plus a sort: the
/// first counts targets per similarity to find the cut, the second collects the
/// ones above it. Each pass is two array lookups per hexamer, so this stays far
/// cheaper than sorting and needs no per-query allocation proportional to the
/// database.
fn top_hits_collapsed(
    query: &Node,
    database: &Database,
    residue: &[i32; 400],
    scorer: &Scorer,
    table: &KmerSimilarityTable,
    seed: TerminalSeed,
    wanted: usize,
    floor: u16,
) -> Vec<Hit> {
    let (front_scores, end_scores, front_reachable, end_reachable) =
        query_tables(query, database, residue, table, seed);

    let mut counts = [0u64; BINS];
    for (index, &(front_id, end_id, _)) in database.hexamers.iter().enumerate() {
        let similarity = scorer.combine_termini(
            front_scores[front_id as usize] as u16,
            end_scores[end_id as usize] as u16,
        );
        if similarity >= floor {
            let start = database.hexamer_offsets[index] as usize;
            let end = database.hexamer_offsets[index + 1] as usize;
            counts[similarity as usize] += (end - start) as u64;
        }
    }
    let cut = cutoff_for(&counts, wanted as u64, floor);

    let mut hits: Vec<Hit> = Vec::with_capacity(wanted.min(1 << 16));
    for (index, &(front_id, end_id, _)) in database.hexamers.iter().enumerate() {
        let similarity = scorer.combine_termini(
            front_scores[front_id as usize] as u16,
            end_scores[end_id as usize] as u16,
        );
        if similarity < cut {
            continue;
        }
        let seed_hit =
            front_reachable[front_id as usize] && end_reachable[end_id as usize];
        let start = database.hexamer_offsets[index] as usize;
        let end = database.hexamer_offsets[index + 1] as usize;
        for &member in &database.hexamer_members[start..end] {
            hits.push(Hit { target: member, similarity, seed_hit });
        }
    }
    // Ties at the cut can overshoot `wanted`; rank first, then truncate, so the
    // kept rows are the best ones rather than whichever came first.
    hits.sort_unstable_by(|a, b| b.similarity.cmp(&a.similarity).then(a.target.cmp(&b.target)));
    hits.truncate(wanted);
    hits
}

/// Top hits of one query, scoring every pair. The only route in an alignment
/// mode, where similarity reads the whole peptide and the collapse does not
/// apply. A bounded heap keeps the cost at O(D) plus a log N factor paid only
/// when a candidate beats the current worst kept hit - negligible beside a
/// constrained alignment at roughly 5.5 microseconds.
fn top_hits_exhaustive(
    query: &Node,
    database_nodes: &[Node],
    buckets: &[Vec<u32>],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    seed: TerminalSeed,
    wanted: usize,
    floor: u16,
    hits_only: bool,
) -> Vec<Hit> {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;

    let retrieved = retrieve_candidates(None, query, buckets, table, seed);
    let mut is_hit = vec![false; database_nodes.len()];
    for id in &retrieved {
        is_hit[*id as usize] = true;
    }

    let mut heap: BinaryHeap<Reverse<(u16, u32)>> = BinaryHeap::with_capacity(wanted + 1);
    let mut consider = |index: usize, heap: &mut BinaryHeap<Reverse<(u16, u32)>>| {
        let similarity = reported_similarity(scorer, query, &database_nodes[index]) as u16;
        if similarity < floor {
            return;
        }
        if heap.len() < wanted {
            heap.push(Reverse((similarity, index as u32)));
        } else if let Some(Reverse((worst, _))) = heap.peek() {
            if similarity > *worst {
                heap.pop();
                heap.push(Reverse((similarity, index as u32)));
            }
        }
    };
    if hits_only {
        for id in &retrieved {
            consider(*id as usize, &mut heap);
        }
    } else {
        for index in 0..database_nodes.len() {
            consider(index, &mut heap);
        }
    }

    let mut hits: Vec<Hit> = heap
        .into_iter()
        .map(|Reverse((similarity, target))| Hit {
            target,
            similarity,
            seed_hit: is_hit[target as usize],
        })
        .collect();
    hits.sort_unstable_by(|a, b| b.similarity.cmp(&a.similarity).then(a.target.cmp(&b.target)));
    hits
}

/// The similarity this run reports, taken from the scoring mode.
///
/// Modes that align report the constrained full-alignment similarity of
/// Section 8, which reads every residue including the core. Modes that do not
/// align report the terminal k-mer similarity of Section 6, which reads only the
/// six terminal residues. Either way it is one component, not the mode's
/// eligibility rule: the anchor-combination component is not folded in, so a
/// high value does not by itself mean the pair would be clustered together.
#[inline]
fn reported_similarity(scorer: &Scorer, a: &Node, b: &Node) -> usize {
    let scores = scorer.scores(a, b);
    if scorer.mode().uses_alignment() {
        scores.alignment as usize
    } else {
        scores.terminal_kmer as usize
    }
}

/// One row computed pair by pair: every query against every database peptide,
/// retrieval by the production index, similarity by the production scorer.
///
/// This is the only way to compute an alignment-mode row, because the
/// constrained alignment reads the whole peptide and so cannot be answered from
/// a database collapsed onto terminal hexamers. In a k-mer mode it is instead the
/// reference `--verify-search` checks the collapsed path against.
///
/// When `hits_only` is set, pairs the seed does not retrieve are never scored and
/// the whole-database columns are left undefined. That is the difference between
/// roughly one tenth of the work and all of it.
fn exhaustive_row(
    query: &Node,
    database_nodes: &[Node],
    buckets: &[Vec<u32>],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    seed: TerminalSeed,
    hits_only: bool,
) -> Row {
    let mut all = [0u64; BINS];
    let mut hits = [0u64; BINS];
    let mut total: u64 = 0;
    let mut n_hits: u64 = 0;

    // `None` because the query is not a member of the database index: a query
    // that also occurs in the database must count as its own hit.
    let retrieved = retrieve_candidates(None, query, buckets, table, seed);

    if hits_only {
        for id in &retrieved {
            let node = &database_nodes[*id as usize];
            let similarity = reported_similarity(scorer, query, node);
            hits[similarity] += node.frequency;
            n_hits += node.frequency;
        }
    } else {
        // One pass over the database, with a membership mask, so no pair is
        // aligned twice. At roughly 3.8 microseconds per constrained alignment
        // that second pass would be a straight doubling of the run.
        let mut is_hit = vec![false; database_nodes.len()];
        for id in &retrieved {
            is_hit[*id as usize] = true;
        }
        for (index, node) in database_nodes.iter().enumerate() {
            let similarity = reported_similarity(scorer, query, node);
            all[similarity] += node.frequency;
            total += node.frequency;
            if is_hit[index] {
                hits[similarity] += node.frequency;
                n_hits += node.frequency;
            }
        }
    }
    let scale = 1000.0;
    Row {
        max_all: maximum(&all) / scale,
        q99_all: quantile(&all, total, 0.99) / scale,
        q98_all: quantile(&all, total, 0.98) / scale,
        q97_all: quantile(&all, total, 0.97) / scale,
        median_all: quantile(&all, total, 0.5) / scale,
        n_hits,
        max_hit: maximum(&hits) / scale,
        q99_hit: quantile(&hits, n_hits, 0.99) / scale,
        q98_hit: quantile(&hits, n_hits, 0.98) / scale,
        q97_hit: quantile(&hits, n_hits, 0.97) / scale,
        median_hit: quantile(&hits, n_hits, 0.5) / scale,
    }
}

/// Every field of two rows agrees, treating NaN as equal to NaN so that empty
/// hit sets compare cleanly.
fn rows_agree(fast: &Row, slow: &Row) -> bool {
    let same = |a: f64, b: f64| (a.is_nan() && b.is_nan()) || a == b;
    fast.n_hits == slow.n_hits
        && same(fast.max_all, slow.max_all)
        && same(fast.q99_all, slow.q99_all)
        && same(fast.q98_all, slow.q98_all)
        && same(fast.q97_all, slow.q97_all)
        && same(fast.median_all, slow.median_all)
        && same(fast.max_hit, slow.max_hit)
        && same(fast.q99_hit, slow.q99_hit)
        && same(fast.q98_hit, slow.q98_hit)
        && same(fast.q97_hit, slow.q97_hit)
        && same(fast.median_hit, slow.median_hit)
}

fn field(value: f64) -> String {
    if value.is_nan() {
        String::new()
    } else {
        format!("{value:.6}")
    }
}

/// Settings for `--store-search`.
pub struct StoreParams {
    pub path: PathBuf,
    /// Rows kept per query, best first.
    pub max_hits: usize,
    /// Quantised similarity floor; hits below it are never written.
    pub floor: u16,
}

/// Queries per parallel batch. Bounds peak memory at `chunk * max_hits` rows
/// rather than `queries * max_hits`: at 500k queries and 100 hits the whole
/// result would be 50 M rows held at once, which is gigabytes for no reason.
const STORE_CHUNK: usize = 8192;

/// Write the per-query hit table.
///
/// Rows carry the target peptide inline, so the file needs no dictionary and no
/// reader. That costs roughly five times what a target index would, which at
/// this scale is the right trade: 1.5 GB against 0.3 GB for 500k queries at 100
/// hits, and the text form goes straight into pandas or awk.
#[allow(clippy::too_many_arguments)]
fn write_hits(
    query_nodes: &[Node],
    database_nodes: &[Node],
    database: &Database,
    buckets: &[Vec<u32>],
    residue: &[i32; 400],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    seed: TerminalSeed,
    collapsible: bool,
    hits_only: bool,
    params: &StoreParams,
) -> Result<(), DynError> {
    let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, File::create(&params.path)?);
    writeln!(writer, "query\ttarget\tsimilarity\trank\tseed_hit")?;
    let mut written: u64 = 0;
    for chunk in query_nodes.chunks(STORE_CHUNK) {
        let batch: Vec<Vec<Hit>> = chunk
            .par_iter()
            .map(|query| {
                if collapsible {
                    top_hits_collapsed(
                        query, database, residue, scorer, table, seed,
                        params.max_hits, params.floor,
                    )
                } else {
                    top_hits_exhaustive(
                        query, database_nodes, buckets, table, scorer, seed,
                        params.max_hits, params.floor, hits_only,
                    )
                }
            })
            .collect();
        for (query, hits) in chunk.iter().zip(&batch) {
            let name = String::from_utf8_lossy(&query.sequence);
            for (rank, hit) in hits.iter().enumerate() {
                writeln!(
                    writer,
                    "{}\t{}\t{:.3}\t{}\t{}",
                    name,
                    String::from_utf8_lossy(&database_nodes[hit.target as usize].sequence),
                    hit.similarity as f64 / 1000.0,
                    rank + 1,
                    if hit.seed_hit { 1 } else { 0 }
                )?;
                written += 1;
            }
        }
    }
    writer.flush()?;
    eprintln!(
        "[store] wrote {written} hit rows to {}",
        params.path.display()
    );
    Ok(())
}

/// Run the search and write one CSV row per query peptide.
#[allow(clippy::too_many_arguments)]
pub fn run(
    database_path: &Path,
    query_path: &Path,
    output: &Path,
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    seed: TerminalSeed,
    strict: bool,
    verify: bool,
    hits_only: bool,
    store: Option<StoreParams>,
) -> Result<(usize, usize, usize), DynError> {
    let (database_nodes, database_stats) = load_nodes(database_path, strict)?;
    let (query_nodes, query_stats) = load_nodes(query_path, strict)?;
    if database_nodes.is_empty() || query_nodes.is_empty() {
        return Err("both the database and the query set must contain usable peptides".into());
    }
    let database = build_database(&database_nodes, seed, store.is_some());
    let residue = normalized_residue_scores();

    // The collapsed path is available only when the reported similarity reads
    // nothing but the terminal residues. A constrained alignment reads the whole
    // peptide, core included, so two peptides sharing a terminal hexamer are no
    // longer interchangeable and the collapse would be wrong, not merely slow.
    let collapsible = !scorer.mode().uses_alignment() && !hits_only;
    let buckets = if collapsible && !verify {
        Vec::new()
    } else {
        build_exact_index(&database_nodes, seed)
    };

    let rows: Vec<Row> = if collapsible {
        query_nodes
            .par_iter()
            .map(|query| search_one(query, &database, &residue, table, seed))
            .collect()
    } else {
        query_nodes
            .par_iter()
            .map(|query| {
                exhaustive_row(
                    query,
                    &database_nodes,
                    &buckets,
                    table,
                    scorer,
                    seed,
                    hits_only,
                )
            })
            .collect()
    };

    if verify {
        if !collapsible {
            return Err("--verify-search checks the collapsed path against an \
                        all-against-all reference; there is no collapsed path to check \
                        in an alignment mode or under --search-hits-only"
                .into());
        }
        let disagreements: usize = query_nodes
            .par_iter()
            .zip(&rows)
            .filter(|(query, fast)| {
                let slow = exhaustive_row(
                    query,
                    &database_nodes,
                    &buckets,
                    table,
                    scorer,
                    seed,
                    false,
                );
                !rows_agree(fast, &slow)
            })
            .count();
        if disagreements > 0 {
            return Err(format!(
                "--verify-search: {disagreements} of {} query rows disagree with the \
                 all-against-all reference",
                query_nodes.len()
            )
            .into());
        }
        eprintln!(
            "[verify] all {} query rows reproduce the all-against-all reference exactly",
            query_nodes.len()
        );
    }

    if let Some(params) = store {
        write_hits(
            &query_nodes,
            &database_nodes,
            &database,
            &buckets,
            &residue,
            table,
            scorer,
            seed,
            collapsible,
            hits_only,
            &params,
        )?;
    }

    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, File::create(output)?);
    writeln!(
        writer,
        "peptide,n_hits,max,quantile99,quantile98,quantile97,median,\
         hits_max,hits_quantile99,hits_quantile98,hits_quantile97,hits_median"
    )?;
    for (query, row) in query_nodes.iter().zip(&rows) {
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{},{}",
            String::from_utf8_lossy(&query.sequence),
            row.n_hits,
            field(row.max_all),
            field(row.q99_all),
            field(row.q98_all),
            field(row.q97_all),
            field(row.median_all),
            field(row.max_hit),
            field(row.q99_hit),
            field(row.q98_hit),
            field(row.q97_hit),
            field(row.median_hit)
        )?;
    }
    writer.flush()?;
    let _ = (database_stats, query_stats);
    Ok((
        query_nodes.len(),
        database_nodes.len(),
        database.hexamers.len(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hist(values: &[(usize, u64)]) -> [u64; BINS] {
        let mut h = [0u64; BINS];
        for &(value, count) in values {
            h[value] += count;
        }
        h
    }

    #[test]
    fn quantile_matches_linear_interpolation() {
        // values 0,1,2,3 each once: numpy.quantile([0,1,2,3], 0.5) == 1.5
        let h = hist(&[(0, 1), (1, 1), (2, 1), (3, 1)]);
        assert!((quantile(&h, 4, 0.5) - 1.5).abs() < 1e-9);
        // q=0 and q=1 hit the ends exactly
        assert!((quantile(&h, 4, 0.0) - 0.0).abs() < 1e-9);
        assert!((quantile(&h, 4, 1.0) - 3.0).abs() < 1e-9);
    }

    #[test]
    fn weights_behave_like_repeats() {
        // [5,5,5,9] as weights must match the same multiset spelled out
        let weighted = hist(&[(5, 3), (9, 1)]);
        let spelled = hist(&[(5, 1), (5, 1), (5, 1), (9, 1)]);
        for q in [0.0, 0.25, 0.5, 0.97, 0.99, 1.0] {
            assert!((quantile(&weighted, 4, q) - quantile(&spelled, 4, q)).abs() < 1e-9);
        }
    }

    #[test]
    fn terminus_floors_a_negative_sum_at_zero() {
        let residue = normalized_residue_scores();
        // W against three residues it scores poorly with: the sum is negative,
        // so the terminus must contribute exactly zero rather than a negative.
        let w = 17usize; // index of W in ARNDCQEGHILKMFPSTWYV
        let d = 3usize; // D
        let left = w * 400 + w * 20 + w;
        let right = d * 400 + d * 20 + d;
        assert_eq!(terminus_score(&residue, left, right), 0);
    }

    #[test]
    fn cutoff_keeps_at_least_the_requested_number() {
        // 5 targets at 900, 10 at 800, 100 at 100.
        let mut counts = [0u64; BINS];
        counts[900] = 5;
        counts[800] = 10;
        counts[100] = 100;
        // 3 wanted: 900 alone already supplies 5, so the cut is 900.
        assert_eq!(cutoff_for(&counts, 3, 0), 900);
        // 5 wanted: still 900, exactly enough.
        assert_eq!(cutoff_for(&counts, 5, 0), 900);
        // 6 wanted: must drop to 800, which brings the total to 15.
        assert_eq!(cutoff_for(&counts, 6, 0), 800);
        // More than exist: the cut falls to the floor and everything is kept.
        assert_eq!(cutoff_for(&counts, 10_000, 0), 0);
        // The floor is never crossed, however many are wanted.
        assert_eq!(cutoff_for(&counts, 10_000, 500), 500);
    }

    #[test]
    fn identical_termini_score_one() {
        let residue = normalized_residue_scores();
        let trimer = 1usize * 400 + 2 * 20 + 3;
        assert_eq!(terminus_score(&residue, trimer, trimer), 1000);
    }
}

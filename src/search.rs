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
use std::path::Path;

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
}

fn build_database(nodes: &[Node], seed: TerminalSeed) -> Database {
    let mut front_map: HashMap<usize, u32> = HashMap::new();
    let mut end_map: HashMap<usize, u32> = HashMap::new();
    let mut front_trimers: Vec<usize> = Vec::new();
    let mut end_trimers: Vec<usize> = Vec::new();
    let mut hexamer_map: HashMap<(u32, u32), usize> = HashMap::new();
    let mut hexamers: Vec<(u32, u32, u64)> = Vec::new();
    let mut total: u64 = 0;

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
        match hexamer_map.entry((front_id, end_id)) {
            std::collections::hash_map::Entry::Occupied(entry) => {
                hexamers[*entry.get()].2 += weight;
            }
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(hexamers.len());
                hexamers.push((front_id, end_id, weight));
            }
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

    Database {
        front_trimers,
        end_trimers,
        hexamers,
        front_keys,
        end_keys,
        total,
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

fn search_one(
    query: &Node,
    database: &Database,
    residue: &[i32; 400],
    table: &KmerSimilarityTable,
    seed: TerminalSeed,
) -> Row {
    // Score this query's own termini against every distinct database terminus
    // once, then reuse the two tables for every hexamer.
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

/// The same row computed the slow, obvious way: every query against every
/// database peptide, retrieval by the production index, similarity by the
/// production scorer. Used only by `--verify-search`, which exists to prove that
/// collapsing the database onto hexamers changes no number.
fn brute_force_row(
    query: &Node,
    database_nodes: &[Node],
    buckets: &[Vec<u32>],
    table: &KmerSimilarityTable,
    scorer: &Scorer,
    seed: TerminalSeed,
) -> Row {
    let mut all = [0u64; BINS];
    let mut hits = [0u64; BINS];
    let mut total: u64 = 0;
    let mut n_hits: u64 = 0;
    for node in database_nodes {
        let similarity = scorer.scores(query, node).terminal_kmer as usize;
        all[similarity] += node.frequency;
        total += node.frequency;
    }
    // `None` because the query is not a member of the database index; a query
    // that also occurs in the database must count as its own hit.
    for id in retrieve_candidates(None, query, buckets, table, seed) {
        let node = &database_nodes[id as usize];
        let similarity = scorer.scores(query, node).terminal_kmer as usize;
        hits[similarity] += node.frequency;
        n_hits += node.frequency;
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
) -> Result<(usize, usize, usize), DynError> {
    let (database_nodes, database_stats) = load_nodes(database_path, strict)?;
    let (query_nodes, query_stats) = load_nodes(query_path, strict)?;
    if database_nodes.is_empty() || query_nodes.is_empty() {
        return Err("both the database and the query set must contain usable peptides".into());
    }
    let database = build_database(&database_nodes, seed);
    let residue = normalized_residue_scores();

    let rows: Vec<Row> = query_nodes
        .par_iter()
        .map(|query| search_one(query, &database, &residue, table, seed))
        .collect();

    if verify {
        let buckets = build_exact_index(&database_nodes, seed);
        let disagreements: usize = query_nodes
            .par_iter()
            .zip(&rows)
            .filter(|(query, fast)| {
                let slow =
                    brute_force_row(query, &database_nodes, &buckets, table, scorer, seed);
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
    fn identical_termini_score_one() {
        let residue = normalized_residue_scores();
        let trimer = 1usize * 400 + 2 * 20 + 3;
        assert_eq!(terminus_score(&residue, trimer, trimer), 1000);
    }
}

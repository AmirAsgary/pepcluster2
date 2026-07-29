use crate::model::Node;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub type DynError = Box<dyn std::error::Error + Send + Sync>;

#[derive(Default, Clone, Copy)]
pub struct FastaStats {
    pub records: u64,
    pub accepted: u64,
    pub skipped: u64,
}

#[inline]
pub fn is_canonical(aa: u8) -> bool {
    matches!(
        aa,
        b'A' | b'R'
            | b'N'
            | b'D'
            | b'C'
            | b'Q'
            | b'E'
            | b'G'
            | b'H'
            | b'I'
            | b'L'
            | b'K'
            | b'M'
            | b'F'
            | b'P'
            | b'S'
            | b'T'
            | b'W'
            | b'Y'
            | b'V'
    )
}

#[inline]
pub fn aa_code(aa: u8) -> u8 {
    match aa {
        b'A' => 0,
        b'R' => 1,
        b'N' => 2,
        b'D' => 3,
        b'C' => 4,
        b'Q' => 5,
        b'E' => 6,
        b'G' => 7,
        b'H' => 8,
        b'I' => 9,
        b'L' => 10,
        b'K' => 11,
        b'M' => 12,
        b'F' => 13,
        b'P' => 14,
        b'S' => 15,
        b'T' => 16,
        b'W' => 17,
        b'Y' => 18,
        b'V' => 19,
        _ => 255,
    }
}

pub fn extract_anchor(sequence: &[u8]) -> Option<[u8; 6]> {
    if sequence.len() < 8 || sequence.iter().any(|aa| !is_canonical(*aa)) {
        return None;
    }
    let n = sequence.len();
    Some([
        sequence[0],
        sequence[1],
        sequence[2],
        sequence[n - 3],
        sequence[n - 2],
        sequence[n - 1],
    ])
}

/// Plausible terminal anchor pairs, in this fixed order:
/// N1-C1, N1-C2, N1-C3, N2-C2, N2-C3, N3-C3.
/// A bit is retained only when the two residues are at least six sequence
/// positions apart. All six hypotheses are valid for lengths >= 9.
pub fn combination_mask(sequence_length: usize) -> u8 {
    const COMBINATIONS: [(usize, usize); 6] = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];
    let mut mask = 0u8;
    for (bit, (front, back_local)) in COMBINATIONS.iter().copied().enumerate() {
        let back = sequence_length - 3 + back_local;
        if back >= front + 6 {
            mask |= 1 << bit;
        }
    }
    mask
}

pub fn scan_fasta<F>(path: &Path, strict: bool, mut callback: F) -> Result<FastaStats, DynError>
where
    F: FnMut(&str, &[u8], [u8; 6], u8) -> Result<(), DynError>,
{
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(1024 * 1024, file);
    let mut line = String::new();
    let mut header = String::new();
    let mut sequence = Vec::<u8>::new();
    let mut stats = FastaStats::default();

    let finish_record = |header: &str,
                         sequence: &[u8],
                         stats: &mut FastaStats,
                         callback: &mut F|
     -> Result<(), DynError> {
        if header.is_empty() {
            return Ok(());
        }
        stats.records += 1;
        if let Some(anchor) = extract_anchor(sequence) {
            callback(header, sequence, anchor, combination_mask(sequence.len()))?;
            stats.accepted += 1;
        } else if strict {
            return Err(format!(
                "invalid FASTA record '{}': sequence must contain at least eight canonical amino acids",
                header
            )
            .into());
        } else {
            stats.skipped += 1;
        }
        Ok(())
    };

    loop {
        line.clear();
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            finish_record(&header, &sequence, &mut stats, &mut callback)?;
            break;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            finish_record(&header, &sequence, &mut stats, &mut callback)?;
            header.clear();
            header.push_str(rest.trim());
            sequence.clear();
        } else {
            if header.is_empty() {
                return Err("FASTA sequence encountered before the first header".into());
            }
            sequence.extend(trimmed.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    Ok(stats)
}

pub fn load_nodes(path: &Path, strict: bool) -> Result<(Vec<Node>, FastaStats), DynError> {
    let mut counts: HashMap<Vec<u8>, u64> = HashMap::new();
    let stats = scan_fasta(path, strict, |_header, sequence, _anchor, _mask| {
        *counts.entry(sequence.to_vec()).or_insert(0) += 1;
        Ok(())
    })?;
    let mut nodes: Vec<Node> = counts
        .into_iter()
        .map(|(sequence, frequency)| {
            let anchor = extract_anchor(&sequence).expect("validated peptide sequence");
            let combination_mask = combination_mask(sequence.len());
            Node {
                sequence_codes: sequence.iter().copied().map(aa_code).collect(),
                sequence,
                codes: anchor.map(aa_code),
                anchor,
                combination_mask,
                frequency,
            }
        })
        .collect();
    nodes.sort_unstable_by(|a, b| a.sequence.cmp(&b.sequence));
    Ok((nodes, stats))
}

pub fn clean_field(value: &str) -> String {
    value
        .chars()
        .map(|c| {
            if c == '\t' || c == '\n' || c == '\r' {
                ' '
            } else {
                c
            }
        })
        .collect()
}

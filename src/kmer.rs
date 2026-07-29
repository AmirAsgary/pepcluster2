use crate::fasta::DynError;
use crate::scoring::normalized_residue_scores;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

pub const N_DIMERS: usize = 400;
const N_SCORES: usize = N_DIMERS * N_DIMERS;
const MAGIC: &[u8; 8] = b"PC2K2S01";

/// Threshold-independent normalized BLOSUM similarity for all ordered
/// 2-mer pairs, plus the thresholded neighbour lists used during indexing.
/// The binary table is cached on disk and loaded once per process.
pub struct KmerSimilarityTable {
    neighbours: Vec<Vec<u16>>,
}

fn rounded_half(sum: i32) -> i16 {
    if sum >= 0 {
        ((sum + 1) / 2) as i16
    } else {
        ((sum - 1) / 2) as i16
    }
}

fn generate_scores() -> Vec<i16> {
    let residue = normalized_residue_scores();
    let mut scores = vec![0i16; N_SCORES];
    for first in 0..N_DIMERS {
        let a = first / 20;
        let b = first % 20;
        for second in 0..N_DIMERS {
            let c = second / 20;
            let d = second % 20;
            scores[first * N_DIMERS + second] =
                rounded_half(residue[a * 20 + c] + residue[b * 20 + d]);
        }
    }
    scores
}

fn write_table(path: &Path, scores: &[i16]) -> Result<(), DynError> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let temporary = path.with_extension(format!("tmp-{}", std::process::id()));
    let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(&temporary)?);
    writer.write_all(MAGIC)?;
    writer.write_all(&(N_DIMERS as u32).to_le_bytes())?;
    writer.write_all(&0u32.to_le_bytes())?;
    for score in scores {
        writer.write_all(&score.to_le_bytes())?;
    }
    writer.flush()?;
    fs::rename(temporary, path)?;
    Ok(())
}

fn neighbours_from_bytes(bytes: &[u8], threshold_q: i16) -> Result<Vec<Vec<u16>>, DynError> {
    if bytes.len() != 16 + N_SCORES * 2
        || &bytes[..8] != MAGIC
        || u32::from_le_bytes(bytes[8..12].try_into()?) != N_DIMERS as u32
    {
        return Err("invalid PepCluster2 k-mer similarity table".into());
    }
    let mut neighbours = vec![Vec::new(); N_DIMERS];
    for first in 0..N_DIMERS {
        for second in 0..N_DIMERS {
            let offset = 16 + (first * N_DIMERS + second) * 2;
            let score = i16::from_le_bytes([bytes[offset], bytes[offset + 1]]);
            if score >= threshold_q {
                neighbours[first].push(second as u16);
            }
        }
    }
    Ok(neighbours)
}

#[cfg(unix)]
struct ReadOnlyMap {
    address: *mut std::ffi::c_void,
    length: usize,
}

#[cfg(unix)]
impl ReadOnlyMap {
    fn open(path: &Path) -> Result<Self, DynError> {
        use std::os::fd::AsRawFd;
        unsafe extern "C" {
            fn mmap(
                address: *mut std::ffi::c_void,
                length: usize,
                protection: i32,
                flags: i32,
                descriptor: i32,
                offset: isize,
            ) -> *mut std::ffi::c_void;
        }
        let file = File::open(path)?;
        let length = file.metadata()?.len() as usize;
        let address = unsafe {
            mmap(
                std::ptr::null_mut(),
                length,
                1, // PROT_READ
                2, // MAP_PRIVATE
                file.as_raw_fd(),
                0,
            )
        };
        if address as isize == -1 {
            return Err(std::io::Error::last_os_error().into());
        }
        Ok(Self { address, length })
    }

    fn bytes(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.address.cast::<u8>(), self.length) }
    }
}

#[cfg(unix)]
impl Drop for ReadOnlyMap {
    fn drop(&mut self) {
        unsafe extern "C" {
            fn munmap(address: *mut std::ffi::c_void, length: usize) -> i32;
        }
        unsafe {
            munmap(self.address, self.length);
        }
    }
}

fn read_neighbours(path: &Path, threshold_q: i16) -> Result<Vec<Vec<u16>>, DynError> {
    #[cfg(unix)]
    {
        let mapped = ReadOnlyMap::open(path)?;
        return neighbours_from_bytes(mapped.bytes(), threshold_q)
            .map_err(|_| format!("invalid PepCluster2 k-mer table: {}", path.display()).into());
    }
    #[cfg(not(unix))]
    {
        let bytes = fs::read(path)?;
        neighbours_from_bytes(&bytes, threshold_q)
            .map_err(|_| format!("invalid PepCluster2 k-mer table: {}", path.display()).into())
    }
}

impl KmerSimilarityTable {
    pub fn open_or_create(path: &Path, seed_threshold: f64) -> Result<Self, DynError> {
        if !(0.0..=1.0).contains(&seed_threshold) {
            return Err("--kmer-seed-threshold must be between 0 and 1".into());
        }
        if !path.exists() {
            let scores = generate_scores();
            write_table(path, &scores)?;
        }
        let threshold_q = (seed_threshold * 1000.0).round() as i16;
        let neighbours = read_neighbours(path, threshold_q)?;
        Ok(Self { neighbours })
    }

    #[inline]
    pub fn neighbours(&self, dimer: usize) -> &[u16] {
        &self.neighbours[dimer]
    }

    #[inline]
    pub fn are_neighbours(&self, first: usize, second: usize) -> bool {
        self.neighbours[first]
            .binary_search(&(second as u16))
            .is_ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expected_table_size_and_neighbour_symmetry() {
        let scores = generate_scores();
        assert_eq!(scores.len(), 160_000);
        for a in 0..400 {
            for b in 0..400 {
                assert_eq!(scores[a * 400 + b], scores[b * 400 + a]);
            }
        }
    }

    #[test]
    fn hk_and_hs_are_half_similar() {
        // AA order: ARNDCQEGHILKMFPSTWYV
        let h = 8;
        let k = 11;
        let s = 15;
        let hk = h * 20 + k;
        let hs = h * 20 + s;
        assert_eq!(generate_scores()[hk * 400 + hs], 500);
    }
}

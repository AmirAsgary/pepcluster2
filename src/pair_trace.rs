use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::{Mutex, OnceLock};

static TRACE: OnceLock<Vec<Mutex<HashSet<u64>>>> = OnceLock::new();

#[inline]
fn packed(a: u32, b: u32) -> Option<u64> {
    if a == b {
        None
    } else {
        let (u, v) = if a < b { (a, b) } else { (b, a) };
        Some(((u as u64) << 32) | v as u64)
    }
}

pub fn enable(worker_threads: usize) {
    let slots = worker_threads.max(1) + 1;
    let _ = TRACE.set((0..slots).map(|_| Mutex::new(HashSet::new())).collect());
}

#[inline]
fn slot() -> Option<&'static Mutex<HashSet<u64>>> {
    let trace = TRACE.get()?;
    let index = rayon::current_thread_index().unwrap_or(trace.len() - 1);
    trace.get(index)
}

#[inline]
pub fn record(a: u32, b: u32) {
    let Some(value) = packed(a, b) else { return };
    let Some(slot) = slot() else { return };
    slot.lock()
        .expect("pair trace mutex poisoned")
        .insert(value);
}

pub fn record_many(pairs: &[(u32, u32)]) {
    let Some(slot) = slot() else { return };
    let mut values = slot.lock().expect("pair trace mutex poisoned");
    values.extend(pairs.iter().filter_map(|&(a, b)| packed(a, b)));
}

pub fn write(path: &Path) -> std::io::Result<u64> {
    let mut values = Vec::<u64>::new();
    if let Some(trace) = TRACE.get() {
        for slot in trace {
            values.extend(
                slot.lock()
                    .expect("pair trace mutex poisoned")
                    .iter()
                    .copied(),
            );
        }
    }
    values.sort_unstable();
    values.dedup();
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, File::create(path)?);
    writer.write_all(b"PC2PAIR1")?;
    writer.write_all(&(values.len() as u64).to_le_bytes())?;
    for value in &values {
        writer.write_all(&value.to_le_bytes())?;
    }
    writer.flush()?;
    Ok(values.len() as u64)
}

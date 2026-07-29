#[derive(Clone, Debug)]
pub struct Node {
    /// Canonical full peptide sequence. Full-sequence storage is required by
    /// constrained alignment modes and makes node identity independent of the
    /// terminal blocking representation.
    pub sequence: Vec<u8>,
    pub sequence_codes: Vec<u8>,
    pub anchor: [u8; 6],
    pub codes: [u8; 6],
    /// Bit i marks whether anchor-combination hypothesis i is geometrically
    /// possible for the peptide length represented by this node.
    pub combination_mask: u8,
    pub frequency: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Edge {
    pub u: u32,
    pub v: u32,
    /// Mode-specific ranking weight quantized to 0..=1000. In combined modes
    /// this is the combined score; in separate mode it is the weakest
    /// threshold margin.
    pub weight: u16,
}

#[derive(Clone, Copy, Debug)]
pub struct Neighbor {
    pub node: u32,
    pub weight: u16,
}

pub struct Graph {
    pub offsets: Vec<usize>,
    pub neighbors: Vec<Neighbor>,
    pub edge_count: u64,
}

impl Graph {
    #[inline]
    pub fn neighbors(&self, node: u32) -> &[Neighbor] {
        let i = node as usize;
        &self.neighbors[self.offsets[i]..self.offsets[i + 1]]
    }

    pub fn edge_weight(&self, a: u32, b: u32) -> Option<u16> {
        if a == b {
            return Some(1000);
        }
        self.neighbors(a)
            .binary_search_by_key(&b, |n| n.node)
            .ok()
            .map(|i| self.neighbors(a)[i].weight)
    }
}

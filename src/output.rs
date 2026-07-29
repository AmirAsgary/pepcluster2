use crate::edge_store::read_edge;
use crate::fasta::{clean_field, scan_fasta, DynError};
use crate::graph::Clustering;
use crate::model::Node;
use crate::scoring::Scorer;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

pub fn anchor_string(anchor: &[u8; 6]) -> String {
    String::from_utf8_lossy(anchor).into_owned()
}

fn cluster_name(cluster: u32, width: usize) -> String {
    format!("PC2_{:0width$}", cluster + 1, width = width)
}

pub fn write_cluster_outputs(
    input: &Path,
    output_dir: &Path,
    strict: bool,
    nodes: &[Node],
    scorer: &Scorer,
    clustering: &Clustering,
    min_cluster_size: usize,
    write_cluster_fastas: bool,
) -> Result<(), DynError> {
    let fasta_dir = output_dir.join("cluster_fastas");
    if write_cluster_fastas && fasta_dir.exists() {
        return Err(format!(
            "{} already exists; use a new output directory or remove it before rerunning",
            fasta_dir.display()
        )
        .into());
    }

    struct Record {
        header: String,
        sequence: Vec<u8>,
        anchor: [u8; 6],
        combination_mask: u8,
        node_id: u32,
    }

    let cluster_count = clustering.representatives.len();
    let mut records = Vec::<Record>::new();
    let mut cluster_sizes = vec![0usize; cluster_count];
    let mut unique_anchor_counts = vec![0usize; cluster_count];
    for &cluster in &clustering.cluster_of {
        unique_anchor_counts[cluster as usize] += 1;
    }
    let mut representative_record: Vec<Option<usize>> = vec![None; cluster_count];
    scan_fasta(
        input,
        strict,
        |header, sequence, anchor, combination_mask| {
            let node_id = nodes
                .binary_search_by(|node| node.sequence.as_slice().cmp(sequence))
                .map_err(|_| "internal error: peptide missing from node table")?
                as u32;
            let cluster = clustering.cluster_of[node_id as usize] as usize;
            cluster_sizes[cluster] += 1;
            let record_id = records.len();
            records.push(Record {
                header: clean_field(header),
                sequence: sequence.to_vec(),
                anchor,
                combination_mask,
                node_id,
            });
            if node_id == clustering.representatives[cluster] {
                let replace = representative_record[cluster].is_none_or(|old| {
                    records[record_id]
                        .sequence
                        .cmp(&records[old].sequence)
                        .then(records[record_id].header.cmp(&records[old].header))
                        .is_lt()
                });
                if replace {
                    representative_record[cluster] = Some(record_id);
                }
            }
            Ok(())
        },
    )?;

    let mut cluster_order: Vec<usize> = (0..cluster_count).collect();
    cluster_order.sort_unstable_by(|&a, &b| {
        cluster_sizes[b]
            .cmp(&cluster_sizes[a])
            .then(clustering.representatives[a].cmp(&clustering.representatives[b]))
    });
    let mut output_rank = vec![0u32; cluster_count];
    for (rank, &cluster) in cluster_order.iter().enumerate() {
        output_rank[cluster] = rank as u32;
        if representative_record[cluster].is_none() {
            return Err(format!("representative peptide missing for cluster {cluster}").into());
        }
    }
    let width = cluster_count.max(1).to_string().len().max(6);

    let mut summary = BufWriter::with_capacity(
        1024 * 1024,
        File::create(output_dir.join("cluster_summary.tsv"))?,
    );
    writeln!(
        summary,
        "cluster_id\trepresentative_anchor\trepresentative_peptide\tsize\tunique_sequences"
    )?;
    let mut representatives = BufWriter::with_capacity(
        1024 * 1024,
        File::create(output_dir.join("cluster_representatives.tsv"))?,
    );
    writeln!(
        representatives,
        "cluster_id\trepresentative_anchor\tfrequency\trepresentative_header\trepresentative_sequence"
    )?;
    for (rank, &cluster) in cluster_order.iter().enumerate() {
        let rep = clustering.representatives[cluster];
        let record = &records[representative_record[cluster].unwrap()];
        let name = cluster_name(rank as u32, width);
        writeln!(
            summary,
            "{}\t{}\t{}\t{}\t{}",
            name,
            anchor_string(&nodes[rep as usize].anchor),
            String::from_utf8_lossy(&record.sequence),
            cluster_sizes[cluster],
            unique_anchor_counts[cluster]
        )?;
        writeln!(
            representatives,
            "{}\t{}\t{}\t{}\t{}",
            name,
            anchor_string(&nodes[rep as usize].anchor),
            nodes[rep as usize].frequency,
            record.header,
            String::from_utf8_lossy(&record.sequence)
        )?;
    }
    summary.flush()?;
    representatives.flush()?;

    let mut anchors = BufWriter::with_capacity(
        4 * 1024 * 1024,
        File::create(output_dir.join("anchor_clusters.tsv"))?,
    );
    writeln!(
        anchors,
        "cluster_id\trepresentative_sequence\tsequence\tfrequency\tterminal_kmer_similarity\talignment_similarity\tanchor_combination_similarity\tcombined_score\trepresentative_ranking_weight"
    )?;
    let mut node_order: Vec<usize> = (0..nodes.len()).collect();
    node_order.sort_unstable_by_key(|&node| {
        let cluster = clustering.cluster_of[node] as usize;
        (
            output_rank[cluster],
            node as u32 != clustering.representatives[cluster],
            nodes[node].anchor,
            nodes[node].combination_mask,
        )
    });
    for node_id in node_order {
        let node = &nodes[node_id];
        let cluster = clustering.cluster_of[node_id];
        let rep = clustering.representatives[cluster as usize];
        let scores = scorer.scores(node, &nodes[rep as usize]);
        writeln!(
            anchors,
            "{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
            cluster_name(output_rank[cluster as usize], width),
            String::from_utf8_lossy(&nodes[rep as usize].sequence),
            String::from_utf8_lossy(&node.sequence),
            node.frequency,
            scores.terminal_kmer as f64 / 1000.0,
            scores.alignment as f64 / 1000.0,
            scores.anchor_combination as f64 / 1000.0,
            scores.combined as f64 / 1000.0,
            scores.ranking_weight as f64 / 1000.0
        )?;
    }
    anchors.flush()?;

    let mut assignments = BufWriter::with_capacity(
        8 * 1024 * 1024,
        File::create(output_dir.join("clusters.tsv"))?,
    );
    writeln!(
        assignments,
        "cluster_id\trepresentative_anchor\trepresentative_peptide\tpeptide_header\tsequence\tanchor\tgeometry_mask\talignment_similarity\tanchor_combination_similarity\tcombined_score\trepresentative_ranking_weight"
    )?;
    let mut record_order: Vec<usize> = (0..records.len()).collect();
    record_order.sort_unstable_by(|&a, &b| {
        let ca = clustering.cluster_of[records[a].node_id as usize] as usize;
        let cb = clustering.cluster_of[records[b].node_id as usize] as usize;
        output_rank[ca]
            .cmp(&output_rank[cb])
            .then(records[a].sequence.cmp(&records[b].sequence))
            .then(records[a].header.cmp(&records[b].header))
    });
    for &record_id in &record_order {
        let record = &records[record_id];
        let cluster = clustering.cluster_of[record.node_id as usize] as usize;
        let rep = clustering.representatives[cluster];
        let scores = scorer.scores(&nodes[record.node_id as usize], &nodes[rep as usize]);
        writeln!(
            assignments,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
            cluster_name(output_rank[cluster], width),
            anchor_string(&nodes[rep as usize].anchor),
            String::from_utf8_lossy(&records[representative_record[cluster].unwrap()].sequence),
            record.header,
            String::from_utf8_lossy(&record.sequence),
            anchor_string(&record.anchor),
            record.combination_mask,
            scores.alignment as f64 / 1000.0,
            scores.anchor_combination as f64 / 1000.0,
            scores.combined as f64 / 1000.0,
            scores.ranking_weight as f64 / 1000.0
        )?;
    }
    assignments.flush()?;

    let mut fasta_count = 0usize;
    if write_cluster_fastas {
        fs::create_dir_all(&fasta_dir)?;
        let mut current_rank = u32::MAX;
        let mut fasta_writer: Option<BufWriter<File>> = None;
        for &record_id in &record_order {
            let record = &records[record_id];
            let cluster = clustering.cluster_of[record.node_id as usize] as usize;
            let rank = output_rank[cluster];
            if rank != current_rank {
                if let Some(mut writer) = fasta_writer.take() {
                    writer.flush()?;
                }
                current_rank = rank;
                if cluster_sizes[cluster] >= min_cluster_size {
                    let name = cluster_name(rank, width);
                    fasta_writer = Some(BufWriter::new(File::create(
                        fasta_dir.join(format!("{name}.fasta")),
                    )?));
                    fasta_count += 1;
                }
            }
            if let Some(writer) = fasta_writer.as_mut() {
                writeln!(
                    writer,
                    ">{}\n{}",
                    record.header,
                    String::from_utf8_lossy(&record.sequence)
                )?;
            }
        }
        if let Some(mut writer) = fasta_writer {
            writer.flush()?;
        }
    }

    let singleton_count = cluster_sizes.iter().filter(|size| **size == 1).count();
    let mut sorted_sizes = cluster_sizes.clone();
    sorted_sizes.sort_unstable();
    let median = if sorted_sizes.is_empty() {
        0.0
    } else if sorted_sizes.len() % 2 == 1 {
        sorted_sizes[sorted_sizes.len() / 2] as f64
    } else {
        let i = sorted_sizes.len() / 2;
        (sorted_sizes[i - 1] + sorted_sizes[i]) as f64 / 2.0
    };
    let mean = records.len() as f64 / cluster_count.max(1) as f64;
    let maximum = cluster_sizes.iter().copied().max().unwrap_or(0);
    let report = format!(
        concat!(
            "PEPCLUSTER2 RUN SUMMARY\n",
            "=======================\n",
            "Accepted peptides:       {}\n",
            "Unique peptide sequences: {}\n",
            "Clusters:                {}\n",
            "Singleton clusters:      {}\n",
            "Non-singleton clusters:  {}\n",
            "Mean cluster size:       {:.3}\n",
            "Median cluster size:     {:.3}\n",
            "Largest cluster:         {}\n",
            "Per-cluster FASTA files: {} (minimum size {}; requested {})\n"
        ),
        records.len(),
        nodes.len(),
        cluster_count,
        singleton_count,
        cluster_count - singleton_count,
        mean,
        median,
        maximum,
        fasta_count,
        min_cluster_size,
        write_cluster_fastas
    );
    fs::write(output_dir.join("run_summary.txt"), report)?;
    Ok(())
}

pub fn write_provisional_clusters(
    output_file: &Path,
    nodes: &[Node],
    clustering: &Clustering,
) -> Result<(), DynError> {
    let mut sizes = vec![0usize; clustering.representatives.len()];
    for &cluster in &clustering.cluster_of {
        sizes[cluster as usize] += 1;
    }
    let mut writer = BufWriter::new(File::create(output_file)?);
    writeln!(writer, "node_id\tanchor\tgeometry_mask\tprovisional_cluster\tprovisional_representative_node\tprovisional_cluster_size\tassigned_to_nonsingleton")?;
    for (node, data) in nodes.iter().enumerate() {
        let cluster = clustering.cluster_of[node] as usize;
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            node,
            anchor_string(&data.anchor),
            data.combination_mask,
            cluster,
            clustering.representatives[cluster],
            sizes[cluster],
            sizes[cluster] > 1
        )?;
    }
    writer.flush()?;
    Ok(())
}

/// Minimal deterministic scientific output for high-replicate validation.
/// It deliberately avoids rescanning the FASTA and creating per-cluster files.
pub fn write_compact_outputs(
    output_dir: &Path,
    nodes: &[Node],
    scorer: &Scorer,
    clustering: &Clustering,
) -> Result<(), DynError> {
    let cluster_count = clustering.representatives.len();
    let width = cluster_count.max(1).to_string().len().max(6);
    let mut writer = BufWriter::with_capacity(
        4 * 1024 * 1024,
        File::create(output_dir.join("node_clusters.tsv"))?,
    );
    writeln!(
        writer,
        "cluster_id\trepresentative_sequence\tsequence\tfrequency\tterminal_kmer_similarity\talignment_similarity\tanchor_combination_similarity\tcombined_score\trepresentative_ranking_weight"
    )?;
    for (node_id, node) in nodes.iter().enumerate() {
        let cluster = clustering.cluster_of[node_id] as usize;
        let representative = clustering.representatives[cluster];
        let scores = scorer.scores(node, &nodes[representative as usize]);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
            cluster_name(cluster as u32, width),
            String::from_utf8_lossy(&nodes[representative as usize].sequence),
            String::from_utf8_lossy(&node.sequence),
            node.frequency,
            scores.terminal_kmer as f64 / 1000.0,
            scores.alignment as f64 / 1000.0,
            scores.anchor_combination as f64 / 1000.0,
            scores.combined as f64 / 1000.0,
            scores.ranking_weight as f64 / 1000.0
        )?;
    }
    writer.flush()?;
    Ok(())
}

pub fn write_edges(
    edge_file: &Path,
    output_file: &Path,
    nodes: &[Node],
    scorer: &Scorer,
) -> Result<(), DynError> {
    let mut reader = BufReader::with_capacity(4 * 1024 * 1024, File::open(edge_file)?);
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, File::create(output_file)?);
    writeln!(writer, "sequence_a\tsequence_b\tterminal_kmer_similarity\talignment_similarity\tanchor_combination_similarity\tcombined_score\trepresentative_ranking_weight")?;
    while let Some(edge) = read_edge(&mut reader)? {
        let scores = scorer.scores(&nodes[edge.u as usize], &nodes[edge.v as usize]);
        writeln!(
            writer,
            "{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
            String::from_utf8_lossy(&nodes[edge.u as usize].sequence),
            String::from_utf8_lossy(&nodes[edge.v as usize].sequence),
            scores.terminal_kmer as f64 / 1000.0,
            scores.alignment as f64 / 1000.0,
            scores.anchor_combination as f64 / 1000.0,
            scores.combined as f64 / 1000.0,
            scores.ranking_weight as f64 / 1000.0
        )?;
    }
    writer.flush()?;
    Ok(())
}

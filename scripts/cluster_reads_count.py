#!/usr/bin/env python3
import os
import sys
from pathlib import Path
import re

def main(cluster_dir, fastq_folder, rawseqids_dir):
    cluster_path = Path(cluster_dir)
    fastq_path = Path(fastq_folder)
    rawseqids_path = Path(rawseqids_dir)

    if not cluster_path.is_dir():
        sys.exit(f"Error: Cluster directory '{cluster_dir}' does not exist.")
    if not fastq_path.is_dir():
        sys.exit(f"Error: FASTQ folder '{fastq_folder}' does not exist.")
    if not rawseqids_path.is_dir():
        sys.exit(f"Error: Rawseqids folder '{rawseqids_dir}' does not exist.")

    sample_seqids = {}
    for sample_file in fastq_path.glob('*.fastq'):
        sample_name = sample_file.stem
        rawseq_file = rawseqids_path / f"rawseqids_{sample_name}"

        if rawseq_file.is_file():
            with rawseq_file.open('r') as f:
                sample_seqids[sample_name] = set(line.strip() for line in f)
        else:
            print(f"Warning: '{rawseq_file}' not found, skipping.")
    
    # Processing clusters
    cluster_file_pattern = re.compile(r"^cluster\d+$")

    for cluster_file in cluster_path.iterdir():
        if cluster_file.is_file() and cluster_file_pattern.match(cluster_file.name):
            cluster_name = cluster_file.name
            output_file = Path(f"{cluster_name}_seq_count")

            with cluster_file.open('r') as f:
                cluster_seqs = {line[1:].strip() for line in f if line.startswith('>')}

            with output_file.open('w') as out_f:
                for sample_name, seqids in sample_seqids.items():
                    count = len(cluster_seqs & seqids)
                    out_f.write(f"{sample_name}\t{count}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: cluster_reads_count.py <cluster_dir> <fastq_folder> <rawseqids_dir>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])

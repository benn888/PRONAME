#!/usr/bin/env python3

import os
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
import pandas as pd

# Retrieving the 'HQ_fastq' environment variable
input_fasta = os.getenv('HQ_fastq')
if input_fasta is None:
    raise ValueError("The environment variable 'HQ_fastq' is not defined.")

# Other path
cluster_tsv = "output_cluster.tsv"

# Output directory
output_dir = Path("mmseqs2_clusters")
output_dir.mkdir(exist_ok=True)

# Reading the clustering file (centroids, sequences)
df = pd.read_csv(cluster_tsv, sep='\t', header=None, names=["Centroid", "Sequence"])

# Creating a cluster dictionary
clusters = defaultdict(list)
for _, row in df.iterrows():
    clusters[row["Centroid"]].append(row["Sequence"])

# Load all sequences from fasta file into memory
seq_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))

# Writing fasta files
for i, (centroid_id, seq_ids) in enumerate(clusters.items(), 1):
    cluster_path = output_dir / f"cluster{i}"
    centroid_path = output_dir / f"centroid_cluster{i}.fasta"

    # Writing all sequences from the cluster
    cluster_seqs = [seq_dict[sid] for sid in seq_ids if sid in seq_dict]
    SeqIO.write(cluster_seqs, cluster_path, "fasta")

    # Writing only centroid sequence
    if centroid_id in seq_dict:
        SeqIO.write(seq_dict[centroid_id], centroid_path, "fasta")

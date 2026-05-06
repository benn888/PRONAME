#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO

def delete_small_clusters(cluster_folder, min_reads):
    """
    Delete cluster files containing less than 'min_reads' sequences.
    By default, min_reads = 2 => removes singletons.
    """
    for file in os.listdir(cluster_folder):
        # Keep only files named as "clusterX" (without extension)
        if not file.startswith("cluster") or '.' in file:
            continue

        file_path = os.path.join(cluster_folder, file)

        with open(file_path, "r") as f:
            records = list(SeqIO.parse(f, "fasta"))

        if len(records) < min_reads:
            # Delete the file clusterX
            os.remove(file_path)
            print(f"The file {file} has been deleted because it contains only {len(records)} sequence(s).")

            # Build the corresponding centroid name
            centroid_file = f"centroid_{file}.fasta"
            centroid_path = os.path.join(cluster_folder, centroid_file)

            # Remove the centroid if it exists
            if os.path.exists(centroid_path):
                os.remove(centroid_path)
        else:
            print(f"The file {file} contains {len(records)} sequences and was not deleted.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove cluster fasta files containing fewer than a given number of sequences."
    )
    parser.add_argument("cluster_folder", help="Path to the folder containing fasta files")
    parser.add_argument(
        "--min-reads",
        type=int,
        default=2,
        help="Minimum number of sequences required to keep a cluster (default: 2 => remove singletons)."
    )

    args = parser.parse_args()

    if args.min_reads < 1:
        raise ValueError("The --min-reads value must be >= 1.")

    delete_small_clusters(args.cluster_folder, args.min_reads)

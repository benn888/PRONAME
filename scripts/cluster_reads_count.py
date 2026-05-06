#!/usr/bin/env python3
import sys
from pathlib import Path
import re


def normalize_id(s: str):
    """
    Normalise un identifiant de séquence pour le comptage :
    - retire '>' ou '@' en début de ligne
    - garde uniquement le premier champ avant les espaces/tabulations
    - coupe aussi avant un éventuel ';' (ex. vsearch ;size=XX;)
    """
    s = s.strip()
    if not s:
        return None

    if s[0] in (">", "@"):
        s = s[1:]

    # Premier champ avant espace/tab
    s = s.split()[0]

    # Si vsearch / autres ajoutent ;size=XX; etc.
    s = s.split(";")[0]

    return s


def main(cluster_dir, fastq_folder, rawseqids_dir):
    cluster_path = Path(cluster_dir)
    fastq_path = Path(fastq_folder)
    rawseqids_path = Path(rawseqids_dir)

    if not cluster_path.is_dir():
        sys.exit(f"Error: Cluster directory '{cluster_dir}' does not exist.")
    if not fastq_path.is_dir():
        sys.exit(f"Error: FASTQ folder '{fastq_folder}' does not exist.")
    if not rawseqids_path.is_dir():
        sys.exit(f"Error: Rawseqids directory '{rawseqids_dir}' does not exist.")

    # 1) Charger les ID par échantillon à partir de Rawseqids
    sample_seqids = {}

    for fastq_file in sorted(fastq_path.glob("*.fastq")):
        sample_name = fastq_file.stem
        raw_file = rawseqids_path / f"rawseqids_{sample_name}"

        if not raw_file.is_file():
            # On ignore poliment les FASTQ sans fichier Rawseqids
            continue

        ids = set()
        with raw_file.open("r") as f:
            for line in f:
                nid = normalize_id(line)
                if nid:
                    ids.add(nid)

        sample_seqids[sample_name] = ids

    if not sample_seqids:
        sys.exit("Error: No Rawseqids files could be loaded. Check 'Rawseqids/' and fastq filenames.")

    # 2) Comptage par cluster
    # On prend tous les fichiers dont le nom commence par 'cluster'
    cluster_files = sorted(cf for cf in cluster_path.iterdir() if cf.is_file() and cf.name.startswith("cluster"))

    if not cluster_files:
        sys.exit(f"Error: No cluster files found in '{cluster_dir}' (expected names like 'cluster1', 'cluster2', ...).")

    for cluster_file in cluster_files:
        cluster_name = cluster_file.stem  # 'cluster1' même si 'cluster1.fasta'

        # Récupérer les IDs des séquences du cluster
        cluster_seqs = set()
        with cluster_file.open("r") as f:
            for line in f:
                if line.startswith(">"):
                    nid = normalize_id(line)
                    if nid:
                        cluster_seqs.add(nid)

        # Si cluster vide, on met quand même des zéros pour cohérence
        output_file = Path(f"{cluster_name}_seq_count")
        with output_file.open("w") as out_f:
            for sample_name, seqids in sample_seqids.items():
                count = len(cluster_seqs & seqids)
                out_f.write(f"{sample_name}\t{count}\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: cluster_reads_count.py <cluster_dir> <fastq_folder> <rawseqids_dir>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])

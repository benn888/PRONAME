#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO

def delete_singletons(cluster_folder):
    # Browsing through all the files of the cluster folder
    for file in os.listdir(cluster_folder):
        file_path = os.path.join(cluster_folder, file)
        # Reading the fasta file
        with open(file_path, "r") as f:
            records = list(SeqIO.parse(f, "fasta"))
            
        # Checking the number of sequences
        if len(records) <= 1:
            # Deleting the file if it contains only one sequence
            os.remove(file_path)
            print(f"The file {file} has been deleted because it contains only one sequence.")
        else:
            print(f"The file {file} contains multiple sequences and was not deleted.")

if __name__ == "__main__":
    # Definition of the argument parser
    parser = argparse.ArgumentParser(description="Remove fasta files containing only one sequence.")
    parser.add_argument("cluster_folder", help="Path to the folder containing fasta files")
    
    # Parse command-line arguments
    args = parser.parse_args()
    
    # Call the delete_singletons function with the folder path provided by the user
    delete_singletons(args.cluster_folder)

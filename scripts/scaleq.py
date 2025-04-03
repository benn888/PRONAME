#!/usr/bin/env python3

"""
This script processes a FASTQ file to extract read lengths and average quality scores,
then generates a scatter plot (read length vs. quality) along with marginal histograms.
It is designed to provide a visualization similar to NanoPlot.

Key features:
- Computes average Phred quality scores using the NanoPlot method.
- Produces a scatter plot with transparency to highlight density variations.
- Includes histograms of read lengths and quality distributions.
- Automatically adjusts axis scaling and tick spacing for clarity.

Usage example:
    scaleq.py --fastq example.fastq

Requirements:
- Python 3
- NumPy
- Matplotlib
- pysam

Help:
Use the following command to display help:
    scaleq.py --help
"""

import sys
import os

def print_help():
    """Displays help information and exits."""
    help_text = """
Usage: scaleq.py --fastq <path_to_fastq_file>

Options:
  --fastq   Path to the FASTQ file to analyze (required).
  --help    Show this help message and exit.

Description:
This script processes a FASTQ file to extract read lengths and average quality scores,
then generates a scatter plot (read length vs. quality) along with marginal histograms.
It is designed to provide a visualization similar to NanoPlot.

Example:
  scaleq.py --fastq example.fastq
"""
    print(help_text)
    sys.exit(0)

# Check for help argument before importing other packages
if "--help" in sys.argv:
    print_help()

import numpy as np
import matplotlib.pyplot as plt
import argparse
import pysam
from math import log, ceil, floor

def errs_tab(n):
    """Generates a list of error probabilities for Phred quality scores."""
    return [10 ** (q / -10) for q in range(n + 1)]

# Precompute the error probability table up to Phred 128
ERRS_TABLE = errs_tab(128)

def ave_qual(quals, tab=ERRS_TABLE):
    """Computes the average base call quality using NanoPlot's approach."""
    if len(quals) > 0:
        return -10 * log(sum(tab[q] for q in quals) / len(quals), 10)
    else:
        return 0

def parse_fastq_lengths_qualities(fastq_file):
    """Extracts read lengths and average Phred quality scores from a FASTQ file."""
    lengths = []
    qualities = []

    with pysam.FastxFile(fastq_file) as fh:
        for entry in fh:
            if entry.get_quality_array() is None:
                continue  # Skip reads with missing quality scores

            seq_length = len(entry.sequence)
            quality_scores = list(entry.get_quality_array())
            avg_quality = ave_qual(quality_scores)  # Correct quality computation
            
            lengths.append(seq_length)
            qualities.append(avg_quality)

    return np.array(lengths), np.array(qualities)

def plot_length_vs_quality(fastq_file):
    """Generates a Length vs. Quality scatter plot with marginal histograms."""
    # Extract data
    lengths, qualities = parse_fastq_lengths_qualities(fastq_file)

    if len(lengths) == 0 or len(qualities) == 0:
        print("Error: The FASTQ file contains no usable data.")
        return

    # Determine axis limits
    max_quality = ceil(max(qualities) / 5) * 5  # Round up to the nearest multiple of 5
    max_length = ceil(max(lengths) / 100) * 100  # Round up to the nearest hundred

    # Define tick spacing with rounded hundred increments
    tick_spacing = max(100, round(max_length / 10, -2))

    # Generate rectangular figure
    fig = plt.figure(figsize=(12, 6))
    grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)

    main_ax = fig.add_subplot(grid[1:, :-1])  # Scatter plot (center)
    hist_x = fig.add_subplot(grid[0, :-1], sharex=main_ax)  # Length histogram (top)
    hist_y = fig.add_subplot(grid[1:, -1], sharey=main_ax)  # Quality histogram (right)

    # Scatter plot with small, semi-transparent points
    main_ax.scatter(lengths, qualities, alpha=0.1, s=2, c="royalblue")  # Semi-transparent for density visualization
    main_ax.set_xlabel("Read Length")
    main_ax.set_ylabel("Average Quality Score")
    main_ax.grid(True, linestyle="solid", color="lightgray", alpha=0.2)
    main_ax.set_xticks(np.arange(0, max_length + 1, tick_spacing))
    main_ax.set_yticks(np.arange(0, max_quality + 1, 5))
    main_ax.set_xlim(0, max_length)
    main_ax.set_ylim(0, max_quality)

    # Format tick labels without decimals
    main_ax.tick_params(axis='both', which='both', labelsize=10)
    main_ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
    main_ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{int(y)}"))

    # Length histogram (X-axis) with adjusted bins
    hist_x.hist(lengths, bins=100, color="royalblue", alpha=0.7, density=True)
    hist_x.set_ylabel("Density")
    hist_x.grid(True, linestyle="solid", color="lightgray", alpha=0.2)
    hist_x.tick_params(axis='x', labelbottom=False)

    # Quality histogram (Y-axis) with adjusted bins
    hist_y.hist(qualities, bins=100, orientation="horizontal", color="royalblue", alpha=0.7, density=True)
    hist_y.set_xlabel("Density")
    hist_y.grid(True, linestyle="solid", color="lightgray", alpha=0.2)
    hist_y.tick_params(axis='y', labelleft=False)

    # Save the plot
    plt.savefig("length_vs_quality_plot.png", dpi=300, bbox_inches="tight")

def main():
    parser = argparse.ArgumentParser(description="Generate a Length vs. Quality scatter plot from a FASTQ file.")
    parser.add_argument("--fastq", required=True, help="Path to the FASTQ file to analyze.")
    args = parser.parse_args()

    # Check if the file exists
    if not os.path.isfile(args.fastq):
        print(f"Error: The file {args.fastq} does not exist.")
        return

    # Generate the plot
    plot_length_vs_quality(args.fastq)

if __name__ == "__main__":
    main()

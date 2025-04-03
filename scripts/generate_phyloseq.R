#!/usr/bin/env Rscript
library(phyloseq)
library(Biostrings)
library(ape)

# Loading input files
args <- commandArgs(trailingOnly = TRUE)
otu_table_file <- args[1]
tax_table_file <- args[2]
rep_seqs_file <- args[3]
phy_tree_file <- ifelse(length(args) > 3, args[4], NULL)

# Creating the OTU table
otu_table_data <- read.table(otu_table_file, header = TRUE, row.names = 1, sep = "\t")
otu_table <- otu_table(as.matrix(otu_table_data), taxa_are_rows = TRUE)

# Creating the taxonomy table
tax_table_data <- read.table(tax_table_file, header = TRUE, row.names = 1, sep = "\t")
tax_table <- tax_table(as.matrix(tax_table_data))

# Creating the reference sequences
rep_seqs <- readDNAStringSet(rep_seqs_file)

# Optionally loading the phylogenetic tree
if (!is.null(phy_tree_file)) {
  phy_tree <- read.tree(phy_tree_file)
}

# Combining into a phyloseq object
physeq <- phyloseq(otu_table, tax_table, refseq(rep_seqs))

if (!is.null(phy_tree_file)) {
  physeq <- merge_phyloseq(physeq, phy_tree(phy_tree))
}

# Saving phyloseq object to an RDS file
saveRDS(physeq, file = "phyloseq_object.rds")


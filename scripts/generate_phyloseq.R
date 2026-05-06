#!/usr/bin/env Rscript
library(phyloseq)
library(Biostrings)
library(ape)

args <- commandArgs(trailingOnly = TRUE)

otu_file  <- args[1]
tax_file  <- args[2]
rep_file  <- args[3]
tree_file <- args[4]
meta_file <- ifelse(length(args) >= 5, args[5], NA)

# OTU table
otu_tab <- read.table(otu_file,
                      header = TRUE,
                      row.names = 1,
                      sep = "\t",
                      check.names = FALSE,
                      comment.char = "")
OTU <- otu_table(as.matrix(otu_tab), taxa_are_rows = TRUE)

# Taxonomy table (headerless TSV: FeatureID<TAB>TaxonomyString)
tax_raw <- read.table(tax_file,
                      header = FALSE,
                      sep = "\t",
                      quote = "",
                      comment.char = "",
                      stringsAsFactors = FALSE)
colnames(tax_raw)[1:2] <- c("FeatureID", "Taxonomy")
rownames(tax_raw) <- tax_raw$FeatureID
tax_mat <- as.matrix(tax_raw["Taxonomy"])
TAX <- tax_table(tax_mat)

# Reference sequences
REF <- refseq(readDNAStringSet(rep_file))

# Build phyloseq (without tree first)
physeq <- phyloseq(OTU, TAX, REF)

# Optional tree
if (!is.na(tree_file) && file.exists(tree_file)) {
  tr <- read.tree(tree_file)
  physeq <- merge_phyloseq(physeq, phy_tree(tr))
}

# Optional metadata
if (!is.na(meta_file) && file.exists(meta_file)) {
  meta <- read.table(meta_file,
                     header = TRUE,
                     row.names = 1,
                     sep = "\t",
                     check.names = FALSE,
                     comment.char = "")
  physeq <- merge_phyloseq(physeq, sample_data(meta))
}

saveRDS(physeq, file = "phyloseq_object.rds")

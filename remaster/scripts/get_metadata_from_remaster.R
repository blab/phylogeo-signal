library(treeio) # for read.beast file
library(tidyverse)

source('scripts/utils_remaster.R')

args <- commandArgs(trailingOnly = T)
mutation_scenario <- as.character(args[1])
migration_scenario <- as.character(args[2])
tree_file_path <- as.character(args[3])
output_file_metadata <- as.character(args[4])

# Read transmission tree (with associated node characteristics)
phylo_tree <- read.beast(file = tree_file_path)
  
# Extract metadata
write_metadata_from_phylo_tree(output_file_metadata, phylo_tree)

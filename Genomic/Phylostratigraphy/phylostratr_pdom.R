library(phylostratr)
library(tidyverse)

setwd("/home/benjamin/Documents/LF_2020_repo/Genomic/Phylostratigraphy/Pdom")

weights=uniprot_weight_by_ref()
focal_taxid = '743375'
strata =
  # Get stratified relatives represented in UniProt
  uniprot_strata(focal_taxid, from=2) %>%
  # Select a diverse subset of 5 or fewer representatives from each stratum.
  strata_apply(f=diverse_subtree, n=5, weights=weights) %>%
  # Use a prebuilt set of prokaryotic species
  use_recommended_prokaryotes %>%
  # Add yeast and humans
  add_taxa(c('4932', '9606')) %>%
  # Download genomes, storing the filenames
  uniprot_fill_strata

# LF data aren;t on UniProt, so add them manually
strata@data$faa[['743375']] = '/home/benjamin/Documents/LF_2020_repo/Genomic/Polistes_dominula.faa'

# for testing purposes, take just 1% of genes
#strata@data$faa[['85445']] <- thin_fasta(strata@data$faa[['85445']], 100)

# Plot tree for included species
strata %>% strata_convert(target='all', to='name') %>% sort_strata %>% plot

# BLAST against each target genome (this will take a few hours)
strata <- strata_blast(strata, makedb_args=list(verbose=T), blast_args=list(nthreads=8)) %>% strata_besthits
# Merge results into a single table
results <- merge_besthits(strata)
phylostrata <- stratify(results)
table(phylostrata$mrca_name)

#save
write.csv(phylostrata, file = "/home/benjamin/Documents/LF_2020_repo/Genomic/Phylostratigraphy/phylostrata_pdom.csv", row.names=F)

phylostrata_pdom = read.csv("/home/benjamin/Documents/LF_2020_repo/Genomic/Phylostratigraphy/phylostrata_pdom.csv")

phylostrata_pdom

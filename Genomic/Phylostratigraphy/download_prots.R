#get libraries
basic_libraries <- c("tidyverse",
                     "devtools")
for (lib in basic_libraries) {
  if (require(package = lib, character.only = TRUE)) {
    print("Successful")
  } else {
    print("Installing")
    
    install.packages(lib, repos="http://cran.r-project.org")
    library(lib, character.only = TRUE )
  }
}

#install phylostratr from github
install_github("arendsee/phylostratr")
library(phylostratr)

# the below will download all the relevant fastas from UniProt
weights=uniprot_weight_by_ref()
focal_taxid = '85445' # LF is the focal taxon
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

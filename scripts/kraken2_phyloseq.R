# Function to convert preprocessed data and metadata into a phyloseq object
kraken2_phyloseq <- function(kraken2_preprocess_output, metadata, keep_bacteria_only = TRUE) {
  # Get data ready for phyloseq creation
  otu_mat <- as.matrix(kraken2_preprocess_output$kraken2_otu_table)
  tax_mat <- as.matrix(kraken2_preprocess_output$kraken2_tax_table)
  samples_df <- metadata
  OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX <- tax_table(tax_mat)
  samples <- sample_data(samples_df)
  sample_names(samples) <- colnames(otu_mat)
  
  # Import data into a phyloseq object
  bact_shotgun <- phyloseq(OTU, TAX, samples)
  
  # Retain only bacterial taxa
  bact_shotgun <- subset_taxa(bact_shotgun, kingdom == 'Bacteria')
  
  # Return phyloseq
  return(bact_shotgun)
}
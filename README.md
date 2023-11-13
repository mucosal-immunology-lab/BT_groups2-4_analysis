# Multi-omic data analysis of groups 2-4 of the Breathing Together study

Here we present the data processing and subsequent multi-omic integration and analysis pipeline utilised in analysing bronchial brush and bronchoalveolar lavage fluid from the Breathing Together Cohort groups 2-4. 

For context, the groups are as follows:

* **Group 2**: healthy preschool-aged controls
* **Group 3**: preschool children with severe wheeze
* **Group 4**: school-aged children diagnosed with asthma

## RNA sequencing data processing

Raw RNA sequencing FASTQ files were processed using the NF-CORE rnaseq pipeline (version 3.10.1). Briefly, reads underwent initial quality control (QC) with FastQC, UMIs were extracted with UMI-tools, adapters removed, and quality trimming performed with Trim Galore. Genomic contaminants and ribosomal components were removed with BBSplit and SortMeRNA respectively. Reads were aligned and quantified via a combination of STAR and Salmon tools, and sorted and indexed using SAMtools, and dereplicated with UMI-tools. An abundance matrix was produced, representing scaled and normalised transcripts per million reads, from which a re-estimated counts table was generated for downstream analysis. Gene data was filtered for a minimum of 100 reads per sampling group using the edgeR package (version 3.42.0) filterByExpr function, followed by voom normalisation with limma (version 3.56.0), and removal of non-protein-coding genes. Voom-normalised data matrix was saved for multi-omic integration.
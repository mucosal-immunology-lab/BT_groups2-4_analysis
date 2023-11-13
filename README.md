# Multi-omic data analysis of groups 2-4 of the Breathing Together study

Here we present the data processing and subsequent multi-omic integration and analysis pipeline utilised in analysing bronchial brush and bronchoalveolar lavage fluid from the Breathing Together Cohort groups 2-4. 

For context, the groups are as follows:

* **Group 2**: healthy preschool-aged controls
* **Group 3**: preschool children with severe wheeze
* **Group 4**: school-aged children diagnosed with asthma

- [Multi-omic data analysis of groups 2-4 of the Breathing Together study](#multi-omic-data-analysis-of-groups-2-4-of-the-breathing-together-study)
  - [Host RNA sequencing data processing](#host-rna-sequencing-data-processing)
  - [Bacterial DNA sequencing data processing](#bacterial-dna-sequencing-data-processing)


## Host RNA sequencing data processing

Raw RNA sequencing FASTQ files were processed using the [NF-CORE rnaseq pipeline](https://nf-co.re/rnaseq/3.12.0) (version 3.10.1) as described. Briefly, reads underwent initial quality control (QC) with FastQC, UMIs were extracted with UMI-tools, adapters removed, and quality trimming performed with Trim Galore. Genomic contaminants and ribosomal components were removed with BBSplit and SortMeRNA respectively. Reads were aligned and quantified via a combination of STAR and Salmon tools, and sorted and indexed using SAMtools, and dereplicated with UMI-tools. An abundance matrix was produced, representing scaled and normalised transcripts per million reads, from which a re-estimated counts table was generated for downstream analysis.

An RMarkdown file with the R-based data processing is provided [here](./01_HostTranscriptomics.Rmd).  Gene data was filtered for a minimum of 100 reads per sampling group using the [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) package (version 3.42.0) filterByExpr function, followed by voom normalisation with [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) (version 3.56.0), and removal of non-protein-coding genes. Voom-normalised data matrix was saved for multi-omic integration.

<img src="./assets/nf-core-rnaseq_metro_map_grey.png">

## Bacterial DNA sequencing data processing

Raw shotgun sequencing reads were processed as previously described using the [Sunbeam pipeline](https://sunbeam.readthedocs.io/en/stable/) for adaptor trimming, quality control, host genome decontamination, assembly of contiguous sequences, and co-assembly. Contig taxonomy was estimated using [Kraken 2](https://ccb.jhu.edu/software/kraken2/). 

Individual FASTQ files were aligned and mapped to the co-assembly using [bowtie2](https://github.com/BenLangmead/bowtie2), per-sample contig counts generated using a custom bash script, then combined into a counts matrix in R. 

To decontaminate the contigs, any sequence found in the extraction negative controls was removed. Using the ratio of reads before and after decontamination, samples with >90% non-contaminant reads were retained. A decontaminated co-assembly FASTQ file was generated from the remaining contigs. 

Functional KEGG orthology (KO) was assigned using the [GhostKOALA](https://www.kegg.jp/ghostkoala/) tool provided by the Kyoto Encyclopedia of Genes and Genomes, then agglomerated using taxonomy to produce taxa-KO pairs with their associated counts. Using the R [microbiome](https://bioconductor.org/packages/release/bioc/html/microbiome.html) package (version 1.22.0), the dataset was filtered using a detection threshold of 1 in at least 10% of samples via the core function, and centralised log-ratio (CLR)-normalised via the transform function. The final dataset was saved as a matrix for multi-omic data integration.
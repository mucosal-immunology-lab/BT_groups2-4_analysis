# Multi-omic data analysis of groups 2-4 of the Breathing Together study

Here we present the data processing and subsequent multi-omic integration and analysis pipeline utilised in analysing bronchial brush and bronchoalveolar lavage fluid from the Breathing Together Cohort groups 2-4. 

For context, the groups are as follows:

* **Group 2**: healthy preschool-aged controls
* **Group 3**: preschool children with severe wheeze
* **Group 4**: school-aged children diagnosed with asthma

- [Multi-omic data analysis of groups 2-4 of the Breathing Together study](#multi-omic-data-analysis-of-groups-2-4-of-the-breathing-together-study)
  - [Host RNA sequencing data processing](#host-rna-sequencing-data-processing)
    - [Cluster steps](#cluster-steps)
    - [R steps](#r-steps)
  - [Bacterial DNA sequencing data processing](#bacterial-dna-sequencing-data-processing)
    - [Cluster steps](#cluster-steps-1)


## Host RNA sequencing data processing

### Cluster steps

Raw RNA sequencing FASTQ files were processed using the [NF-CORE rnaseq pipeline](https://nf-co.re/rnaseq/3.12.0) (version 3.10.1) as described. Briefly, reads underwent initial quality control (QC) with FastQC, UMIs were extracted with UMI-tools, adapters removed, and quality trimming performed with Trim Galore. Genomic contaminants and ribosomal components were removed with BBSplit and SortMeRNA respectively. Reads were aligned and quantified via a combination of STAR and Salmon tools, and sorted and indexed using SAMtools, and dereplicated with UMI-tools. An abundance matrix was produced, representing scaled and normalised transcripts per million reads, from which a re-estimated counts table was generated for downstream analysis.

### R steps

An RMarkdown file with the R-based data processing is provided [here](./01_HostTranscriptomics.Rmd).  Gene data was filtered for a minimum of 100 reads per sampling group using the [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) package (version 3.42.0) filterByExpr function, followed by voom normalisation with [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) (version 3.56.0), and removal of non-protein-coding genes. Voom-normalised data matrix was saved for multi-omic integration.

<img src="./assets/nf-core-rnaseq_metro_map_grey.png">

## Bacterial DNA sequencing data processing

Raw shotgun sequencing reads were processed as previously described using the [Sunbeam pipeline](https://sunbeam.readthedocs.io/en/stable/) for adaptor trimming, quality control (QC), host genome decontamination, assembly of contiguous sequences, and co-assembly. Contig taxonomy was estimated using [Kraken 2](https://ccb.jhu.edu/software/kraken2/).

### Cluster steps

These steps were performed on the [M3 MASSIVE cluster](https://www.massive.org.au/). The following Sunbeam commands will run the previous steps including the QC prior to assembly and co-assembly.

```bash
# Assembly with Megahit - assembles reads to produce longer contiguous sequences
sunbeam run --configfile sunbeam_config.yml --cluster "sbatch --job-name=sunbeam_all_assembly --account=of33 --time=04:00:00 --mem-per-cpu=8G --ntasks=1 --cpus-per-task=20 --partition=genomics --qos=genomics" -j 8 -w 60 -p all_assembly --max-jobs-per-second 1 --keep-going

# Coassembly - produces a single file with all contiguous sequences
sunbeam run --configfile sunbeam_config.yml --cluster "sbatch --job-name=sunbeam_all_coassembly --account=of33 --time=04:00:00 --mem-per-cpu=8G --ntasks=1 --cpus-per-task=20 --partition=genomics --qos=genomics" -j 8 -w 60 -p --use-conda all_coassemble --rerun-incomplete --max-jobs-per-second 1 --keep-going

# Create a new folder for Anvi'o-related tasks
# Then copy the co-assembly file
mkdir -p anvio; cd anvio
cp ../sunbeam_output/assembly/coassembly/all_final_contigs.fa .

# Reformat fasta file and filter for contigs >500
anvi-script-reformat-fasta all_final_contigs.fa -o contigs-fixed.fa -l 500 --simplify-names
```

Individual FASTQ files were aligned and mapped to the co-assembly using a [custom script](./map_contigs.sh) built around [bowtie2](https://github.com/BenLangmead/bowtie2). Per-sample contig counts generated using a second [custom bash script](./contig_counts.sh).

```bash
# Build a bowtie2 index from the co-assembly
bowtie2-build contigs-fixed.fa assembly

# Run the contig mapping step for each of your samples
bash map_contigs.sh
```

The output from this is a set of `.txt` files with contig names (`RNAME`, i.e. read name) and counts that were saved to the local computer to be combined into a counts matrix in R.

```bash
count  RNAME
      2 c_000000000054
      2 c_000000000185
      7 c_000000000398
     90 c_000000000445
     28 c_000000000616
      2 c_000000000650
      2 c_000000000748
      2 c_000000000755
      4 c_000000000958
# etc.
```

The following steps are provided in an [RMarkdown file](./02_Metagenomics.Rmd). To decontaminate the contigs, any sequence found in the extraction negative controls was removed. Using the ratio of reads before and after decontamination, samples with >90% non-contaminant reads were retained. A decontaminated co-assembly FASTQ file was generated from the remaining contigs. This clean co-assembly was used to build an Anvi'o contigs database, from which the amino acid sequences for each gene call could be retrieved and used as input for functional KEGG orthology (KO) assignment using the [GhostKOALA](https://www.kegg.jp/ghostkoala/) tool provided by the Kyoto Encyclopedia of Genes and Genomes. The relevant steps for building the database and parsing GhostKOALA assignments can here found in more detail [here](https://github.com/mucosal-immunology-lab/microbiome-analysis/wiki/Anvio-pipeline).

The dataset was then agglomerated using taxonomy to produce taxa-KO pairs with their associated counts. Using the R [microbiome](https://bioconductor.org/packages/release/bioc/html/microbiome.html) package (version 1.22.0), the dataset was filtered using a detection threshold of 1 in at least 10% of samples via the core function, and centralised log-ratio (CLR)-normalised via the transform function. The final dataset was saved as a matrix for multi-omic data integration.
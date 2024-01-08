#!/bin/bash

# A simple loop to serially map all samples.

# How many threads should each mapping task use?
NUM_THREADS=8

# Create a directory for mapping results if it doesn't exist
mkdir -p mapping

# This assumes that your main Sunbeam folder is the parent folder 
# of your current "anvio" folder, and that you have not renamed 
# the samples.csv file created by sunbeam init
for sample in `awk -F, 'BEGIN {OFS=","} { print $1 }' ../samples.csv`
do
    if [ "$sample" == "sample" ]; then continue; fi
    
    # you need to make sure you "ls 01_QC/*QUALITY_PASSED_R1*" returns R1 files for all your samples in samples.csv
    R1s=`ls ../sunbeam_output/qc/decontam/$sample*_1* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    R2s=`ls ../sunbeam_output/qc/decontam/$sample*_2* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    
    echo "Running:" $R1s
    echo "Running:" $R2s
    
    bowtie2 --threads $NUM_THREADS -x ./assembly -1 $R1s -2 $R2s --no-unal -S mapping/$sample.sam
    samtools view -F 4 -bS mapping/$sample.sam > mapping/$sample-RAW.bam
    anvi-init-bam mapping/$sample-RAW.bam -o mapping/$sample.bam
    rm mapping/$sample.sam mapping/$sample-RAW.bam
done

# Mapping is done, and we no longer need bowtie2-build files
rm *.bt2

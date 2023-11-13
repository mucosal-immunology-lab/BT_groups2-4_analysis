#!/bin/bash

# Create a directory for mapping results if it doesn't exist
mkdir -p mapping/alignment_counts

# Generate contig counts
for FILE_NAME in `ls mapping/*.bam`
do
    FILE="${FILE_NAME##*/}"
    rm -f "mapping/alignment_counts/${FILE}_RNAME.txt"
    samtools view "mapping/${FILE}" > "mapping/alignment_counts/tmp.txt"
    echo "count  RNAME" >> "mapping/alignment_counts/${FILE}_RNAME.txt"
    cut -f 3 "mapping/alignment_counts/tmp.txt" | sort | uniq -c >> "mapping/alignment_counts/${FILE}_RNAME.txt"
    rm "mapping/alignment_counts/tmp.txt"
done
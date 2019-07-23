#!/bin/bash

# Merge multiple bam files from same individual using the short ID prefix from
# previous steps. This step uses an underscore to separate wildcard
# Output file will then be re-indexed

samples="/data2/redpolls/sample_name_unders-pt1.txt"
outdir="/data2/redpolls/merged_bam_files"

while read -r ID; do
echo "merging: "$ID
samtools merge "$outdir"$ID"merged.bam" $ID"_"*".bam"
samtools index "$outdir"$ID"merged.bam"

done<"$samples"

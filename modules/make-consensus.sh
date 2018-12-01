#! /bin/bash

inputvcf="/data2/PracticeGenomes/vcf_files/merged_vcf_files/minimum_set_filter3.recode.vcf.gz"
reference="/data2/PracticeGenomes/Gfortis_genome/GeoFor_1.0_genomic.fna"
samples="/data2/PracticeGenomes/cfc1.txt"
outdir="/data2/PracticeGenomes/consensus_seqs/"

while read -r samp; do
echo "writing consensus for sample" $samp

bcftools consensus --sample $samp -f $reference \
-o $outdir$samp"_consensus.fa" $inputvcf

done<"$samples"

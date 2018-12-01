#!/bin/bash

InputDir="/data2/rosyfinches/vcf_files/0p/"
OutDir="/data2/rosyfinches/vcf_files/0p/"
OutfilePrefix="0p_AUS"
SampleNames="/data2/rosyfinches/sample_names/sample_names_BCRF.txt"
InVCF="0p_AUS_snps_only-variants_biallelic.vcf"


while read -r SAMP; do

vcftools --vcf $InputDir$InVCF --indv $SAMP --freq --out $SAMP"_temp_allele_freqs"

awk 'FNR>1 && $3<3 {print $5}' $SAMP"_temp_allele_freqs.frq" > $SAMP"_temp_only_alleles"

sed s/[ACTG]:// $SAMP"_temp_only_alleles" | tr '\n' ',' | sed s/^/$SAMP","/ | sed s/','$/'\n'/ >> $OutDir$OutfilePrefix"_allele_freq_table.csv"

rm $SAMP"_temp_allele_freqs.frq"
rm $SAMP"_temp_allele_freqs.log"
rm $SAMP"_temp_only_alleles"

done<"$SampleNames"



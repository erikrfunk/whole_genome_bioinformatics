#!/bin/bash

# Call snps in samtools

ref="/data2/rosyfinches/HouseFinch/final.assembly.homo.fa"
bamdir="/data2/rosyfinches/sorted_bam_files/"

ID="rosyfinches" # This will be used as a prefix for the output file

echo "making a pileup file for" $ID >> log
#can also add the -R flag joined with a scaffold list to subset and parralel
bcftools mpileup -Ou -f $ref --ignore-RG -a AD,ADF,DP,SP,INFO/AD,INFO/ADF \
"$bamdir"*.bam | bcftools call -mv > "$ID"_raw_variants.vcf
#echo "removing all lines with two comment marks" >> log
#grep -v "##" "$ID"_snps_indels.vcf > "$ID"_snps_indels_short.vcf
echo "filtering low quality snps (<100)" >> log
awk '$1~/^#/ || $6 > 100 {print $0}' > \
"$ID"_snps_indels_filtered.vcf "$ID"_raw_variants.vcf
echo "add the header and check length of column 4 and 5 to make sure
they are snp type variants" >> log
awk '$1~/^#/ || length($4)==1 && length($5)==1 {print $0}'> \
"$ID"_snps_filtered.vcf "$ID"_snps_indels_filtered.vcf

# Call snps in samtools

ref="/data2/rosyfinches/HouseFinch/final.assembly.homo.fa"
samplenames="/data2/rosyfinches/sample_names.txt"
bamdir="/data2/rosyfinches/sorted_bam_files"

# while read -r ID; do
ID="rosyfinches"

echo "making a pileup file for" $ID >> log
samtools mpileup -uf $ref "$bamdir"/*_sorted_RGadded_dupmarked.bam | bcftools call -cv > "$ID"_snps_indels.vcf
echo "removing all lines with two comment marks" >> log
grep -v "##" "$ID"_snps_indels.vcf > "$ID"_snps_indels_short.vcf
echo "filtering low quality snps (<100)" >> log
awk '$6 > 100 {print $0}' > "$ID"_snps_indels_filtered.vcf "$ID"_snps_indels_short.vcf
echo "add the header and check length of column 4 and 5 to make sure they are snp type variants" >> log
awk 'FNR==1 || length($4)==1 && length($5)==1 {print $0}'> "$ID"_snps_filtered.vcf "$ID"_snps_indels_filtered.vcf

# done<"$samplenames"

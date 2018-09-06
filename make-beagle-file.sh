#Make a beagle-format file of genotype likelihoods for input to NGSadmix
#Input requires vcf file 

variants="/data2/rosyfinches/vcf_files/rosyfinches_all_genotyped_filtered_snps.vcf.gz"
outdir="/data2/rosyfinches/Admix"

#First create a list of all the scaffolds and then remove duplicates
#awk '$1!~/^[#]/ {print $1}' $variants > scaffold_list.txt
#awk '!seen[$0]++' scaffold_list.txt > final_scaffold_list.txt
rm scaffold_list.txt

echo "number of scaffolds:"
wc -l final_scaffold_list.txt

while read -r CHROM; do
vcftools --gzvcf $variants --out "$outdir"temp --BEAGLE-PL --chr $CHROM
cat "$outdir"temp.BEAGLE.PL >> "$outdir"rosyfinches_BEAGLElikes_whole_genome.PL
done < final_scaffold_list.txt
rm "$outdir"temp

#This renames all the scaffolds in the vcf file to be numerical
#awk '$1~/^[#]/ {print $0;} $1!~/^[#]/ {split($1,a,"_"); $1=a[2]; print $0}' $vcf > /data2/rosyfinches/snps_from_samtools/rosyfinches_snps_samtools_final_TEST_SCAFF_REFORM.vcf

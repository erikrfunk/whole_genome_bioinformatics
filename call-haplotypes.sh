# input name of sample and path to sorted, duplicate marked, RG added bam file

samplenames="/data2/rosyfinches/sample_names_batch1_subset1.txt"

while read -r sample; do
gatk HaplotypeCaller \
-R /data2/rosyfinches/HouseFinch/final.assembly.homo.fa \
-I /data2/rosyfinches/sorted_bam_files/"$sample"_sorted_RGadded_dupmarked.bam \
-O /data2/rosyfinches/vcf_files/"$sample".g.vcf \
-ERC GVCF
done<"$samplenames"


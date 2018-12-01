# Shell script to sort bam files according to a sample identifier 
# and prepare them for haplotype calling, including adding SM tag in RG, marking duplicates, and indexing

# This requries a txt file with a list of sample names you wish to sort bam files for (samplenames)

samplenames="/data2/rosyfinches/sample_names_batch2.txt"

while read -r SAMPLE; do

samtools sort -T temp -@ 6 \
-o /data2/rosyfinches/sorted_bam_files/"$SAMPLE"_sorted.bam \
/data2/rosyfinches/bam_files/$SAMPLE.bam

picard-tools AddOrReplaceReadGroups \
I=/data2/rosyfinches/sorted_bam_files/"$SAMPLE"_sorted.bam \
O=/data2/rosyfinches/sorted_bam_files/"$SAMPLE"_sorted_RGadded.bam \
RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$SAMPLE

picard-tools MarkDuplicates \
I=/data2/rosyfinches/sorted_bam_files/"$SAMPLE"_sorted_RGadded.bam \
O=/data2/rosyfinches/sorted_bam_files/"$SAMPLE"_sorted_RGadded_dupmarked.bam \
M=/data2/rosyfinches/sorted_bam_files/"$SAMPLE"_dupmarked_metrics.txt

samtools index \
/data2/rosyfinches/sorted_bam_files/"$SAMPLE"_sorted_RGadded_dupmarked.bam

rm /data2/rosyfinches/sorted_bam_files/"$SAMPLE"_sorted.bam
rm /data2/rosyfinches/sorted_bam_files/"$SAMPLE"_sorted_RGadded.bam

done<"$samplenames"

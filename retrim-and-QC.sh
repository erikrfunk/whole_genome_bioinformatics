# shell script for pre-trim QC, trimming, and post-trim QC
# pre and post trim QC files go to separate directories

filename="sample_names.txt"

echo "Re-trim and post-trim QC for"$filename>>retrim_and_QC_log.txt

while read -r ID; do
echo "Beginning trimming for "$ID>>retrim_and_QC_log.txt
TrimmomaticPE -threads 6 \
/data2/rosyfinches/fastqs/"$ID"_read1.fastq /data2/rosyfinches/fastqs/"$ID"_read2.fastq \
-baseout /data2/rosyfinches/fastqs/"$ID"_trimmed.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90>>retrim_and_QC_log.txt

echo "Beginning post-trim QC for "$ID>>retrim_and_QC_log.txt
fastqc -t 6 \
/data2/rosyfinches/fastqs/"$ID"_trimmed_1P.fq.gz \
/data2/rosyfinches/fastqs/"$ID"_trimmed_2P.fq.gz \
--outdir /data2/rosyfinches/post_trim_QC_files/
echo $ID" post-trim QC done">>retrim_and_QC_log.txt

done<"$filename"

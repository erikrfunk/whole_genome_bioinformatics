# Assemble reads with reference genome using bwa

ref="/data2/rosyfinches/HouseFinch/final.assembly.homo.fa"
seqs="/data2/rosyfinches/shell_scripts/sample_names.txt"

echo "Beginning alignment for " $seqs>>bwa_alignment_log.txt
echo "indexing reference">>bwa_alignment_log.txt

bwa index $ref

while read -r ID; do
echo "aligning " $ID >> bwa_alignment_log.txt

bwa mem -t 6 $ref /data2/rosyfinches/fastqs/"$ID"_trimmed_1P.fq.gz \
/data2/rosyfinches/fastqs/"$ID"_trimmed_2P.fq.gz | \
samtools view -b -o /data2/rosyfinches/bam_files/$ID.bam -S

echo "sam file piped into samtools view to convert to .bam">>bwa_alignment_log.txt

done<"$seqs"

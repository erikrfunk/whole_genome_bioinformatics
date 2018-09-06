# Assemble reads with reference genome using bwa

ref="/data2/rosyfinches/mt_scaffolds/Larctoa_complete_mitochondrion.fasta"
seqs="/data2/rosyfinches/sample_names/sample_names_batch1_short.txt"

echo "Beginning alignment for " $seqs>>bwa_alignment_log.txt
echo "indexing reference">>bwa_alignment_log.txt

bwa index $ref

while read -r ID; do
echo "aligning " $ID >> bwa_alignment_log.txt

bwa mem -t 12 $ref /data2/rosyfinches/fastqs/"$ID"_trimmed_1P.fq.gz /data2/rosyfinches/fastqs/"$ID"_trimmed_2P.fq.gz | samtools view -b -o /data2/rosyfinches/mt_scaffolds/mt_bam_files/$ID.bam -S

echo "sam file piped into samtools view to convert to .bam">>bwa_alignment_log.txt

done<"$seqs"


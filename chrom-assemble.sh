# Map consensus sequence scaffolds to chromosomes from reference individually

Ref="/data2/rosyfinches/mt_scaffolds/Larctoa_complete_mitochondrion.fasta"
PathToSeqs="/data2/rosyfinches/consensus_sequences/"
SampleNames="/data2/rosyfinches/sample_names/sample_names.txt"

bwa index $Ref

while read -r SAMPLE; do

bwa mem -t 6 $Ref $PathToSeqs$SAMPLE"_consensus.fa" | samtools view -b -o /data2/rosyfinches/consensus_sequences/chrom_assemblies/"$SAMPLE"_mt_genome.bam -S

done<"$SampleNames"


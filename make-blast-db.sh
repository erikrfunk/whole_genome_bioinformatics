samplelist="/data2/rosyfinches/sample_names.txt"
workingdir="/data2/rosyfinches/snps_from_samtools/consensus_sequences/"

while read -r SAMP; do

makeblastdb -in $workingdir$SAMP"_consensus.fa" -dbtype nucl

done<"$samplelist"

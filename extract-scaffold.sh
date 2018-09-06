samplenames="/data2/rosyfinches/whole_genome_scripts/sample_1.txt"
scaffoldnames="/data2/rosyfinches/whole_genome_scripts/scaffolds_to_extract.txt"
pathtoconsensus="/data2/rosyfinches/consensus_sequences/"
outputdir="/data2/rosyfinches/"

while read -r SAMPLE; do
samplefile=$pathtoconsensus$SAMPLE"_consensus.fa"
echo "extracting from sample" $SAMPLE

while read -r SCAFFOLD; do
echo "extracting scaffold_"$SCAFFOLD

location=$(echo ">scaffold_"$SCAFFOLD )
firstline=$(awk -v a="$location" '$0==a {print FNR}' $samplefile)
echo "starting at line" $firstline

nextscaffold=$(( SCAFFOLD+1 ))
echo "next scaffold is scaffold_"$nextscaffold

secondline=$(awk -v b="$nextscaffold" '$0==">scaffold_"b {print FNR}' $samplefile)
echo "ending before line" $secondline

awk -v a="$firstline" -v b=$secondline 'FNR>=a && FNR<b {print $0}' $samplefile >> $outputdir$SAMPLE"_scaffold"$SCAFFOLD".txt"

done<"$scaffoldnames"
done<"$samplenames"

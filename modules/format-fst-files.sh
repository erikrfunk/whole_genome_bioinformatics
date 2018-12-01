# Script to format files output from weir and Cockerham's Fst for use in qqman
# Provide the prefix for all the comparison files needed reformatting
# Script then reads in each file, removes NaN rows, and calls an R script to add snp column

comparisons="/data2/rosyfinches/whole_genome_scripts/comparisons_list.txt"

while read -r COMP; do
#FstFile="/data2/rosyfinches/snps_from_samtools/Fst_files/"$COMP".weir.fst"
#outdir="/data2/rosyfinches/snps_from_samtools/Fst_files/"

# Or for windowed, comment out above, and uncomment below
FstFile="/data2/rosyfinches/Fst_files/"$COMP".windowed.weir.fst"
outdir="/data2/rosyfinches/Fst_files/"
echo "Beginning formatting for" $COMP
echo "Removing rows with NaNs and renaming scaffolds to numerical values"
grep -v "nan" $FstFile | sed 's/scaffold_// ' > "$outdir"temp.fst

Rscript /data2/rosyfinches/whole_genome_scripts/manhattan_plot_format.R


echo "Writing file to" $outdir
mv "$outdir"temp2.fst "$outdir""$COMP"_20kb_formatted.windowed.weir.fst

done<"$comparisons"

rm "$outdir"temp.fst

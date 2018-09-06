inputvcf="/data2/rosyfinches/vcf_files/rosyfinches_all_genotyped_filtered_snps.vcf.gz"
reference="/data2/rosyfinches/HouseFinch/final.assembly.homo.fa"
samples="/data2/rosyfinches/sample_names.txt"
outdir="/data2/rosyfinches/consensus_sequences/"

while read -r samp; do
echo "writing consensus for sample" $samp

bcftools consensus --sample $samp -f $reference \
-o $outdir$samp"_consensus.fa" $inputvcf

done<"$samples"

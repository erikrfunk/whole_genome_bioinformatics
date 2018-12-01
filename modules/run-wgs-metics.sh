SampleNames="/data2/rosyfinches/sample_names/sample_names.txt"
PathToBams="/data2/rosyfinches/sorted_bam_files/"

while read -r SAMPLE; do

gatk CollectWgsMetrics -I $PathToBams$SAMPLE"_sorted_RGadded_dupmarked.bam" -O $SAMPLE"_wgs_metrics.txt" -R "/data2/rosyfinches/HouseFinch/final.assembly.homo.fa"

done<"$SampleNames"

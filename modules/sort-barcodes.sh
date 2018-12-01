# Shell script to demultiplex from fastq header

barcodes="/data2/rosyfinches/barcodes/additional_barcodes.txt"
forward="/data2/rosyfinches/Lane4/ROFI_Lane_4_TKD180302537-W-1_HL57NCCXY_L3_1.fq"
reverse="/data2/rosyfinches/Lane4/ROFI_Lane_4_TKD180302537-W-1_HL57NCCXY_L3_2.fq"

count=1
while read -r MATCH; do
echo "This line reads" $MATCH
bc=$(awk -v counter=$count 'NR==counter {print $1}' $barcodes) 
echo "Barcode is" $bc
sample_ID=$(awk -v counter=$count 'NR==counter {print $2}' $barcodes)
echo "Writing $bc to" $sample_ID"read.fastq"

grep --no-group-separator -A 3 "$bc" $forward >> "/data2/rosyfinches/fastqs/"$sample_ID"_read1.fastq"
grep -A 3 --no-group-separator "$bc" $reverse >> "/data2/rosyfinches/fastqs/"$sample_ID"_read2.fastq"

count=$[$count + 1]
done <"$barcodes"

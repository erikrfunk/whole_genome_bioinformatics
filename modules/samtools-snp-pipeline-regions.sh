#!/bin/bash

if [ $# -lt 1 ]
  then
    echo "Call variants using bcftools mpileup.

    [-i] Path to directory of sorted bam files
    [-r] Reference genome

    OPTIONAL ARGUMENTS

    [-g] Regions file to limit the variant call
    [-f] Logical to filter - if yes, type T. Filters pass only high quality snps
    [-o] Output prefix"

  else
    while getopts i:r:g:f:o: option
    do
    case "${option}"
    in
    i) bamdir=${OPTARG};;
    r) ref=${OPTARG};;
    g) regions=${OPTARG};;
    f) filter=${OPTARG};;
    o) ID=${OPTARG};;
    esac
    done

    bamdir="${bamdir:-sorted_bam_files/}"
    regions="${regions:-FALSE}"
    filter="${filter:-F}"
    ID="${ID:-output}"

    echo "making a pileup file for" $ID >> log
    # Check if a regions file is provided then call mpileup
    if [ $regions == FALSE ]
      then
        bcftools mpileup -Ou -f $ref--ignore-RG -a AD,ADF,DP,SP,INFO/AD,INFO/ADF \
        "$bamdir"*.bam | bcftools call -mv > "$ID"_raw_variants.vcf
      else
        bcftools mpileup -Ou -f $ref --ignore-RG --regions-file $regions \
        -a AD,ADF,DP,SP,INFO/AD,INFO/ADF \
        "$bamdir"*.bam | bcftools call -mv > "$ID"_raw_variants.vcf
    fi

    # Check if filtering is required then either filter or pass
    if [ $filter == T]
      then
        echo "filtering low quality snps (<100) for" $ID >> log
        awk '$1~/^#/ || $6 > 100 {print $0}' > \
        "$ID"_filtered_variants.vcf "$ID"_raw_variants.vcf
        echo "checking the length of column 4 and 5 to make sure
        they are snp type variants for " $ID >> log
        awk '$1~/^#/ || length($4)==1 && length($5)==1 {print $0}'> \
        "$ID"_filtered_snps.vcf "$ID"_filtered_variants.vcf
      else
        echo "Not filtering" $ID >> log
    fi
fi

# Call snps in samtools
# Erik's original script from which the above was modified
#ref="/data2/rosyfinches/HouseFinch/final.assembly.homo.fa"
#bamdir="/data2/rosyfinches/sorted_bam_files/"
#regions="/data5/meadowlarks/scaffolds1-2000.txt"
#ID="rosyfinches" # This will be used as a prefix for the output file

#echo "making a pileup file for" $ID >> log
#can also add the -R flag joined with a scaffold list to subset and parralel
#bcftools mpileup -Ou -f $ref--ignore-RG --regions-file $regions -a AD,ADF,DP,SP,INFO/AD,INFO/ADF \
#"$bamdir"*.bam | bcftools call -mv > "$ID"_raw_variants.vcf
#echo "removing all lines with two comment marks" >> log
#grep -v "##" "$ID"_snps_indels.vcf > "$ID"_snps_indels_short.vcf
#echo "filtering low quality snps (<100)" >> log
#awk '$1~/^#/ || $6 > 100 {print $0}' > \
#"$ID"_snps_indels_filtered.vcf "$ID"_snps_indels_short.vcf
#echo "add the header and check length of column 4 and 5 to make sure
#they are snp type variants" >> log
#awk '$1~/^#/ || length($4)==1 && length($5)==1 {print $0}'> \
#"$ID"_snps_filtered.vcf "$ID"_snps_indels_filtered.vcf

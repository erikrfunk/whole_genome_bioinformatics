#!/bin/bash

#Merge gvcfs, genotype, and filter variants

if [ $# -lt 1 ]
  then
    echo "Merge individual gvcf files, genotype, and filter variants.
    Default filter settings come from GATK best practices.

    [-i] Sample list
    [-r] Reference genome
    [-p] Path to gvcfs
    [-o] Name of the cohort

    OPTIONAL ARGUMENTS

    [-q] QD filter - default is < 2.0
    [-f] FS filter - default is > 60.0
    [-m] MQ filter - default is < 40.0
    [-s] MQRankSum - default is < -12.5
    [-e] ReadPosRankSum - default is < -8.0
    [-n] Filter set name - default is 'GATK_recommended_filters'"

  else
    while getopts i:r:p:o:q:f:m:s:e:n: option
    do
    case "${option}"
    in
    i) names=${OPTARG};;
    r) ref=${OPTARG};;
    p) gvcf_path=${OPTARG};;
    o) cohort=${OPTARG};;
    q) qd=${OPTARG};;
    f) fs=${OPTARG};;
    m) mq=${OPTARG};;
    s) mqranksum=${OPTARG};;
    e) readposranksum=${OPTARG};;
    n) filter_set_name=${OPTARG};;

    esac
    done

    gvcf_path="${gvcf_path:-vcf_files/}"
    qd="${qd:-< 2.0}"
    fs="${fs:-> 60.0}"
    mq="${mq:-< 40.0}"
    mqranksum="${mqranksum:-< -12.5}"
    readposranksum="${readposrank:-< -8.0}"
    filter_set_name="${filter_set_name:-GATK_recommended_filters}"

    files=""
    while read -r sample; do
    files="$files-V $gvcf_path$sample.g.vcf.gz "
    done<"$names"

    echo "Combining gvcf files."
    gatk CombineGVCFs \
    -R $ref \
    $files \
    -O $gvcf_path$cohort.g.vcf.gz

    echo "Genotyping cohort vcf"
    gatk GenotypeGVCFs \
    -R $ref \
    -V $gvcf_path$cohort.g.vcf.gz \
    -O $cohort.vcf.gz \

    echo "Applying hard filters"
    gatk VariantFiltration \
    -R $ref \
    -V $cohort.vcf.gz \
    -O "$cohort"_filtered.vcf.gz \
    --filter-expression "QD $qd || FS $fs || MQ $mq || MQRankSum $mqranksum || ReadPosRankSum $readposranksum" \
    --filter-name $filter_set_name

fi

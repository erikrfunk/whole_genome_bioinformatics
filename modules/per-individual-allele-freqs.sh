#!/bin/bash

if [ $# -lt 1 ]
  then
    echo "Calculates per-individual allele frequencies using VCFtools.

    [-i] Input VCF file
    [-o] Output file prefix
    [-s] Sample names - text file with one sample per line "

  else
    while getopts i:o:s: option
    do
    case "${option}"
    in
    i) input_vcf=${OPTARG};;
    o) output_prefix=${OPTARG};;
    s) sample_names=${OPTARG};;
    esac
    done

    while read -r SAMP; do

    vcftools --vcf $input_vcf --indv $SAMP --freq --out $SAMP"_temp_allele_freqs"

    awk 'FNR>1 && $3<3 {print $5}' $SAMP"_temp_allele_freqs.frq" > $SAMP"_temp_only_alleles"

    sed s/[ACTG]:// $SAMP"_temp_only_alleles" | tr '\n' ',' | sed s/^/$SAMP","/ | sed s/','$/'\n'/ >> $output_prefix"_allele_freq_table.csv"

    rm $SAMP"_temp_allele_freqs.frq"
    rm $SAMP"_temp_allele_freqs.log"
    rm $SAMP"_temp_only_alleles"

    done<"$sample_names"
fi

#!/bin/bash

if [ $# -lt 1 ]
  then
    echo "Divide a vcf.gz file into individual scaffolds using vcf tools
    A new directory will be created to house these
    This will loop through a list of scaffolds that should be separated out into
    individual vcf files.

    [-i] input VCF file - needs to be gzipped
    [-s] list of scaffolds - one per line as they appear in the vcf file

    OPTIONAL ARGUMENTS
    [-o] Output directory
    [-r] Switch to indicate the provided scaffolds should instead be removed
         and only the remaining scaffolds will be written (to a single vcf)"
  else
    while getopts i:s:ro: option
    do
    case "${option}"
    in
    i) vcf=${OPTARG};;
    s) scaffolds=${OPTARG};;
    r) reverse=1;;
    o) outdir=${OPTARG};;
    esac
    done

    outdir="${outdir:-individual_vcfs/}"

    if [$outdir == individuals_vcfs/]
      then
        mkdir individuals_vcfs
    fi

    if [$reverse == 1]
      then
        chr_arg = ""
        while read -r ID; do
          chr_arg += "--not-chr $ID "
        done<$scaffolds
        vcftools --gzvcf $vcf --out $outdir"remaining-scaffolds" $chr_arg \
        --recode

      else
        while read -r ID; do
          vcftools --gzvcf $vcf --out $outdir$ID --chr $ID --recode
        done<$scaffolds
    fi
fi

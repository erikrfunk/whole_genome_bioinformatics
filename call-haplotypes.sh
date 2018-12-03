#!/bin/bash

# calls variants using gatk4 HaplotypeCaller

if [ $# -lt 1 ]
  then
    echo "Calls variants using GATK 4 HaplotypeCaller in the -ERC GVCF mode.

    [-i] Sample list
    [-r] Reference genome
    [-p] Path to sorted bam files - If multiple bam files per individual are in
         this directory (i.e. a sorted and unsorted) a bam suffix needs to be
         provided using the '-s' flag (see optional arguments below)

    OPTIONAL ARGUMENTS

    [-o] Output directory
    [-s] Suffix of sorted bam file - default is '.bam'"

  else
    while getopts i:r:p:o:s: option
    do
    case "${option}"
    in
    i) names=${OPTARG};;
    r) ref=${OPTARG};;
    p) bam_path=${OPTARG};;
    o) outdir=${OPTARG};;
    s) bam_suffix=${OPTARG};;

    esac
    done

  outdir="${outdir:-vcf_files/}"
  bam_path="${bam_path:-sorted_bam_files/}"
  bam_suffix="${bam_suffix:-.bam}"
  if [ ! -d $outdir ]
    then
      mkdir vcf_files/
  fi

  echo "Creating sequence dictionary for HaplotypeCaller"
  prefix=$(cat <<< $ref | sed -r 's/^(.*)\.\w+$/\1/')
  picard-tools CreateSequenceDictionary \
  R=$ref \
  O=$prefix.dict

  echo "Indexing reference"
  samtools faidx $ref

  echo "Calling Haplotypes"
  while read -r sample; do
  gatk HaplotypeCaller \
  -R $ref \
  -I $bam_path"$sample"*"$bam_suffix" \
  -O $outdir"$sample".g.vcf.gz \
  -ERC GVCF
  done<"$names"

fi

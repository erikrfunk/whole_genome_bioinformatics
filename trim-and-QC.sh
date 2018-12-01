#!/bin/bash

# shell script for pre-trim QC, trimming, and post-trim QC
# pre and post trim QC files go to separate directories

if [ $# -lt 1 ]
  then
    echo "Runs a pre- and post-trim qc using fastqc for samples provided by an input file. Both a pre- and post-trim directory will be created.
    Then trims sequences using TrimmomaticPE. Changes to trimming parameters should be made if necessary.
    Following trimming, a post-trim qc will also be run using fastqc.

    [-i] Sample list - sample name should match the prefix of the fastq files
    [-f] Path to fastq files - should contain files in the format 'SAMPLENAME_read1.fastq' and 'SAMPLENAME_read2.fastq'

    OPTIONAL ARGUMENTS

    [-a] File with adapter sequences for Trimmomatic. Defaults to the TruSeq3.fa file in this repository"

  else
    while getopts i:f:a: option
    do
    case "${option}"
    in
    i) filename=${OPTARG};;
    f) fastqs_path=${OPTARG};;
    a) adapters=${OPTARG};;
    esac
    done

    adapters="${adapters:-TruSeq3-PE.fa}"

    echo "Trim and QC for "$filename>>trim_and_QC_log.txt
    mkdir pre_trim_QC_files
    mkdir post_trim_QC_files

    while read -r ID; do
    echo "Beginning pre-trim QC for "$ID>>trim_and_QC_log.txt
    fastqc -t 6 \
    $fastqs_path"$ID"_read1.fastq \
    $fastqs_path"$ID"_read2.fastq \
    --outdir pre_trim_QC_files/
    echo $ID" pre-trim QC done" >> trim_and_QC_log.txt

    echo "Beginning trimming for "$ID>>trim_and_QC_log.txt
    TrimmomaticPE -threads 6 \
    $fastqs_path"$ID"_read1.fastq $fastqs_path"$ID"_read2.fastq \
    -baseout $fastqs_path"$ID"_trimmed.fq.gz \
    ILLUMINACLIP:$adapters:1:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90>>trim_and_QC_log.txt

    echo "Beginning post-trim for "$ID>>trim_and_QC_log.txt
    fastqc -t 6 \
    $fastqs_path"$ID"_trimmed_1P.fq.gz \
    $fastqs_path"$ID"_trimmed_2P.fq.gz \
    --outdir post_trim_QC_files/
    echo $ID" post-trim QC done">>trim_and_QC_log.txt

    done<"$filename"
fi

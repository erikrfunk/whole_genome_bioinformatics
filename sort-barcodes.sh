#!/bin/bash

# Shell script to demultiplex from illumina pair-end fastq header

if [ $# -lt 1 ]
  then
    echo "Demultiplex from illumina paired-end fastq header.
    [-f] forward read file
    [-r] reverse read file
    [-b] barcode file - this file should contain the barcodes in column 1, as they appear
    in the fastq header, and the corresponding sample ID in column 2. See barcodes.txt

    OPTIONAL ARGUMENTS

    [-o] output directory - default will make new directory called 'fastqs' in current directory"

  else
    while getopts f:r:b:o: option
    do
    case "${option}"
    in
    f) forward=${OPTARG};;
    r) reverse=${OPTARG};;
    b) barcodes=${OPTARG};;
    o) outdir=${OPTARG};;
    esac
    done

    outdir="${outdir:-fastqs/}"
    if [ $outdir == fastqs/ ]
    then
      mkdir fastqs
    fi

    count=1
    while read -r MATCH; do
    echo "This line reads" $MATCH
    bc=$(awk -v counter=$count 'NR==counter {print $1}' $barcodes)
    echo "Barcode is" $bc
    sample_ID=$(awk -v counter=$count 'NR==counter {print $2}' $barcodes)
    echo "Writing $bc to" $sample_ID"read.fastq"

    grep --no-group-separator -A 3 "$bc" $forward >> $outdir$sample_ID"_read1.fastq"
    grep -A 3 --no-group-separator "$bc" $reverse >> $outdir$sample_ID"_read2.fastq"

    count=$[$count + 1]
    done <"$barcodes"

fi

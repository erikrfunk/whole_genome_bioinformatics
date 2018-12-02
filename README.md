This repository contains the shell scripts used for whole genome sequence analysis from demultiplexing to variant calls as well as a variety of other analyses related to population and phylogenomics.

Notes on the first iteration and a description of the pipeline on the rosy-finch data is located in rosyfinch_notes.md

Note: Most all code here is invoked in shell scripts that have not yet been formatted to take input directly. Instead, input files and directories need to be changed in the bash script - sorry
The call of any script without arguments will show a short description and the arguments list.

The workflow is broken up into ## sections, with each breakpoint used for validation or to provide modularity.

Required programs:
fastqc
Trimmomatic


------------------------------------------------------------------------------------------------------

**Step 1:** Demultiplexing \
This step can be done with any freely available demultiplexing program you prefer. The script provided here is a simple approach based on a grep search, though could be modified fairly easily to allow mismatches. Default is to make a new directory called "fastqs" in the current directory.

    sort-barcodes.sh -f -r -b

ARGUMENTS \
[-f] forward read file \
[-r] reverse read file \
[-b] barcode file - this file should contain the barcodes in column 1, as    they appear \
     in the fastq header, and the corresponding sample ID in column 2. See barcodes.txt

OPTIONAL ARGUMENTS \
[-o] output directory - default will make new directory called 'fastqs' in current directory

**Step 2:** Trim and Quality Control \
This step trims sequences using TrimmomaticPE, and performs quality control using fastqc on fastq files both before and after trimming. This step requires a text file with all sample names wanting to be included. These sample names should match the prefixes of the fastq files. This script also requires the path to the fastq files is provided.

Two new directories will be created, one for each quality control run. This allows the user to evaluate the effectiveness of trimming prior to moving on to the next step. An optional argument can be used to pass adapter sequences to trimmomatic. This will default to the TruSeq3.txt file, so if no arguments are passed using this flag, the TruSeq3 file needs to be in the directory. At this time, any additional adjustments to the trimmomatic settings or trim settings need to be made manually in the script, in the code block containing the TrimmomaticPE function call.

    trim-and-QC.sh -i -f

ARGUMENTS \
[-i] Sample list \
[-f] Path to fastq files - should contain files in the format 'SAMPLENAME_read1.fastq' and 'SAMPLENAME_read2.fastq'

OPTIONAL ARGUMENTS \
[-a] File with adapter sequences for Trimmomatic

**Step 3:** Assembly and Prep for Variant Calling \
This step aligns fastq reads to a reference genome using bwa mem.
Bam file is then sorted, duplicates are marked, and file is indexed using
samtools and picard-tools. After sorting, the sample name is added as a RG (Read Group) tag. This seems to be necessary if wanting to call variants with GATK later on. The default operation for this program is to retrieve fastq files from a directory, which defaults to the 'fastqs' directory created previously, so this should be called from the parent directory. This may change in the future but for now will stay as this. The directory for both the bam and sorted bam files can be specified, but see below.

    align-and-sort.sh -i -r

ARGUMENTS
[-i] Sample list
[-r] Reference genome

OPTIONAL ARGUMENTS

[-t] Number of threads to use
[-p] Path to trimmed fastqs - the default is a directory called 'fastqs' as
     produced from the initial sorting
[-b] Output directory for bam files - default is to make a directory
     called 'bam_files'
[-s] Output directory for sorted bam files - default is to make a
     directory called 'sorted_bam_files'

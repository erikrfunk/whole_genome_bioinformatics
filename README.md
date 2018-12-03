This repository contains the shell scripts used for whole genome sequence analysis from demultiplexing to variant calls as well as a variety of other analyses related to population and phylogenomics.

Notes on the first iteration and a description of the pipeline on the rosy-finch data is located in rosyfinch_notes.md

The call of any script without arguments will show a short description and the arguments list.

The workflow is broken up into 5 sections, with each breakpoint used for validation or to provide modularity.

Required programs:
fastqc
Trimmomatic
bwa
samtools
picard-tools
GATK v4

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

ARGUMENTS \
[-i] Sample list \
[-r] Reference genome

OPTIONAL ARGUMENTS \
[-t] Number of threads to use \
[-p] Path to trimmed fastqs - the default is a directory called 'fastqs' as
     produced from the initial sorting \
[-b] Output directory for bam files - default is to make a directory
     called 'bam_files' \
[-s] Output directory for sorted bam files - default is to make a
     directory called 'sorted_bam_files'

**Step 4:** Call Variants \
This step uses HaplotypeCaller from GATK ver.4 in -ERC GVCF mode. This step is isolated from the rest of the merging and filtering of variants due to the time constraint HaplotypeCaller may impose. Becuase HaplotypeCaller call variants for one individual at a time, is has been necessary for my own analyses to manually parallelize this process by splitting up the sample list and running multiple jobs at once. By separating this step out, variants for all individuals can be called faster. Whenever done, the complete sample list can be used again to merge all g.vcf files and filter in the following step. This step begins by preparing the reference sequence dictionary and index for gatk. For the time being, this step requires the input of the path to the sorted bam files. If multiple bam files per individual are in this directory (i.e. a sorted and unsorted), a bam suffix needs to be provided using the '-s' flag.

    call-haplotypes.sh -i -r -p

ARGUMENTS \
[-i] Sample list - may be truncated to include only some individuals \
[-r] Reference genome \
[-p] Path to sorted bam files - this defauls to 'sorted_bam_files' created in the previous step

OPTIONAL ARGUMENTS \
[-o] Output directory for vcf files - default is to make a directory called 'vcf_files' \
[-s] Suffix of sorted bam - default is '.bam'

**Step 5:** Merge, Genotype, and Filtering \
This step merges the gvcfs from HaplotypeCaller into a cohort gvcf using GATK CombineGVCFs. This file is then genotyped with GenotypeGVCFs, and finally filtered with VariantFiltration. The path to the gvcf files is required, but defaults to the vcf_files directory created previously. The cohort gvcf is placed in the vcf_files directory as well, but the final genotyped, filtered cohort vcf is placed in the parent directory. For filtration, each of the standard filters can be set, but the gatk recommended filters are all set as defaults. If a different filter expression is desired, this should be run separately after as an individual filter step (see gatk4 user guide for VariantFiltration). Depth and minor allele frequency should be filtered for manually depending on the dataset. Finally, a name for the filter set is required, but defaults to 'GATK_recommended_filters'.

    merge-gvcfs.sh -i -r -p -o

ARGUMENTS \
[-i] Sample list \
[-r] Reference genome \
[-p] Path to gvcfs \
[-o] Name of the cohort

OPTIONAL ARGUMENTS \
[-q] QD filter - default is < 2.0 \
[-f] FS filter - default is > 60.0 \
[-m] MQ filter - default is < 40.0 \
[-s] MQRankSum - default is < -12.5 \
[-e] ReadPosRankSum - default is < -8.0 \
[-n] Filter set name - default is 'GATK_recommended_filters'

**Final Considerations** \
The pipeline in its current state ends here with the first set of variant filters. The next step will almost certainly be to apply additional filters depending on the dataset and its intended use.

This may include filtering for variant type (SNP vs INDEL), depth, or minor allele frequency. Many different programs can carry out these filters, including gatk as used above. See the documentation on 'VariantFiltration', 'SelectVariants', and the manual on the package 'vcftools'.

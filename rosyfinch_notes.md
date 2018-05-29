# Comments on analysis of Rosy-Finch genome sequences
-------------------------------------------------------

# 1. Demultiplexing using custom bash script

/data2/rosyfinches/whole_genome_scripts/sort_barcodes.sh

Example input: ROFI_Lane_1_TKD180302534-1_HL57NCCXY_L5_1.fq

This essentially greps the sequence header for the barcode (given through an input file) and writes it to the fastq files corresponding to the sample name (given in column two of input file)

Four barcodes still undetermined (one a probable match with single bp mismatch)

# 2. Trim and QC 

performing pre-trim QC, trim, and post-trim QC using fastqc and trimmomatic in custom shell script

/data2/rosyfinches/whole_genome_scripts/trim-and-QC.sh

progress saved to

/data2/rosyfinches/whole_genome_scripts/trim_and_sort_QC_log.txt

used TruSeq3-PE.fa sequences for adapters to trim \
added the sequences from Winnie for this library prep

after trimming, QC showed continued levels of Nextera Transposase sequence contamination. \
Adding this nextera transposase sequences to the TruSeq3-PE file
	
	>Transposase_1
	AGATGTGTATAAGAGACAG

No difference. Trying these:
	
	>Trans1
	TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
	>Trans1_rc
	CTGTCTCTTATACACATCTGACGCTGCCGACGA
	>Trans2
	GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
	>Trans2_rc
	CTGTCTCTTATACACATCTCCGAGCCCACGAGAC

Also changed allowed mismatch in adapter to 2. QC after both of these changes looked better and removed adapter contamination 

May think about a second set with softer filter parameters. \
RF1_read1 pre trim = 14,473,504 total sequences at 150bp length \
RF1_trimmed_1P = 11,192,354 total sequences at 90-150bp length

# 3. Aligning reads

For this pass, I am aligning reads to the House Finch genome using bwa \
Processing these using the shell script

/data2/rosyfinches/whole_genome_scripts/bwa-assemble.sh

# 4. Processing alignments for variant calling

For this step, bam files are sorted, SM field is added as sample name in RG field of header, duplicates are marked, and final file is indexed \
Processed using the shell script

/data2/rosyfinches/whole_genome_scripts/merge-bam-files.sh

Processing these in two batches \
Batch1: \
RF1-35 (excluding RF7 which is still missing data) \
Batch2: \
RF36-68

# 5. Variant calling 

Performing variant calling using gatk HaplotypeCaller

first creating a sequence dictionary for the reference

	picard-tools CreateSequenceDictionary R=/data2/rosyfinches/HouseFinch/final.assembly.homo.fa \
	O=/data2/rosyfinches/HouseFinch/final.assembly.homo.dict

then index the reference fasta

	samtools faidx /data2/rosyfinches/HouseFinch/final.assembly.homo.fa

Finally, processing samples in loops of five using shell script \
/data2/rosyfinches/whole_genome_scripts/call-haplotypes.sh \
with samples sorted by

/data2/rosyfinches/sample_names_batch1_subset1.txt

# Comments on analysis of Rosy-Finch genome sequences
-------------------------------------------------------

# 1. Demultiplexing using custom bash script

	sort_barcodes.sh

Example input: ROFI_Lane_1_TKD180302534-1_HL57NCCXY_L5_1.fq

This essentially greps the sequence header for the barcode (given through an input file) and writes it to the fastq files corresponding to the sample name (given in column two of input file)

Four barcodes still undetermined (one a probable match with single bp mismatch)

# 2. Trim and QC

performing pre-trim QC, trim, and post-trim QC using fastqc and trimmomatic in custom shell script

	trim-and-QC.sh

progress saved to trim_and_sort_QC_log.txt

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

	bwa-assemble.sh

# 4. Processing alignments for variant calling

For this step, bam files are sorted, SM field is added as sample name in RG field of header, duplicates are marked, and final file is indexed \
Processed using the shell script

	merge-bam-files.sh

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
	call-haplotypes.sh
with samples sorted by

/data2/rosyfinches/sample_names_batch1_subset1.txt

Merging gvcfs from all samples into cohort gvcf

	gatk CombineGVCFs \
	-R reference.fasta \
	-V SAMPLE1.g.vcf.gz \
	-V SAMPLES2.g.vcf.gz \
	-V ... \
	-O cohort_name.g.vcf.gz

Genotyping the merged gvcf

	gatk GenotypeGVCFs \
	-R reference.fasta \
	-V MERGED.g.vcf.gz \
	-O OUTPUT.vcf.gz

# 6. Filtering Variants

Filtering variants using hard cutoffs recommended by BroadInstitute

	gatk VariantFiltration \
	-R reference.fasta \
	-V input.vcf.gz \
	-O output.vcf.gz \
	--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filter-name "GATK_recommended_filters"

Note: The decimals are required so that values are recognized as floats!

	gatk SelectVariants \
	-R reference.fasta \
	-V input.vcf.gz \
	--select-type-to-include SNP \
	--exclude-filtered \
	--exclude-non-variants

Get total number of variants using

	gatk CountVariants \
	-V input.vcf.gz

For this dataset 50249125 snps were retained

Should also look at the histogram of per site depth \
Want to filter out low and high tails that might be copy number variants

Can using gatk VariantFiltration again to filter variants that might fall within these tails\
Going to use depth cutoffs of 2 and 8 based on mean coverage

	gatk VariantFiltration \
	-R reference.fasta \
	-V input.vcf.gz \
	-O output.vcf.gz \
	--genotype-filter-expression "DP < 2 || DP > 8" \
	--genotype-filter-name "added_depth_filters"

May need to re-run SelectVariants with the --exclude-filtered flag set
Total number of variants did not change, but this may be becuase the filter is genotype level, not variant across all samples.

Also tried this same filter with vcftools --minDP and --maxDP. This also returned the same number of variants. Should verify this worked by calculating a depth by genotype table using vcftools.

	vcftools --gzvcf INPUT \
	--geno-depth

What this filter seem to be doing is setting the genotype to ./., but retaining all other information about the sample variant, so the high depth numbers remain on the vcf, and thus will continue to be counted as part of the total number of variants. This may also skew depth statistics, but this I am not sure about

Ultimately, I am going to make two datasets, one with the above filters, and one where only 25 percent missing data is allowed. See optional filters below.

Additional optional filters for minor allele frequence less than 0.1 \
Allow for 25% missing data by site  

	vcftools --gzvcf INPUT --out OUTPREFIX \
	--maf 0.1 \
	--max-missing 0.75 \
	--recode # Writes a new vcf

New snp counts for all genotypes and for 75p about 6million and about 100k \
Quite low for the 75p set but will check analyses anyway

# 7. Phylogenetic Analyses

Selecting a random fraction of snps for use in SVDquartets. \
Note, this may not be necessary depending on the number of variants left after filtering.

	gatk SelectVariants \
	-V input.vcf.gz \
	--select-random-fraction PERCENT

Verify the number of snps by rerunning the CountVariants command from above

Convert the vcf to a phylip using

	python vcf2phylip.py -i FILENAME -m MIN_SAMPLES_LOCUS -o OUTGROUP -f writes FASTA -n writes NEXUS

Constructed phylogenetic tree using SVDquartets, run in paup command line

	./Paup \
	SVDQuartets evalQuartets=all nquartets=20000 bootstrap=standard \
	nreps=500 nthreads=12 \
	treeFile=BOOTTREES_OUTPUT.tre

I then like to just open up the bootstrap trees on the desktop version to make my consensus tree, but this can be done command line (I think this is right)

	contree /majrul cutoff=50 le50

# 8. Extracting Gene Sequences

Re-ran whole pipeline using a complete mitochonrion as reference \
Then constructed a blastn database of all the consensus mt sequences

make-blastdb-and-seq-extract.sh

Ran this using the L. arctoa reference sequence for cyt b

Then make this new file a searchable database

	makeblastdb -in INFILE -parse_seqids -dbtype nucl

Then run blastdbcmd using the previously constructed seq_ids

	blastdbcmd -db NEW_DATABASE -dbtype nucl -entry_batch SEQ_IDS -out OUTFILE -outfmt %f


# 9. Make per indiviual allele frequencies

vcftools --vcf INPUT --freq --out OUTPUT

Format this file for input into R using

	per-individual-allele-freqs.sh


# 10. Assess introgression/ILS with dstats

The following sets up files for a D statistic that averages allele frequencies across populations of individuals \
First, filter the vcf for only biallelic snps

	vcftools --vcf INPUT.vcf --out OUT-PREFIX --min-alleles 2 --max-alleles 2 --recode

Convert this to a simplified "Genotype" vcf file using Simon Martin's parseVCF.py code

	python ~/genomics_general/VCF_processing/parseVCF.py -i input.vcf.gz --many more optional arguments | gzip > output.geno.gz

From this geno file, calculate frequencies. Again see Simon Martin's github for details on the arguments and file formats

	python ~/genomics_general/freq.py -g input.gz -p -p -p --popsFile --target derived -o output

Finally, take this frequency output and run the abba_baba.R script

Or, can use an individual input and use the Dfoil approach from John Pease

forth coming bash script for this

Can also use the ABBABABAwindows.py script to calculate D, and Fd across windows.\
Its recommended that you keep window size large enough to include about 100 snps. See the 'sites used' column to ensure this.

#Notes on rosy-finch snps using the samtools pipeline
-------------------------------------------------------

# 1. Demultiplexed using the same script as standard pipeline notes

# 2. Aligned, indexed, sorted, called variants, and filtered snps using *samtools-snp-pipeline.sh \
Aligned to the House Finch genome \
Counted snps that passed filtering using

	wc -l rosyfinces_snps_filtered.vcf

Reported 41573240 snps

Consider zipping at some point in this pipeline, or at least the final output file

# 3. PCA of snps

still working on installation of R packages onto server 

Currently working with SNPRelate \
Transferred the filtered vcf file to personal machine

When converting the vcf to GDS format, an error was thrown: Fewer columns than expected \
This error occurs for the same snp as described below (Scaffold_118, 1016359)\
The following line appears to code for the same snp, but has more columns (more samples) \
I think this may have come from the program getting interrupted, and restarting oon this snp \
Going to try deleting the line 

Made a copy of the vcf file: *rosyfinches_snps_filtered_copy2.vcf* \
Removed the line this snp occurs on: 22211484

	sed -e '22211484d' rosyfinches_snps_filtered_copy2.vcf > rosyfinches_snps_filtered_copy2_error_snp_removed.vcf 

Continuing to test file using vcf2gprobs in step 4 below. See there for details \
File was fixed and succfully ran through SNPRelate. R script on personal computer
\
Final snp count was ~20.9 million
 
Then created a second vcf file with just NA taxa to re-run PCA 

	vcftools --vcf FILE --out OUTPUT PREFIX \
	--keep FILE WITH SAMPLE NAMES \
	--recode # This is necessary to get a new vcf!

Final snp count was 15,538,763

saved as rosyfinches_snps_samtools_final_NAONLY.vcf

# 4. NGSadmix

Converting .vcf file from samtools to beagle format

	cat rosyfinches_snps_filtered.vcf | java -jar ~/BeagleUtilities/vcf2beagle.jar - rosyfinch_samtools_snps 

vcf2beagle requires the input file called via cat 
	
	cat [input_file]

with the output piped into the java call. This includes two parameters that need to be specified

	vcf2beagle.jar [missing] [prefix] \
where \
	missing = missing allele code in output Beagle version 3 genotypes file
	prefix = prefix for Beagle version 3 output files

Failed. snp at Scaffold_118, 1016359: Expected 76 fields but only found 49\
Trying the program vcf2gprobs.jar: Same error \
Tried a vcf to fasta conversion and it threw an out of range error. \
Removed this line (see step 3 above) and new error thrown: Expected 76 fields, but found 127

It appears two lines were merged including scaffold_118, 1017059 and scaffold_150, 636676 \
Need to instert a line break between these two. Managed to do this using vim \

Checking that all remaining lines contain 76 columns:

	awk 'NF!=76 {print FNR}' rosyfinches_snps_filtered_copy2_error_snp_removed.vcf > error_counts.txt

Removing returned line using: 

	sed -e '22211485d' rosyfinches_snps_filtered_copy2_error_snp_removed.vcf > rosyfinches_snps_samtools_final.vcf

Re-ran the awk check and no lines were returned, running vcf2gprobs again. \
Finished, but program only wrote headers, no probabilites. \
This command only works on the GL field, not the log-scaled PL field \

Using vcftools to format input file with GL field \
	vcftools [--vcf FILE | --gzvcf FILE | --bcf FILE][--out PREFIX][FILTERING OPTIONS][OUTPUT OPTIONS]

Calling as: \
	vcftools --vcf rosyfinches_snps_samtools_final.vcf --out rosyfinches_snps_PL_samtools --BEAGLE-PL --chr CHROMS

Throws an error. In filtering step after mpileup, I removed the metadata which defined the FORMAT fields \
I think this is required for the conversion. \
Added the PL and GT format fields back into the metadata header section.\
Retrying. Same error. Same error. SAME ERROR.
Put --chr in quotes and only inlcuded one scaffold and it worked \
Its possible multiple scaffolds could be inlcuded if scaffold name was numberic and not a string. \
Some internal filtering going on here min and max num alleles equals two, kept ~500k sites

Testing NGSadmix on scaffold_0:

	~/./NGSadmix \
	-likes BEAGLE LIKELIHOOD FILE \
	-K NUM OF ANC POPS \

Default here is minimum 0.05 minor allele frequency

Plotting in R on personal computer with *AdmixPlot.R*

Writing bash script to calculate liklihoods for each scaffold, and concatenate them into single file for NGSadmix \
*make-beagle-file.sh*

NGS admix threw an error with this file. \
Error said a genotype likelihood was zero. Going to try filtering with VCFtools

	vcftools --vcf rosyfinches_BEAGLElikes_whole_genome.PL --out rosyfinches_BEAGLElikes_whole_genome_filtered.PL --max-missing-count 0



# 5. Fst comparisons

Beginning with pairwise comparisons of Fst across species \
Will move in to more specific comparisons later 

Calculating Weir and Cockerham's Fst implemented in vcftools \
Starting with Fst for a per-site basis, then will move into a non-overlapping window of 25kb 

	vcftools --vcf FILE --out PREFIX --weir-fst-pop POP1_SAMPLES.txt --weir-fst-pop POP2_SAMPLES.txt 

Then format the file for plotting by running *format-fst-files.sh*

Averaging across non-overlapping window of 20kb by adding this flag to the end of the above command 

	--fst-window-size 20000 \

adds suffix .windowed.weir.fst


# 6. Creating consensus sequence

Trying with bcftools consensus \
Need to bg zip file first \
Then index the file

	bcftools index INPUTFILE.vcf.gz

Then construct consensus

	bcftools consensus --sample -f REF.fa -o OUTPUT INPUT


# 7. SVDquartets

Slimming down number of variants by randomly selecting fraction of variants

	gatk SelectVariants --select-random-fraction

Then verifying the number of variants selected using
	
	gatk CountVariants

Convert the vcf file to phy using
	
	vcf2phylip.py -i FILENAME -m MIN_SAMPELS_LOCUS -o OUTGROUP -f writes FASTA -n writes NEXUS

Transfer file onto computer where SVDquartets will be run


--------------------------------------
Next Steps
--------------------------------------
snp density in vcftools?

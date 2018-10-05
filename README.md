This repository contains the shell scripts used for whole genome sequence analysis from demultiplexing to variant calls as well as a variety of other analyses related to population and phylogenomics. 

Note: Most all code here is invoked in shell scripts that have not yet been formatted to take input directly. Instead, input files and directories need to be changed in the bash script - sorry


28 May 2018 -- Currently there are two options for calling variants \
GATK HaplotypeCaller is likely best but takes a very long time \
Samtools mpileup is faster, and may be decent enough for at least preliminary results

Notes on the first iteration and a description of the pipeline on the rosy-finch data is located in rosyfinch_notes.md

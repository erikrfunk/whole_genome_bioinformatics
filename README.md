This repository contains the shell scripts used for whole genome sequence analysis from demultiplexing to variant calls as well as a variety of other analyses related to population and phylogenomics.

Notes on the first iteration and a description of the pipeline on the rosy-finch data is located in rosyfinch_notes.md

Note: Most all code here is invoked in shell scripts that have not yet been formatted to take input directly. Instead, input files and directories need to be changed in the bash script - sorry
The call of any script without arguments will show a short description and the arguments list.

The workflow is broken up into ## sections, with each breakpoint used for validation or to provide modularity.

------------------------------------------------------------------------------------------------------

Step 1: Demultiplexing \
This step can be done with any freely available demultiplexing program you prefer. The script provided here is a simple approach based on a grep search, though could be modified fairly easily to allow mismatches. Default is to make a new directory called "fastqs" in the current directory.

ARGUMENTS \
[-f] forward read file \
[-r] reverse read file \
[-b] barcode file - this file should contain the barcodes in column 1, as they appear \
in the fastq header, and the corresponding sample ID in column 2. See barcodes.txt \

OPTIONAL ARGUMENTS \
[-o] output directory - default will make new directory called 'fastqs' in current directory

Step 2: Trim and Quality Control \

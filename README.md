## README ##

# RADseq_Tools #

**RADseq_Tools** is an R pipeline for the exploratory analysis of RADseq data to guide researchers in planning their study and improving their data quality.

RADseq_Tools allows researchers to estimate important information of r the initial experiment design, including: the number of cutsites in a reference genome of a study taxa or closely related study taxa, the density of markers and distribution through the genome, and estimated genome coverage given various restriction enzymes, reference genomes, and standard sequencing machines. Even when there is not reference genome available for a particular species of interest, the user can explore with available genomes in related taxa and use this information to extrapolate the marker properties for their own experiments. In addition, this package can assist with determining the cost of a RADseq experiment by providing estimations of expected coverage, number of samples per lane, and required read throughput. 

# Getting Started
Users will need to download R. Users will also need a genome FASTA file to start the pipeline. There is a test FASTA file provided to play with, but ultimately users will need to provide their own reference genome in the FASTA format. 

## Prerequisites
There are no prerequisites required to run this package, only standard R. 

## Installing
Use the standard install.packages("RADseq_tools").

# Authors
Angel Rivera-Colon & Kira Long


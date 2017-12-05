# testing functions and utilities of the `RADTools` package
setwd("~/mercurial_repo/radseq_tools")

# source of all functions
source('./Functions.R')

# path to FASTA file
path <- './test.fa'
#path <- './test_geno.fa.gz'

# create sequence object
Seqs <- process_fasta(path, 1)
Seqs

# Sequence of the restriction site
renz <- 'CCTGCAGG'
renz <- 'ACGT'


# find all cutsite positions in the genome
all_cuts <- find_cuts(Seqs, renz)
all_cuts

# distance of cutsites
distances_all <- cutsite_distance(all_cuts)
distances_all
hist(distances_all)


cuts_n = 39000
sams = 96
Per_Sample_Coverage(cuts_n, sams)
Per_Sample_Coverage(cuts_n, sams, 'hiseq2500')
Per_Sample_Coverage(cuts_n, sams, 'hiseq4000')


Samples_Per_Lane(cuts_n, 30)
Samples_Per_Lane(cuts_n, 60)
Samples_Per_Lane(cuts_n, 60, 'hiseq4000')



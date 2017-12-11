# testing functions and utilities of the `RADTools` package
#setwd("~/mercurial_repo/radseq_tools")
setwd("~/GitHub/radseq_tools")

# source of all functions
source('./R/Functions.R')

# path to FASTA file
path <- './test.fa'
path <- './test_geno.fa.gz'
path <- './inst/extdata/test_geno.fa.gz'

# create sequence object
Seqs <- process_fasta(path, 1)
Seqs
process_fasta('./test.fa', 2)

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
mean(distances_all)
sd(distances_all)
summary(distances_all)
boxplot(distances_all)

# number of cutsites
number_cutsites(all_cuts)


cuts_n <- 39000
sams <- 96
Per_Sample_Coverage(cuts_n, sams)
is.vector(Per_Sample_Coverage(cuts_n, sams, 'hiseq2500'))
Per_Sample_Coverage(cuts_n, sams, 'hiseq4000')


Samples_Per_Lane(cuts_n, 30)
Samples_Per_Lane(cuts_n, 60)
Samples_Per_Lane(cuts_n, 60, 'hiseq4000')


DNA_Reads_Per_Lane(cuts_n, sams, 30)


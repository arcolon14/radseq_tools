\name{cutsite_distance}
\alias{cutsite_distance}

\title{
Calculate Cutsite Distance Distribution.
}
\description{
Function to calculate inter-cutsite distances. It uses R \code{diff()} to calculate the difference between per-chromosome/scaffold cutsite positions vectors. 
}

\usage{
cutsite_distance(cut_list)
}

\arguments{
  \item{cut_list}{
List containing the per-chromosome/scaffold cutsite position vectors. Argument is designed to be compatible with the output of the \code{find_cuts} function. 
  }
}

\value{
Returns a numeric atomic vector of genome-wide inter-cutsite distances. It is compatible with R vector functions.}

\author{
Angel G. Rivera-Colon
}


\examples{
# Generate prior data
myfasta <- system.file("extdata", 
                       "test_geno.fa.gz",
                       package = "RADseqTools")
mySeqs <- process_fasta(myfasta, 10000)
myEnz <- 'CTGCAG'
cutPos <- find_cuts(mySeqs, myEnz)

#
# For this function
#

# Find cutsite distance distribution in cutPos
cutDist <- cutsite_distance(cutPos)

# Object is compatible with R vector functions
# Obtain summary statistic for the cutsite distribution vector
summary(cutDist)
mean(cutDist)
sd(cutDist)
}


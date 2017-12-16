\name{number_cutsites}
\alias{number_cutsites}

\title{
Determining the number of cutsites
}
\description{
This function determines the total number of cutsites found in the genome by counting the number of cutsites from a list of cutsites made from previous function \code{find_cuts}. Output can serve as an input for functions. 
}
\usage{
number_cutsites(cutsite_list)
}

\arguments{
  \item{cutsite_list}{
List containing the per-chromosome/scaffold cutsite position vectors. Argument is designed to be compatible with the output of the \code{find_cuts} function.
}
}

\value{
Returns a single numeric value, the total number of cutsites in the provided genome. 
}

\author{
Kira Long
}

\seealso{
\code{\link{find_cuts}}
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

#Calculate number of cutsites in cutPos
nCuts <- number_cutsites(cutPos)
}

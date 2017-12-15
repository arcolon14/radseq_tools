\name{find_cuts}
\alias{find_cuts}

\title{
Find Restriction Site Positions in Genomic Sequence
}
\description{
Function to find the positions of user-defined restriction site sequences in genomic sequences of interest. Uses an R \code{gregexpr} function to find matches of a pattern (restriction site sequence) in a target character vector (genomic sequences). Positions of these matches are then saved as output. Character vectors with no matches are defined as missing values (\code{NA}).
}
\usage{
find_cuts(genomic_seq, restriction_site_seq)
}

\arguments{
  \item{genomic_seq}{
Vector containing per-chromosome/scaffold genomic sequence strings. Argument is designed to be compatible with the output of the \code{process_fasta} function. (\code{x} from \code{gregexpr} function)
}
\item{restriction_site_seq}{
String containing the restriction site sequence of interest. (\code{pattern} from \code{gregexpr} function). Sequence can be obtained from the \code{renz} dataset. 
}
}

\value{
Returns a list of vectors. Each list element contains a vector of the per-chromosome/scaffold restriction site positions. It is extracted from the first element from the \code{gregexpr} function output. 

Example:

In chromosome 1, we have cutsites at positions 7, 19, 28 and 50. No cutsites were found in chromosome 2. Cutsites were found in positions 12, 30, 36 and 42 in chromosome 3. 


\code{[[1]]}

\code{[1]  7 19 28 50}


\code{[[2]]}

\code{[1] NA}


\code{[[3]]}

\code{[1] 12 30 36 42}
}


\author{
Angel G. Rivera-Colon
}

\examples{
# Create sequences object
myfasta <- system.file("extdata",
                       "test_geno.fa.gz",
                       package = "RADseqTools")
mySeqs <- process_fasta(myfasta, 10000)

#
#For this function
#

# Define restriction enzyme sequence
data(renz)
myEnz <- renz["sbfI"]

# Find position of cutsites in sequences
cutPos <- find_cuts(mySeqs, myEnz)
}

\seealso{
\code{\link{gregexpr}}
}



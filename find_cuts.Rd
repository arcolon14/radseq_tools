\name{find_cuts}
\alias{find_cuts}

\title{
Find Restriction Site Positions in Genomic Sequence
}
\description{
Function to find the positions of user-defined restriction sites sequences in genomic sequences of interest. Uses an R \code{gregexpr} function to find matches of a pattern (restriction site sequence) in a target character vector (genomic sequences). Positions of these matches are then saved as output. Character vectors with no matches are defined as missing values (\code{NA}).
}
\usage{
find_cuts(genomic_seq, restriction_site_seq)
}

\arguments{
  \item{genomic_seq}{
Vector containing per-chromosome/scaffold genomic sequence strings. Argument is designed to be compatible with the output of the \code{process_fasta} function. (\code{x} from \code{gregexpr} function)
}
\item{restriction_site_seq}{
String containing the restriction site sequence of interest. (\code{pattern} from \code{gregexpr} function)
}
}

\value{
Returns a list of vectors. Each list elements contains a vector of the per-chromosome/scaffold restriction site positions. It is extracted from the first element from the \code{gregexpr} function output. 

Example:

\code{[[1]]}

\code{[1]  7 19 28 50}


\code{[[2]]}

\code{[1]  2 22 33 38 53}


\code{[[3]]}

\code{[1] NA}


\code{[[4]]}

\code{[1] 12 30 36 42}}
}

\author{
Angel G. Rivera-Colon
}

\examples{
# Define restriction enzyme sequence
myEnz <- 'CTGCAG'

# Find position of enzyme sequence in 
# genomic sequence object (mySeqs)
cutPositions <- find_cuts(mySeqs, myEnz)
}

\seealso{
\code{\link{gregexpr}}
}



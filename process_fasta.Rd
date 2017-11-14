\name{process_fasta}
\alias{process_fasta}

\title{
Read and Process FASTA Files. 
}
\description{
Function to read and process FASTA files into single-line sequences, R compatible objects.
}
\usage{
process_fasta(path_to_fasta)
}

\arguments{
  \item{path_to_fasta}{
Path to FASTA file. It can be zipped (\code{.gz}) or unzipped. 
}
}

\value{
Returns vector containing per-Chromosome/Scaffold sequences. 

Example: \code{"ACGTACGT...", "ACGTACGT...", "ACGTACGT...", ...}
}

\author{
Angel G. Rivera-Colon
}

\examples{
# Path to FASTA file
myfasta <- './path/to/fasta/fa'

# Create a sequences object
mysequences <- process_fasta(fasta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

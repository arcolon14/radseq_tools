\name{process_fasta}
\alias{process_fasta}

\title{
Read and Process FASTA Files. 
}
\description{
Function to read and process FASTA files into single-line sequences, R compatible objects.
}
\usage{
process_fasta(path_to_fasta, min_seq_len)
}

\arguments{
  \item{path_to_fasta}{
Path to FASTA file. It can be zipped (\code{.gz}) or unzipped. 
}
\item{min_seq_len}{
(\code{int}) Minimum sequence length. Sequences smaller than this value will not be included in the \code{sequences} object. Default: (\code{5000}). 
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
myfasta <- system.file("extdata", "test_geno.fa.gz", package = "RADseqTools")

# Create a sequences object
mySeqs <- process_fasta(myfasta, 10000)
}

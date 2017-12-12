\name{illumina}
\alias{illumina}
\docType{data}
\title{
Sequencing machine read throughput
}
\description{
Vector containing the read throughput distribution for commonly used sequencing machines. Contains expected values for low, medium and high read throughputs for the Illumina HiSeq2500 and Illumina HiSeq4000 machines. Dataset is used in as part of the \code{Samples_Per_Lane} and \code{Per_Sample_Coverage} functions.
}
\usage{data("sequencing_machines")}
\format{
  2 vectors containing 3 values each. Each value corresponds to a 'low', 'med' and 'hi' throughput.
  \describe{
    \item{\code{hiseq2500}}{\code{c(2.2e8, 3.1e8, 4.0e8)}}
    \item{\code{hiseq4000}}{\code{c(5.0e8, 7.5e8, 1.0e9)}}
  }
}

\examples{
data(illumina)

# calculate expected reads in `n` lanes
n_lanes <- 3
hiseq2500 * n_lanes

# use values as part of the Per_Sample_Coverage function
Per_Sample_Coverage(30000, 96, 'hiseq4000')
}

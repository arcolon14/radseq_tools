\name{renz}
\alias{renz}
\docType{data}
\title{
Restriction enzymes and corresponding cutsite sequences
}
\description{
A vector of commonly used restriction enzymes and their corresponding cutsite sequences. These cutsite sequences can be used in \code{find_cuts} to find the cutsites in a provided genome. This dataset can be used to explore the varying number of cutsites and cutsite density of various restriction enzymes on a genome. Selecting the correct enzyme for a RADseq study is important for the study's success. 
}
\usage{data("renz")}
\format{
  A vector containing 4 restriction enzymes and a corresponding cutesite sequence.
  \describe{
    \item{\code{pstI}}{\code{"CTGCAG"}}
    \item{\code{sbfI}}{\code{"CCTGCAGG"}}
    \item{\code{ecoRI}}{\code{"GAATTC"}}
    \item{\code{hindIII}}{\code{"AAGCTT"}}
  }
}

\examples{
data(renz)

}

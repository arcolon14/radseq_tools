\name{Per_Sample_Coverage}
\alias{Per_Sample_Coverage}

\title{
Estimate coverage per sample
}
\description{
This function will calculate the predicted per-sample coverage given the number of cutsites (from previous \code{number_cutsites()}) in the genome and a user-defined number of samples. The user can also define the Illumina platform to be used, which will determine the range of raw reads used to calculate expected coverage. This function outputs a predicted range of coverages (low, medium, high) given the set range of reads that an Illumina sequencer can give per lane. 
}
\usage{
Per_Sample_Coverage(num_cutsites, num_samples, sequencing_machine = NULL)
}

\arguments{
  \item{num_cutsites}{
Number of cutsites present in the genome. This number can be the output from \code{number_cutsites()}.
}
  \item{num_samples}{
Number of samples that the user wants to sequence in a single Illumina lane.
}
  \item{sequencing_machine}{
Illumina sequencing platform that the user wants to use for sequencing. It requires the loading of the loading of the \code{illumina} dataset. Argument has two parameters, \code{hiseq2500} and \code{hiseq4000}, each has different read throughputs. If no machine is provided, the function will default to \code{hiseq2500}.
}
}

\value{
Returns a vector with an estimate of low, medium, and high coverage.
}

\author{
Kira Long
}

\seealso{
\code{number_cutsites()}
}

\examples{
# define parameters
numberCuts <- 39000
numberSamples <- 96
data(illumina)

#Run arguments in function
Per_Sample_Coverage(numberCuts, numberSamples, "hiseq2500")

}
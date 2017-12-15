\name{DNA_Reads_Per_Lane}
\alias{DNA_Reads_Per_Lane}

\title{
Estimate per-lane reads
}
\description{
This function calculates how many reads are required from the sequencing process in order to obtain a desired coverage for a given number of samples and cutsites. This function is to assist in determining the necessary lanes of sequencing required, and/or the Illumina platform to use for a given experiment.
}
\usage{
DNA_Reads_Per_Lane(num_cutsites,
                   num_samples,
                   desired_coverage = NULL)
}

\arguments{
  \item{num_cutsites}{
Number of cutsites present in the genome. This number can be the output from \code{number_cutsites()}.
}
  \item{num_samples}{
Number of samples that the user wants to sequence in a single Illumina lane. 
}
  \item{desired_coverage}{
Number of the coverage the user wishes to keep. If no value is provided, the default value is \code{30}.
}
}

\value{
Returns a single numberic value.
}

\author{
Kira Long
}

\seealso{
\code{number_cutsites()}
}
\examples{
#Define parameters
numberCuts <- 49000
samples <- 96
coverage <- 40

#Run arguments in Function
DNA_Reads_Per_Lane(numberCuts, samples, coverage)
}

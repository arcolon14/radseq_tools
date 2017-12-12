\name{Samples_Per_Lane}
\alias{Samples_Per_Lane}

\title{
Estimate the number of samples per single Illumina lane
}
\description{
This function calculates the estimated maximum number of samples a user can multiplex in one Illumina lane given the number of cutsites, while maintaining a desired user-defined coverage.
}
\usage{
Samples_Per_Lane(num_cutsites, desired_coverage = NULL, sequencing_machine = NULL)
}

\arguments{
  \item{num_cutsites}{
Number of cutsites present in the genome. This number can be the output from \code{number_cutsites()}
}
  \item{desired_coverage}{
Number of the coverage the user wishes to keep. If no value is provided, the default value is \code{30}.
}
  \item{sequencing_machine}{
Illumina sequencing platform that the user wants to use for sequencing. It requires the loading of the loading of the \code{illumina} dataset.Argument has two parameters, \code{hiseq2500} and \code{hiseq4000}, each has different read throughputs. If no machine is provided, the function will default to \code{hiseq2500}.
}
}

\value{
Returns a vector with an estimate of low, medium, and high estimates of the number of samples a user can fit in a single Illumina lane. 
}

\author{
Kira Long
}

\seealso{
\code{number_cutsites()}
}
\examples{
#Define parameters
numberCuts <- 39000
coverage <- 30
data(illumina)

#Run arguments in Function
Samples_Per_Lane(numberCuts, coverage, "hiseq4000")

}


\name{number_cutsites}
\alias{number_cutsites}

\title{
Determining the number of cutsites
}
\description{
This function determines the total number of cutsites found in the genome by counting the number of cutsites from a list of cutsites made from previous function \code{find_cuts()}. Output can serve as an input for functions. 
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
\code{Per_Sample_Coverage()} \code{Samples_Per_Lane()} \code{DNA_Reads_Per_Lane()}
}
\examples{
#Cutsite position list (myCuts)
myCuts 

#Calculate number of cutsites in myCuts
nCuts <- number_cutsites(myCuts)
}

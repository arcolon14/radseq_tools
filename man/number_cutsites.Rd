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
An object product of \code{find_cuts()}.
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
#Define parameters/get cutsite list
cutsite_list <- list(1,2,3)

#Run argument in function
number_cutsites(cutsite_list)
}

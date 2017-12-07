\name{cutsite_distance}
\alias{cutsite_distance}

\title{
Calculate Cutsite Distance Distribution.
}
\description{
Function to calculate inter-cutsite distances. It uses R \code{diff()} to calculate the difference between per-chromosome/scaffold cutsite positions vectors. 
}

\usage{
cutsite_distance(cut_list)
}

\arguments{
  \item{cut_list}{
List containing the per-chromosome/scaffold cutsite position vectors. Argument is designed to be compatible with the output of the \code{find_cuts} function. 
  }
}

\value{
Returns a numeric atomic vector of genome-wide inter-cutsite distances. It is compatible with R vector functions.}

\author{
Angel G. Rivera-Colon
}


\examples{
# Find cutsite distance distribution in cutsite position list (cutPositions)
cutDistribution <- cutsite_distance(cutPositions)

# Object is compatible with R vector functions
# Obtain summary statistic for the cutsite distribution vector
summary(cutDistribution)

# Plot histogram of distribution
hist(cutDistribution)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

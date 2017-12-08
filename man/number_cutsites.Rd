\name{number_cutsites}
\alias{number_cutsites}

\title{
Determining the number of cutsites
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
number_cutsites(cutsite_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cutsite_list}{
%%     ~~Describe \code{cutsite_list} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cutsite_list) 
{
    count <- 0
    for (i in 1:length(cutsite_list)) {
        if (is.na(cutsite_list[[i]][1] == FALSE)) {
            count <- count + length(cutsite_list[[i]])
        }
    }
    return(count)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{resid.homogenise.lm}
\alias{resid.homogenise.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Homogenise residuals from an lm model to make them exchangeable
}
\description{
  Takes the residuals, and makes them exchangeable.
}
\usage{
resid.homogenise.lm(mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
    The model from which residuals must be homogenised.
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
  Robert Bagchi
  Maintainer: Robert Bagchi <robert.bagchi@uconn.edu>
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
function (mod) 
{
    resids <- resid(mod) * (weights(mod)^0.5)/(1 - hatvalues(mod))^0.5
    mn.resids <- mean(resids)
    resids <- resids - mn.resids
    return(resids)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

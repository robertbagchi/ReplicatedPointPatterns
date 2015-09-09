\name{kfunc.weights.calc}
\alias{kfunc.weights.calc}
\alias{abundance.weights.calc}
\alias{ratio.weights.calc}

\title{
  Calculate kfunction weights
}
\description{
  These functions calcuate the regressioin weights for the models
  following certain options.
}
\usage{
kfunc.weights.calc(pppx, pppy = NULL, r = NULL, correction = NULL, type
  = c("nx", "nx_A", "nx2", "nx2_A", "sqrtnxny", "nxny", "nxny_A",
  "sqrtnxny_A"))
abundance.weights.calc(pppx, pppy = NULL, square = FALSE, Acorr = FALSE,
  r = NULL, correction = NULL)
ratio.weights.calc(pppx, pppy = NULL, r = NULL, correction = "border")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pppx}{
    A ppp object
}
  \item{pppy}{
    A second ppp object for bivariate analyses - defaults to NULL
}
  \item{r}{
    Distances at which to calculate weights. Must start at 0.
}
  \item{correction}{
    Edge correction corresponding to options in spatstat::Kest
}
  \item{type}{
    Type of regression weights. Must be one of "nx", "nx_A", "nx2", "nx2_A", "sqrtnxny", "nxny", "nxny_A", "sqrtnxny_A"
  }
   \item{square}{
     Should the square of nx be used (or nx*ny in the bivariate
  case). Not used by user.
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
function (pppx, pppy = NULL, r = NULL, correction = NULL, type = c("nx", 
    "nx_A", "nx2", "nx2_A", "sqrtnxny", "nxny", "nxny_A", "sqrtnxny_A")) 
{
    if (is.null(r)) 
        stop("you must specify a distance range")
    switch(type, nx = abundance.weights.calc(pppx = pppx, r = r, 
        correction = correction), nx_A = abundance.weights.calc(pppx = pppx, 
        Acorr = TRUE, r = r, correction = correction), nx2 = abundance.weights.calc(pppx = pppx, 
        square = TRUE, r = r, correction = correction), nx2_A = ratio.weights.calc(pppx = pppx, 
        pppy = NULL, r = r, correction = correction), sqrtnxny = sqrt(abundance.weights.calc(pppx = pppx, 
        pppy = pppy, square = TRUE, r = r, correction = correction)), 
        nxny = abundance.weights.calc(pppx = pppx, pppy = pppy, 
            square = TRUE, r = r, correction = correction), nxny_A = ratio.weights.calc(pppx = pppx, 
            pppy = pppy, r = r, correction = correction), sqrtnxny_A = sqrt(ratio.weights.calc(pppx = pppx, 
            pppy = pppy, r = r, correction = correction)), stop(paste(type, 
            "is not a valid selection.", "\nChoose from 'nx', 'nx_A','nx2', 'nx2_A','sqrtnxny', 'nxny', 'nxny_A', 'sqrtnxny_A'")))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

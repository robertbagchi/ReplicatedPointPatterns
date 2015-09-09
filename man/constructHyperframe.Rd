\name{constructHyperframe}
\alias{constructHyperframe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculates the K functions and weights for a hyperframe
}
\description{
  Takes a hyperframe and calcualtes the K functions and the weights for
  that hyperframe given the edge correction and weights type.
}
\usage{
constructHyperframe(hyper, r, correction, pppx='pppx', weights.type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hyper}{
    A hyperframe that must have an element called pppx 
}
  \item{r}{
    distances at which to compute r
}
  \item{correction}{
    The edge correction by name - options are those for Kest in
    spatstat. Defaults to 'border'
  }
  \item{pppx}{
    The name of the element with the ppp object. Experimental.
    }
\item{weights.type}{
   The type of regression weights to be calculated. Options are 'nx',
    'nx_A', 'nx2', 'nx2_A', 'sqrtnxny', 'nxny', 'nxny_A', 'sqrtnxny_A'
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
function (hyper, r, correction, weights.type) 
{
    if (min(r) > 0) 
        r <- c(0, r)
    if (!all(c("pppx") \%in\% names(hyper))) {
        stop("hyperframe object must include 'pppx' element")
    }
    hyper$K <- with.hyperframe(hyper, Kest(pppx, r = r, correction = correction, 
        ratio = TRUE))
    hyper$wts <- with.hyperframe(hyper, list(kfunc.weights.calc(pppx, 
        r = K$r, correction = correction, type = weights.type)))
    minsamp <- sapply(with.hyperframe(hyper, list(kfunc.weights.calc(pppx, 
        r = K$r, correction = correction, type = "nx"))), function(x) min(x[[1]]))
    hyper$minsamp <- minsamp
    return(hyper)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

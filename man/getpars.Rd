\name{getpars}
\alias{getpars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extract parameters from Kfunc lme models
}
\description{
  Used internally to extract the parameters from a kfunclme object.
}
\usage{
getpars(mod, lin.comb.Ct)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
    Model to get parameters from
}
  \item{lin.comb.Ct}{
    Linear combination to be used for predictions
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
function (mod, lin.comb.Ct) 
{
    if (!is.null(mod)) {
        beta.r <- fixef(mod)
        vcov.r <- vcov(mod)
        est.Kmean.r <- as.vector(lin.comb.Ct \%*\% fixef(mod))
        est.Kse.r <- sqrt(diag(lin.comb.Ct \%*\% vcov(mod) \%*\% 
            t(lin.comb.Ct)))
    }
    else {
        beta.r <- rep(NA, ncol(lin.comb.Ct))
        vcov.r <- matrix(NA, ncol = length(beta.r), nrow = length(beta.r))
        est.Kmean.r <- rep(NA, nrow(lin.comb.Ct))
        est.Kse.r <- rep(NA, nrow(lin.comb.Ct))
    }
    return(list(pred.r = est.Kmean.r, se.pred.r = est.Kse.r, 
        beta.r = beta.r, vcov.r = vcov.r))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

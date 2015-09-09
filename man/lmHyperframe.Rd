\name{lmHyperframe}
\alias{lmHyperframe}
\alias{kfunclm}

\alias{print.kfunctionlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  linear model of Kfunctions
}
\description{
  Function to take a hyperframe and fit a linear model to the data
  using all rows and the distance range
  Still very much in development so use carefully
  Not yet generalised to bivariate designs
  Not tested much yet
}
\usage{
lmHyperframe(hyper, r, form, correction = "border", weights.type=NULL,
minsamp = NA, computeK=TRUE,  printwarnings = TRUE)
kfunclm(k, dat, form, weights)

\method{print}{kfunctionlm}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hyper}{
    An object of class hyperframe (spatstat) containing point patterns
    and covariates
  }
  \item{r}{
    distances at which to model - must include only increasing positive
    numbers starting from 0
  }
  \item{form}{
    The names of the explanatory covariates in quotes
  }
  \item{correction}{
    The edge correction by name - options are those for Kest in
    spatstat. Defaults to 'border'
  }
  \item{weights.type}{
    The type of regression weights to be calculated. Options are 'nx',
    'nx_A', 'nx2', 'nx2_A', 'sqrtnxny', 'nxny', 'nxny_A', 'sqrtnxny_A'
    
  }
  \item{minsamp}{
    Defines the minimum number of points that need to be in the central
    part of the plot (i.e. not affected by edge corrections). Defaults
    to NA so no cases are removed.    
  }
  \item{computeK}{
    Controls whether hyperframe should compute K or if this is done
    externally.
  }
  \item{printwarnings}{
    Whether or not to print warnings
  }
  \item{k}{
    Vector of K(r) for all cases for a given distance r.
  }
  \item{dat}{
   Data set with covariates included in formula.
 }
 \item{weights}{
   Vector of weights for all cases for a given distance.
 }
 \item{x}{
   Object to be printed
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
function (hyper, r, form, correction = "border", weights.type, 
    minsamp = 1, printwarnings = TRUE) 
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
    sp.keep <- sapply(hyper$wts, function(x) all(unlist(x[r]) >= 
        minsamp))
    removed.species <- row.names(hyper)[!sp.keep]
    warning(paste("Removed", length(removed.species), "species with insufficient numbers"))
    hyper <- hyper[sp.keep, ]
    dist.keep <- (apply(sapply(hyper$K, function(K) K[[correction]]), 
        1, function(x) var(x)) > 0)
    warning(paste("Not modelling K at distances ", paste(r[!dist.keep], 
        collapse = ", "), "due to zero variance"))
    modr <- match(r[dist.keep], r)
    kmods <- sapply(modr, function(i) {
        kfunclm(k = sapply(hyper$K, function(k) k[[correction]][i]), 
            dat = as.data.frame(hyper, warn = FALSE), form = form, 
            weights = sapply(hyper$wts, function(w) unlist(w)[i]))
    }, simplify = FALSE)
    names(kmods) <- r[dist.keep]
    class(kmods) <- "kfunctionlm"
    return(kmods)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

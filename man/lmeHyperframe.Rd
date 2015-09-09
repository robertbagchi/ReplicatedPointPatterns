\name{lmeHyperframe}
\alias{lmeHyperframe}
\alias{kfuncLme}
\alias{refit.lmek}
\alias{print.kfunctionlme}

\title{
  Fits a set of mixed effects models to a hyperframe object
}
\description{
  Fits a mixed effects model to a given hyperframe object with  K
  function (K) and weights (wts) elements.

  kfunclme and refit.lmek are internally called by lmeHyperframe.
}
\usage{
lmeHyperframe(hyper, r, fixed, random, correlation = NULL, correction =
"border", weights.type, computeK=TRUE, minsamp = NA, printwarnings =
TRUE)

kfuncLme(k, dat, weights, fixed, random, correlation, na.action = na.omit)
refit.lmek(mod, res.r)

\method{print}{kfunctionlme}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hyper}{
    Hyperframe object created with spatstat(as.hyperframe). Must K and
    wts elements.
}
  \item{r}{
    Distances at which to model K. Should include 0
}
  \item{fixed}{
    RHS of the fixed effects formula in quotes, without a tilda
    (e.g. "x1 + x2")
}
  \item{random}{
    RHS of the random effects formula in quotes (e.g. "1|group")
}
  \item{correlation}{
    The correlation argument for models - sent to nlme::lme and must be
    accepted by nlme::lme. 
}
  \item{correction}{
    The edge correction which must correspond to an option in spatstat::Kest
}
  \item{weights.type}{
    The regression weights. Must be one of "nx", "nx_A", "nx2", "nx2_A", "sqrtnxny", "nxny", "nxny_A", "sqrtnxny_A"
  }
  \item{computeK}{
    Should K be computed internally. Defaults to TRUE
    }
  \item{minsamp}{
    The minimum number of points, after edge corrections, for a
    case to be included. Cases with fewer points will be excluded and
    row names returned as an attribute.
  }
  \item{printwarnings}{
    Should warnings about case or distance exclusion be
    reported. Defaults to TRUE
  }
   \item{k}{
     A vector corresponding the K(r) for all cases at a given distance r
}
  \item{dat}{
    Data including the covariates used in the fixed and random parts of
    the model.
  }
   \item{weights}{
     The model weights passed to lme.
   }
   \item{na.action}{
     How to deal with NAs in the data set. Defaults to na.omit.
   }
   \item{mod}{
     An lme model at each distance. Used internally by refit.lmek
   }
   \item{res.r}{
     Vector of residuals. Used internally by refit.lmek
   }
   \item{x}{
     Object to be printed.
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
function (hyper, r, fixed, random, correlation = NULL, correction = "border", 
    weights.type, minsamp = 10, printwarnings = TRUE) 
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
    sp.keep <- sapply(with.hyperframe(hyper, list(kfunc.weights.calc(pppx, 
        r = K$r, correction = correction, type = "nx"))), function(x) all(unlist(x[r]) >= 
        minsamp))
    removed.species <- row.names(hyper)[!sp.keep]
    if (printwarnings) 
        warning(paste("Removed", length(removed.species), "species with insufficient numbers"))
    hyper <- hyper[sp.keep, ]
    dist.keep <- (apply(sapply(hyper$K, function(K) K[[correction]]), 
        1, function(x) var(x)) > 0)
    if (printwarnings) 
        warning(paste("Not modelling K at distances ", paste(r[!dist.keep], 
            collapse = ", "), "due to zero variance"))
    modr <- match(r[dist.keep], r)
    kmods <- sapply(modr, function(i) {
        kfuncLme(k = sapply(hyper$K, function(k) k[[correction]][i]), 
            dat = as.data.frame(hyper, warn = FALSE), fixed = fixed, 
            random = random, correlation = correlation, weights = sapply(hyper$wts, 
                function(w) unlist(w)[i]))
    }, simplify = FALSE)
    names(kmods) <- r[dist.keep]
    return(kmods)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

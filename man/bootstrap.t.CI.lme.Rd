\name{bootstrap.t.CI.lme}
\alias{bootstrap.t.CI.lme}
\alias{bootstrap.compare.lme}
\alias{bootstrap.parallel.lme}
\alias{compare.mods.bootstrap}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bootstrapping for CIs and hypothesis test on kfunclme objects
}
\description{
  Carries out bootstrapping on kfunclme objects to get confidence
  intervals on parameter estimates and predictions.
  
  bootstrap.compare.lme compares two nested model fits across a range of
  distances to determine the loss in explanatory power of a covariate
  and tests the null hypothesis of no effect of the covariate.
}
\usage{
bootstrap.t.CI.lme(mods, lin.comb.Ct, nboot, alpha, ncore = 1, transform
= NULL)

bootstrap.compare.lme(mods, term, dists, nboot, ncore)
bootstrap.parallel.lme(mods, resids, lin.comb.Ct, nboot, ncore = 1)
compare.mods.bootstrap(modH0, modH1, res.r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mods}{
    The models for each distance r.
}
  \item{lin.comb.Ct}{
    The model matrix to be multiplied by the fixed effects to get predictions.
}
  \item{nboot}{
    Number of iterations
}
  \item{alpha}{
    Confidence level.
}
  \item{ncore}{
    Number of cores to use in parallel computations
}
  \item{transform}{
    Experimental and not fully implemented - function to transform the response
  }
  \item{term}{
    Term in model formula to be dropped in simpler, nested, model.
  }
  \item{dists}{
    Distances to be considered in model test.
  }
   \item{resids}{
     List of residuals to be used in bootstrapping. Not usually specified by user.
   }
   \item{modH0}{
     Simple model
   }
   \item{modH1}{
     Complex model
   }
   \item{res.r}{
     Residuals after randomisation.
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
function (mods, lin.comb.Ct, nboot, alpha, ncore = 1, transform = NULL) 
{
    require("abind")
    require("parallel")
    resids <- lapply(mods, residual.homogenise.lme)
    boot.esti <- bootstrap.parallel.lme(mods = mods, resids = resids, 
        lin.comb.Ct = lin.comb.Ct, nboot = nboot, ncore = ncore)
    sample.esti <- sapply(mods, getpars, lin.comb.Ct = lin.comb.Ct, 
        simplify = FALSE)
    boot.pars <- sapply(boot.esti, function(r, est) {
        mapply(function(x, est) {
            t.sim <- (x$beta.r - est$beta.r)/sqrt(diag(x$vcov.r))
            est$beta.r - t.sim * sqrt(diag(est$vcov.r))
        }, x = r, est = est, SIMPLIFY = TRUE)
    }, est = sample.esti, simplify = FALSE)
    boot.fix.cis <- apply(do.call("abind", args = list(what = boot.pars, 
        along = 3)), c(2, 1), quantile, c(alpha/2, 1 - alpha/2))
    sample.fix.cis <- aperm(sapply(sample.esti, function(x) x$beta.r), 
        c(2, 1))
    sample.fix.cis <- array(sample.fix.cis, dim = c(1, dim(sample.fix.cis)))
    modelpars <- abind(list(estimate = sample.fix.cis, boot.fix.cis), 
        along = 1)
    t.score <- lapply(boot.esti, function(bootsamp, obssamp) {
        mapply(function(sim, obs) {
            t.r <- (sim$pred.r - obs$pred.r)/sim$se.pred.r
            return(t.r)
        }, bootsamp, obssamp)
    }, obssamp = sample.esti)
    t.score <- do.call("abind", args = list(what = t.score, along = 3))
    uci <- apply(t.score, c(2, 1), quantile, alpha/2, na.rm = T)
    lci <- apply(t.score, c(2, 1), quantile, 1 - alpha/2, na.rm = T)
    uci <- split(uci, row(uci))
    lci <- split(lci, row(lci))
    CIs <- mapply(function(obs, ucl, lcl) {
        lower.CI <- obs$pred.r - lcl * obs$se.pred.r
        upper.CI <- obs$pred.r - ucl * obs$se.pred.r
        return(list(LCL = lower.CI, UCL = upper.CI))
    }, obs = sample.esti, ucl = uci, lcl = lci, SIMPLIFY = FALSE)
    estimator <- lapply(sample.esti, function(x) x$pred.r)
    lower.CI <- lapply(CIs, function(x) x$LCL)
    upper.CI <- lapply(CIs, function(x) x$UCL)
    bootstrap.CI <- list(estimator = estimator, lower = lower.CI, 
        upper = upper.CI, modelpars = modelpars)
    return(bootstrap.CI)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

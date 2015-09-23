\name{lm.t.boot}
\alias{lm.t.boot}
\alias{lm.boot}
\title{
  calculate t-distribution based confidence intervals for an kfunclm object
}
\description{
  Carries out bootstrapping on a kfunclm object to get confidence intervals.
}
\usage{
lm.t.boot(lmmods, lincomb, nsim, alpha, simple.method = TRUE)
lm.boot(mods, lincomb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmmods}{
    List of linear models fitted to each distance.
}
  \item{lincomb}{
    The model matrix used to multiply the coefficients to get predictions.
}
  \item{nsim}{
    Number of bootstraps.
}
  \item{alpha}{
    Confidence level.
}
  \item{simple.method}{
    Whether to use a simple method based on just the quantiles of the
  predictions or to use an empirical t-distribution. Defaults to TRUE,
  but maybe changed in the future.
}
  \item{mods}{
    The individual models used in the bootstrap.
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
function (lmmods, lincomb, nsim, alpha, simple.method = TRUE) 
{
    require(abind)
    bsci <- replicate(nsim, lm.boot(lmmods, lincomb = lincomb), 
        simplify = FALSE)
    obs.pred <- lapply(lmmods, function(mod) {
        pred <- lincomb \%*\% coef(mod)
        se.pred <- sqrt(diag(lincomb \%*\% vcov(mod) \%*\% t(lincomb)))
        return(list(pred = pred, se.pred = se.pred))
    })
    if (simple.method) {
        cis <- abind(lapply(bsci, function(iter) {
            do.call("cbind", lapply(iter, function(kd) kd$pred))
        }), along = 3)
        cis <- apply(cis, c(2, 1), quantile, c(alpha/2, 1 - alpha/2))
        lower <- cis[1, , ]
        if (class(lower) != "matrix") 
            lower <- as.matrix(lower)
        lower.CI <- split(lower, row(lower))
        upper <- cis[2, , ]
        if (class(upper) != "matrix") 
            upper <- as.matrix(upper)
        upper.CI <- split(upper, row(upper))
    }
    else {
        t.scores <- lapply(bsci, function(sim, obs) {
            do.call("cbind", mapply(function(obs.t, sim.t) {
                t.score.t <- (sim.t$pred - obs.t$pred)/sim.t$se.pred
                t.score.t <- matrix(t.score.t, nc = 1)
                return(t.score.t)
            }, obs.t = obs, sim.t = sim, SIMPLIFY = F))
        }, obs = obs.pred)
        t.scores <- do.call("abind", args = list(t.scores, along = 3))
        uci <- apply(t.scores, c(2, 1), quantile, alpha/2, na.rm = T)
        lci <- apply(t.scores, c(2, 1), quantile, 1 - alpha/2, 
            na.rm = T)
        uci <- split(uci, row(uci))
        lci <- split(lci, row(lci))
        CIs <- mapply(function(obs, ucl, lcl) {
            lower.CI <- obs$pred - lcl * obs$se.pred
            upper.CI <- obs$pred - ucl * obs$se.pred
            return(list(LCL = lower.CI, UCL = upper.CI))
        }, obs = obs.pred, ucl = uci, lcl = lci, SIMPLIFY = FALSE)
        lower.CI <- lapply(CIs, function(x) x$LCL)
        upper.CI <- lapply(CIs, function(x) x$UCL)
    }
    estimator <- lapply(obs.pred, function(x) x$pred)
    return(list(lmK = lmmods, lmKpred = estimator, lower = lower.CI, 
        upper = upper.CI))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

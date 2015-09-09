\name{residual.randomise.lme}
\alias{residual.randomise.lme}
\alias{residual.homogenise.lme}

\title{
  Randomise and homogenise residuals from an lme model and add back heterogeneity
}
\description{
  Functions used to randomise and homogenise residuals of an lme model.
}
\usage{
residual.randomise.lme(mods, resids)
residual.homogenise.lme(mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mods}{
    List of models
  }
  \item{resids}{
    Residuals from model
  }
  \item{mod}{
    Individual model for homogenising the residuals.
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
function (mods, resids) 
{
    level1.resid <- lapply(resids, function(x) x[["level1resids"]])
    N <- max(sapply(level1.resid, length))
    indx <- sample(1:N, replace = T)
    level1.resid.r <- lapply(level1.resid, function(x) return(x[indx]))
    Cmat <- sapply(resids, function(x) attr(x, "zmat"), simplify = FALSE)
    level1.resid.raw.r <- mapply(function(mod, res, Clarge) {
        if (is.null(mod)) 
            return(NULL)
        else {
            wts <- getCovariate(mod$modelStruct$varStruct)[order(order(getGroups(mod)))]
            transform1 <- t(chol(diag(sqrt(wts)) \%*\% Clarge \%*\% 
                diag(sqrt(wts))))
            level1.res.raw.r <- as.vector(transform1 \%*\% res)
            return(level1.res.raw.r)
        }
    }, mod = mods, res = level1.resid.r, Clarge = Cmat, SIMPLIFY = F)
    samp <- lapply(mods, function(mod) {
        if (is.null(mod)) 
            return(NULL)
        else {
            re <- ranef(mod)
            return(lapply(re, function(rj) {
                sample(1:NROW(rj), replace = T)
            }))
        }
    })
    samp <- samp[[1]]
    ranef.r <- mapply(function(mod, res, samp) {
        if (is.null(mod)) 
            return(NULL)
        else {
            ranef.res <- res[-which(names(res) == "level1resids")]
            ranef.res.r <- mapply(function(r, ord) {
                rnew <- as.matrix(r[ord, ])
                rownames(rnew) <- rownames(r)
                return(rnew)
            }, r = ranef.res, ord = samp, SIMPLIFY = FALSE)
            ranef.res.new <- mapply(function(j, r, mod) {
                g <- as.character(getGroups(mod, level = j))
                if (is.null(attr(r, "rownames"))) 
                  rownames(r) <- rownames(ranef(mod, level = j))
                return(r[g, ])
            }, j = as.list(1:length(ranef.res)), r = ranef.res.r, 
                MoreArgs = list(mod = mod), SIMPLIFY = FALSE)
        }
    }, mod = mods, res = resids, MoreArgs = list(samp = samp), 
        SIMPLIFY = FALSE)
    resids <- mapply(function(level1, ranef) {
        list(level1.resid.raw.r = level1, ranef.r = ranef)
    }, level1 = level1.resid.raw.r, ranef = ranef.r, SIMPLIFY = F)
    return(resids)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

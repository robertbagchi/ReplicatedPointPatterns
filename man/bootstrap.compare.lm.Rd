\name{bootstrap.compare.lm}
\alias{bootstrap.compare.lm}

\title{
  Bootstrapping tests of individual terms in the model.
}

\description{
  Uses bootstrapping to test the null hypothesis that a given term in the
  model explains no variance in the second order spatial
  structure.
}

\usage{
bootstrap.compare.lm(mods, term, dists=NULL, nboot)
}

\arguments{
  \item{mods}{
    A list with models for each distance r.
  }
  \item{term}{
    The term to be tested
  }
  \item{dists}{
    The distance range to be used in the tests.
  }
  \item{nboot}{
    Number of iterations in the bootstrap
  }
  
  \details{
  }
  \author{
     Robert Bagchi
  Maintainer: Robert Bagchi <robert.bagchi@uconn.edu>

}
\note{
}

\seealso{
}

\examples{
bootstrap.compare.lm <- function(mods, term, dists=NULL, nboot){
  
  if(is.null(dists))
    dists <- 1:length(mods)
  
  if(length(mods) < length(dists))
    stop('More test distances than modelled distances')
  
  mods <- mods[dists]
  modsH1 <- mods
  
  ## null models
  modsH0 <-  lapply(modsH1, function(mod, term){
    if(!is.null(mod)){
      mod0 <- update(mod, paste('~.-', term))
    }
    ##  if(class(mod)=='try-error'){
    ##     warning("Null model fit did not converge at some distances")
    ##     mod <- NULL
    ##   }
    ## }
    else mod0 <- NULL
    return(mod0)}, term=term)
  
  ## calculate the test statistic
  obs.stat <- mapply(function(mod, mod0){
    if(is.null(mod)|is.null(mod0))
      return(NULL)
    else
      (-2)*logLik(mod0) - (-2)*logLik(mod)}, mod=modsH1, mod0=modsH0)
  
  ## now use a bootstrap to develop the null distribution for this statistic
  boot.stat <- sapply(1:ceiling(nboot*1.2),  function(i, modsH0, modsH1){
    ## First extract the residuals from the model
    resids.H1 <- lapply(modsH1, resid.homogenise.lm)
    samp <- sample(1:max(sapply(resids.H1, length)), replace=T) ## set up sample to be
    # repeated at all distances
    boot.stat <- mapply(function(mod1, mod0, resids, indx){
      k.data.frame <- as.data.frame(mod1$model)
      mod0 <- eval(getCall(mod0))
      mod1 <- eval(getCall(mod1))
      form0 <- formula(mod0)
      form1 <- formula(mod1)
      if(length(resids[indx])==length(mod0$weights)){
          newK <- fitted(mod0) +   resids[indx]/(weights(mod0)^0.5)
        k.data.frame$Kr <- newK
        k.data.frame$wts <- mod0$weights
        if(!any(is.na(newK))){
          modnew0 <- update(lm(form0, weights=wts, data=k.data.frame), Kr~.)
          modnew1 <- update(lm(form1, weights=wts, data=k.data.frame), Kr~.)
          boot.stat <- (-2)*logLik(modnew0) - (-2)*logLik(modnew1)
          return(boot.stat)
        }
        else(return(NA))
      }
      else(return(NA))
      
    }, mod1=modsH1, mod0=modsH0, resids=resids.H1, MoreArgs=list(indx=samp))
  }, modsH0=modsH0, modsH1=modsH1, simplify=TRUE)
  
  boot.stat <- boot.stat[,apply(boot.stat, 2, function(x) all(!is.na(x)))]
  boot.stat
  if(ncol(boot.stat) < nboot)
    warning(paste("only", length(boot.stat), 'completed simulations - increase iterations?'))
  else
    boot.stat <- boot.stat[,1:nboot]
  D.obs <- sum(obs.stat[dists])
  
  D.boot <- apply(boot.stat, 2, function(x, d){ return(sum(x[d]))}, d=dists)
  
  p.val <- (sum(D.obs < D.boot)+1)/(1+length(D.boot))
  p.val
  result <- list(D=D.obs, p=p.val, D.boot=D.boot)
  return(result)}

}

\keyword{ ~kwd1 }

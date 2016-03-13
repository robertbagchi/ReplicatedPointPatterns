########################################################################
## File name: RepSPPFunctions_new.R
## Author: Robert Bagchi (bagchi.r@gmail.com)
## Last modified 18/11/2014 by R Bagchi
## Description Functions for replicated spatial point pattern analysis
##             including both fixed and random effects and functions for
##             predictions and bootstrapped confidence intervals.
##             Supports the paper Bagchi & Illian,
##             "A method for replicated spatial point  pattern analysis
##             in ecology" submitted to MEE
##             Code has been tested in a limited number of situations
##             so use with care - please report any bugs to Robert Bagchi
########################################################################
########################################################################
## ratio covariate - takes a list of point patterns
## Calculates the denominator of the K function
## and returns a matrix of values 
########################################################################
ratio.weights.calc <- function(pppx, pppy=NULL, r=NULL, correction='border'){
  if(is.null(correction))
    stop('you must define a correction for this weights argument')
  if(is.null(pppy))
    Kx <- Kest(pppx, r=r, correction=correction, ratio=TRUE)
  else
    Kx <- Kcross(superimpose(x=pppx, y=pppy, W=Window(pppy)),
                 r=r, correction=correction, ratio=TRUE)
  wts <-  attr(Kx, 'denominator')[[correction]]
  return(wts)}
####################################################
## A function to calculate the weights based on the abundance
## of points only
########################################################
abundance.weights.calc <- function(pppx, pppy=NULL,
                                   square=FALSE, Acorr = FALSE,
                                   r=NULL, correction=NULL){
  ## first need to restrict calculation of n to
  ## points < r away if border correction  
  if(!is.null(correction)) ## if correction - implement the correction
    if(correction=='border')
      pppx.c <- sapply(r, function(ri) pppx[erosion(pppx$window, ri)],
                       simplify=FALSE)
    else
      stop(paste("Correction for", correction,
                 "method not implemented yet")) 
  
  else
    pppx.c <-  sapply(r, function(r) pppx, simplify=FALSE)
  
  wx <- sapply(pppx.c, npoints) ## number of points

  if(square)
    if(!is.null(pppy))
      wx <- wx*npoints(pppy)
  else
    wx <- wx*npoints(pppx)
  if(Acorr)
    if(!is.null(pppy))
      wx <- wx/area.owin(pppy)
  else
    wx <- wx/area.owin(pppx)
  
  return(wx)
}
#############################################################################
## a function for constructing the weights argument from a few simple rules.
## Uses thje functions abundance.weights.calc and ratio.weights.calc
## to do the actual calculations
##########################################################################
kfunc.weights.calc <- function(pppx, pppy=NULL, r=NULL, correction=NULL,
                               type=c('nx', 'nx_A', 'nx2', 'nx2_A',
                                 'sqrtnxny', 'nxny', 'nxny_A', 'sqrtnxny_A')) {
  if(is.null(r))
    stop('you must specify a distance range')
  
  switch(type,
         nx = abundance.weights.calc(pppx=pppx, r=r, correction=correction),
         nx_A = abundance.weights.calc(pppx=pppx, Acorr=TRUE,
           r=r, correction=correction),
         nx2 = abundance.weights.calc(pppx=pppx, square=TRUE,
           r=r, correction=correction),
         nx2_A = ratio.weights.calc(pppx=pppx, pppy=NULL,
           r=r, correction=correction),
         sqrtnxny = sqrt(abundance.weights.calc(pppx=pppx, pppy=pppy,
           square=TRUE,  r=r, correction=correction)),
         nxny = abundance.weights.calc(pppx=pppx, pppy=pppy, square=TRUE,
           r=r, correction=correction),
         nxny_A=ratio.weights.calc(pppx=pppx, pppy=pppy,
           r=r, correction=correction),
         sqrtnxny_A=sqrt(ratio.weights.calc(pppx=pppx, pppy=pppy,
           r=r, correction=correction)),
         stop(paste(type, 'is not a valid selection.',
                    "\nChoose from 'nx', 'nx_A','nx2', 'nx2_A','sqrtnxny', 'nxny', 'nxny_A', 'sqrtnxny_A'"))
         )

}

########################################################################
## Function to do a linear regression of K functions
########################################################################
kfunclm <- function(k, dat, form, weights){
  k.data.frame <- data.frame(K=k, dat, wts=weights) ## extract data
  fixed.formula <- as.formula(paste('K~', form))    ## specify the formula
  ## fit model 
  mod <- eval(substitute(lm(fixed.formula, data=k.data.frame, weights=wts,
                              x=T, y=T), list(fixed.formula=fixed.formula,
                                     k.data.frame=k.data.frame)))
  return(mod)}

########################################################################
## function to homogonise the residuals of a linear model
## to make them exchangable.
########################################################################

resid.homogenise.lm <- function(mod){
  ## extract resids, multiply by sqrt of weights
  ## and divide by hat values - taken from p279 in
  # Davison & Hinkley, 1997
  resids <- resid(mod)*(weights(mod)^0.5)/(1-hatvalues(mod))^0.5
  mn.resids <- mean(resids) # centre the residuals
  resids <- resids-mn.resids
  return(resids)}

################################################################################
## Function to take a hyperframe and fit a linear model to the data
## using all rows and the distance range
## Still very much in development so use carefully
## Not yet generalised to bivariate designs
## Not tested much yet
###############################################################################
constructHyperframe <- function(hyper, r, correction, pppx='pppx', weights.type)
  {
    if(min(r) > 0)
      r <- c(0, r)
    
    if(!(pppx %in% names(hyper)))
      {
        stop("hyperframe object must include 'pppx' element")
      }

    names(hyper)[names(hyper)==pppx] <- 'pppx'
    
    ## calculate the K-functions
    hyper$K <- with.hyperframe(hyper, Kest(pppx, r=r, correction=correction,
                                           ratio=TRUE))
    ## note that we need the "list" function because otherwise hyperframe
    ## throws a fit.
    hyper$wts <- with.hyperframe(hyper,
                                 list(kfunc.weights.calc(pppx, r=K$r,
                                                         correction=correction,
                                                         type=weights.type)))

    minsamp <- sapply(
                     with.hyperframe(hyper,
                                     list(kfunc.weights.calc(pppx, r=K$r,
                                                             correction=correction,
                                                             type='nx'))),
                      function(x) min(x[[1]]))
    hyper$minsamp <- minsamp
 
    return(hyper)
  }



lmHyperframe <- function(hyper, r , form,  correction='border',
                         weights.type=NULL,  minsamp=NA,
                         computeK =TRUE, printwarnings=TRUE){

  if(!all(c('K', 'wts') %in% names(hyper)))
    if(!computeK)
      stop("hyperframe object must include 'K' and 'wts'")
  else
    if(!('pppx' %in% names(hyper)))
      stop("hyperframe object must include 'pppx' entry if 'computeK' is TRUE")
    else
      hyper <- constructHyperframe(hyper=hyper, r=r, correction=correction,
                          pppx='pppx', weights.type=weights.type)
        
  ## find species that have sufficient individuals (define with minsamp)
  ## across the plot for analysis
  if(!is.na(minsamp))
    {
      sp.keep <-sapply(
                       with.hyperframe(hyper,
                                       list(kfunc.weights.calc(pppx, r=K$r,
                                                               correction=correction,
                                                               type='nx'))),
                       function(x) all(unlist(x[r]) >= minsamp))

      removed.species <- row.names(hyper)[!sp.keep]
      ## select these species
      hyper <- hyper[sp.keep,]
      if(printwarnings & length(removed.species) > 0)
        warning(paste("Removed", length(removed.species),
                      "species with insufficient numbers"))
    }
  
  ## Do not model distances where the variance is 0
  dist.keep <-  (apply(sapply(hyper$K, function(K) K[[correction]]), 1,
                       function(x) var(x)) > 0)
  
  warning(paste('Not modelling K at distances ',
                paste(r[!dist.keep], collapse=', '),
                "due to zero variance"))
  modr <-  match(r[dist.keep], r )

  ## fit the models
  kmods <- sapply(modr, function(i)
                  {
                    kfunclm(k=sapply(hyper$K, function(k) k[[correction]][i]),
                            dat=as.data.frame(hyper, warn=FALSE), form=form,
                            weights=sapply(hyper$wts, function(w) unlist(w)[i]))
                  }, simplify=FALSE)
  names(kmods) <- r[dist.keep]
  attr(kmods, 'removed.species') <-  removed.species
  class(kmods) <- 'kfunctionlm'
  return(kmods)}

## method to print kfunction lm objects
print.kfunctionlm <- function(x, ...){
  dists <- as.numeric(names(x))
  cat('linear model fitted to k function \n')
  cat("distances modelled:\n", length(dists), 'distances\n',
      'range =', range(dists), '\n')
  cat('coefficients\n')
  print(rbind(dists, sapply(x, coef)))
}

########################################################################
## bootstrapping function for linear models
########################################################################
lm.boot <- function(mods, lincomb){
  resids <- lapply(mods, resid.homogenise.lm) ## extract exchangable residuals
  samp <- sample(1:max(sapply(resids, length)), replace=T) ## set up sample to be
                                        # repeated at all distances
  ## extract the bootstraped K function for 1 iteration
  preds.all <- mapply(function(mod, resids, lincomb,  indx){
    k.data.frame <- as.data.frame(mod$model)
    mod <- eval(getCall(mod))
    if(length(resids[indx])==length(mod$weights)){
     
      newK <- fitted(mod) +   resids[indx]/(weights(mod)^0.5)

      k.data.frame$Kr <- newK
      if(!any(is.na(newK))){
        modnew <- update(mod, Kr~., weights=mod$weights, data=k.data.frame)
        pred <- lincomb %*% coef(modnew)
        se.pred  <- sqrt(diag(lincomb %*% vcov(mod) %*% t(lincomb)))
      }
    }

    else {
      pred <- rep(NA, length(resids))
      se.pred  <- rep(NA, length(resids))
    }
    return(list(pred=pred, se.pred=se.pred))},
                  mod=mods, resids=resids, MoreArgs=list(indx=samp, lincomb=lincomb),
                      SIMPLIFY=FALSE)
  return(preds.all)
}

lm.t.boot <- function(lmmods, lincomb, nsim, alpha, simple.method=TRUE){
  ###require(abind)
  ## repeats the bootstrap from lm.boot nsim times
  bsci <- replicate(nsim, lm.boot(lmmods, lincomb=lincomb), simplify=FALSE)

  ## extracts the observed k-function and standard error
  obs.pred <- lapply(lmmods, function(mod){
    pred <- lincomb %*% coef(mod)
    se.pred <- sqrt(diag(lincomb %*% vcov(mod) %*% t(lincomb)))
    return(list(pred=pred, se.pred=se.pred))})


  ## simple method to calculate the confidence intervals
  ## as the quantile of the simulations
   if(simple.method){
      cis <- abind(lapply(bsci, function(iter){
          do.call('cbind',lapply(iter, function(kd) kd$pred))}), along=3)
      cis <- apply(cis, c(2, 1), quantile, c(alpha/2, 1-alpha/2))
      lower <-  cis[1,,]
      if(class(lower) != 'matrix') lower <- as.matrix(lower)
      lower.CI <- split(lower, row(lower))
      upper <- cis[2,,]
      if(class(upper) !='matrix') upper <- as.matrix(upper)
      upper.CI <- split(upper, row(upper))

  }
  else
    {
      ## more complex method that takes the difference between each
      ## bootstrap and the fitted, divides by the bootstrap standard error
      ## to get a "t-distribution" that is then multiplied by the observed
      ## standard error and added to the fitted values.
      t.scores <- lapply(bsci, function(sim, obs){
        do.call('cbind',  mapply(function(obs.t, sim.t){
          t.score.t <- (sim.t$pred - obs.t$pred)/sim.t$se.pred
          t.score.t <- matrix(t.score.t, ncol=1)
          return(t.score.t)}, obs.t=obs, sim.t=sim, SIMPLIFY=F))}, obs=obs.pred)
      
      t.scores <- do.call('abind', args=list(t.scores, along=3))

      uci <- apply(t.scores, c(2,1), quantile, alpha/2, na.rm=T)
      lci <- apply(t.scores, c(2,1), quantile, 1-alpha/2, na.rm=T)
      uci <-  split(uci, row(uci))
      lci <- split(lci, row(lci))

      
      CIs <- mapply(function(obs, ucl, lcl){
        lower.CI <- obs$pred - lcl*obs$se.pred
        upper.CI <- obs$pred - ucl*obs$se.pred
        return(list(LCL=lower.CI, UCL=upper.CI))
      },  obs=obs.pred, ucl=uci, lcl=lci, SIMPLIFY=FALSE)
      
      lower.CI <- lapply(CIs, function(x) x$LCL)
      upper.CI <-  lapply(CIs, function(x) x$UCL)
    }
  estimator <-  lapply(obs.pred, function(x) x$pred)

  return(list(lmK=lmmods, lmKpred=estimator,
              lower=lower.CI, upper=upper.CI))
}

###########################################
## Function to test the null hypothesis that
## a given term in the model is not important
###############################################
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


################################################################################
## function to fit mixed effects model to k functions
################################################################################
kfuncLme <- function(k, dat, weights, fixed, random, correlation,  na.action=na.omit)
{

  ## lme's weights are actually the a variance covariate,
  ##  and the inverse of the weights in lm.
  ## So, we need to take their inverse here
  weights <- 1/weights
  ## also need to make sure these weights have a mean of 1
  weights <- weights/mean(weights)
  
    ###require(nlme) ## load required packages
    k.data.frame <- data.frame(K=k, dat, wts = weights) # make data frame

    
    ## fit generalised linear model with unstructured correlation matrix at level 1
    ## using known variance covariate weights
    fixed.formula <- as.formula(paste('K ~', fixed))
    random.formula <- as.formula(paste('~', random))

    correlation <- correlation
    ## Fit model using 'try' to avoid problems of non-convergence at some
    ## distances
    lme.k<- try(eval(substitute(lme(fixed=fixed.formula, random=random.formula,
                                    correlation=correlation, weights=varFixed(~wts),
                                    data=k.data.frame, method="REML",
                                    na.action=naaction),
                                list(fixed.formula=fixed.formula,
                                     random.formula=random.formula,
                                     naaction=na.action
                                     ))), silent=T)
    if(class(lme.k) =='try-error')
        lme.k <- NULL

    return(lme.k)
}

################################################################################
## Function that repeats the k function mixed effects model from kfunclme
## at all distances - not used in paper code but written afterwards to simplify
## fitting process
################################################################################
lmeHyperframe <- function(hyper, r , fixed, random, correlation =NULL,
                          correction='border', weights.type, computeK=TRUE,
                          minsamp=NA, printwarnings=TRUE){
  if(min(r) > 0)
    r <- c(0, r)
  if(!all(c('K', 'wts') %in% names(hyper)))
    if(!computeK)
      stop("hyperframe object must include 'K' and 'wts' if 'computeK is FALSE")
  else
    if(!('pppx' %in% names(hyper)))
      stop("hyperframe object must include 'pppx' entry if 'computeK' is TRUE")
    else
      hyper <- constructHyperframe(hyper=hyper, r=r, correction=correction,
                          weights.type=weights.type)
  
  if(!all(c('pppx') %in% names(hyper)))
    {
      stop("hyperframe object must include 'pppx' element")
    }

  
  ## find species that have sufficient individuals (define with minwt)
  ## across the plot for analysis - do this by working out whether there are at
  ## least minsamp 
  if(!is.na(minsamp))
    {
      sp.keep <-sapply(
                       with.hyperframe(hyper,
                                       list(kfunc.weights.calc(pppx, r=K$r,
                                                               correction=correction,
                                                               type='nx'))),
                       function(x) all(unlist(x[r]) >= minsamp))
      
      removed.species <- row.names(hyper)[!sp.keep]
      
      if(printwarnings & length(removed.species) > 0)
        warning(paste("Removed", length(removed.species),
                      "species with insufficient numbers"))
      
      ## select these species
      hyper <- hyper[sp.keep,]
    }
  else
    removed.species <- NULL
  
  ## Do not model distances where the variance is 0
  dist.keep <-  (apply(sapply(hyper$K, function(K) K[[correction]]), 1,
                      function(x) var(x)) > 0)
  if(printwarnings)
    warning(paste('Not modelling K at distances ',
                  paste(r[!dist.keep], collapse=', '),
                  "due to zero variance"))
  modr <-  match(r[dist.keep], r )

  ## fit the models
  
  kmods <- sapply(modr, function(i)
                  {
                    kfuncLme(k=sapply(hyper$K, function(k) k[[correction]][i]),
                             dat=as.data.frame(hyper, warn=FALSE),
                             fixed=fixed, random=random, correlation=correlation,
                             weights=sapply(hyper$wts, function(w)
                               unlist(w)[i]))
                  }, simplify=FALSE)
  names(kmods) <- r[dist.keep]
  attr(kmods, 'removed.species') <- removed.species
  class(kmods) <- 'kfunctionlme'
  return(kmods)}

## method to print kfunction lm objects
print.kfunctionlme <- function(x, ...){
  dists <- as.numeric(names(x))
  cat('linear model fitted to k function \n')
  cat("distances modelled:\n", length(dists), 'distances\n',
      'range =', range(dists), '\n')
  cat('Fixed effects\n')
  print(rbind(dists, sapply(x, fixef)))
}


################################################################################
## Function that homogenises the residuals from a linear mixed effects model
## to make them exchangable. Handles multiple levels of random effects and
## continuous dependence structures (spatial or temporal)
################################################################################
residual.homogenise.lme <- function(mod){
  if(is.null(mod))
      return(NULL)
  ## Extract grouping levels for each level
  grps <-  lapply(1:length(ranef(mod)), function(i)
                  getGroups(mod, level=i))
  ## pull out the variance covariate (which is the inverse of the
  ## weights, but confusingly labelled wts - sorry!
  if(!is.null(mod$modelStruct$varStruct))
    wts <-   getCovariate(mod$modelStruct$varStruct)[order(order(getGroups(mod)))] 
  ##    wts <- getData(mod)$wts
      
  else wts <- rep(1,mod$dims$N) ## if no weights, set them all as equal -
                                        # this shouldn't happen!
  ## Pull out the variances and co-variances of blups at all levels
  var.comps <- lapply(as.matrix(mod$modelStruct$reStruct), function(x, sig)
                      x*sig^2,
                      sig=mod$sigma)

  ## for the random effects
  Unew <- mapply(function(u, Sigma){
        ## correct for mean
        u <- as.matrix(u)
        u <- apply(u, 2, function(x) x-mean(x))
        ##  calcuate empirical variance-covarance
        S <- (t(u) %*% u)/(NROW(u))
        ## LOWER cholesky decomposition of Sigma
        Lsig <- t(chol(Sigma))
        Ls <- t(chol(S))
        A <-   t(Lsig %*% solve(Ls))
        unew <- u %*% A
        return(unew)
    }, u= ranef(mod), Sigma=var.comps, SIMPLIFY=FALSE)

    ## extract the variance-covariance of the residuals
  Clarge <- diag(mod$dims$N) ## if no correlation just an identity matrix
  ## populate with the correct values extracted from the model if there is
  ## a correlation
    if(!is.null(mod$modelStruct$corStruct)){
        for(i in levels(grps)) {
            Cj <- cov2cor(getVarCov(mod, individuals=i, type='conditional')[[1]])
            Clarge[grps==i, grps==i] <- Cj
        }
    }

  ## note that really we are using the inverse of the wts here,
  ## not the wts - this is because lme needs to be given 
  ## a variance covariate which is the inverse of the weights.
    ## remove correlation in residuals
    transform <- solve(t(chol(diag(sqrt(wts)) %*% Clarge %*%
                              diag(sqrt(wts))))) 
    resids <- as.vector(transform %*% resid(mod))  
    resids <- resids - mean(resids) ## centre
    Ls1 <- sqrt(mean(resids^2)) ## empirical resid variance
    A1 <- mod$sigma/Ls1   ## ratio of modelled to observed
    residNew <- resids*A1  ## make empirical resid variance = to modelled variance 
  result <- c(level1resids=list(residNew), Unew)
  attr(result, 'zmat') <- Clarge
  return(result)
}

################################################################################
# Function that samples transformed level 1 residuals with replacement
## to supply bootstrap replicates
# is called by estimator.bootstrap.distribution
# is called by test.bootstrap.distribution.lin
###############################################################################
residual.randomise.lme <- function(mods, resids)
{
  ## extract the level 1 residuals
  level1.resid <- lapply(resids, function(x) x[['level1resids']])

  ## mods is a list of models, which also contain some extra info.
  ##level1.resid is list of homogenised residuals from models for each distance
  N <- max(sapply(level1.resid, length)) ## extract sample size

  indx <- sample(1:N, replace=T) ## make sample for bootstrap to be replicated
                                        # at all distances

  level1.resid.r <- lapply(level1.resid, function(x) return(x[indx])) ## resample

   ## reapply the inhomogeneities to the data
  Cmat <- sapply(resids, function(x) attr(x, 'zmat'), simplify=FALSE)

  level1.resid.raw.r <- mapply(function(mod, res, Clarge) {
    if(is.null(mod))
      return(NULL)
    else
      {
        ## extract variance covariate
        wts <-   getCovariate(mod$modelStruct$varStruct)[order(order(getGroups(mod)))]
        transform1 <- t(chol( diag(sqrt(wts)) %*% Clarge %*% diag(sqrt(wts)) ))
        level1.res.raw.r <- as.vector(transform1 %*% res) ## apply back transform
        return(level1.res.raw.r)
      }
  }, mod=mods, res=level1.resid.r, Clarge=Cmat, SIMPLIFY=F)

  ## randomise the random effects
  ## first get sample for each level
  samp <- lapply(mods, function(mod){
    if(is.null(mod))
      return(NULL)
    else
      {
        re <- ranef(mod)
        return(lapply(re, function(rj){
          sample(1:NROW(rj), replace=T)}))
      }})
  ## At this point I can think of no good reason why the length of the
  ## vector would differ between distances
  ## this would be possible if some plots were a lot smaller perhaps, but not sure
  ## about the logic of including plots where the maximum distance analysed is
  ## larger than possible with a plot. For now therefore, just taking the first
  ## set of randomisations, but this might need to change in future.
  samp <- samp[[1]]

  ranef.r <- mapply(function(mod, res, samp){
    if(is.null(mod))
      return(NULL)

    else
      {
        ## extract the homogenised random effects
        ranef.res <- res[-which(names(res)=='level1resids')] ## remove level 1
        ranef.res.r <- mapply(function(r, ord){ 
          rnew <- as.matrix(r[ord,])  ## sample ranefs according to the index
          rownames(rnew) <- rownames(r) ## return to original names.
          return(rnew)
        }, r=ranef.res, ord=samp, SIMPLIFY=FALSE )

        ## assign values to replicates
        ranef.res.new <- mapply(function(j, r, mod){
          g <- as.character(getGroups(mod, level=j))  ## pull out the group
                                        # assignments
          if(is.null(attr(r, 'rownames')))
              rownames(r) <- rownames(ranef(mod, level=j)) ## if null, replace with
                                        # original ones
          return(r[g,])  ## extract ranefs for corresponding group (which may have
                                        # been switched.
        }, j=as.list(1:length(ranef.res)), r=ranef.res.r,
                              MoreArgs=list(mod=mod),
                              SIMPLIFY=FALSE)
    }
  }, mod=mods, res=resids, MoreArgs=list(samp=samp), SIMPLIFY=FALSE)

  ## put things together to be returned.
  resids <- mapply(function(level1, ranef) {
    list(level1.resid.raw.r=level1, ranef.r=ranef)
  }, level1=level1.resid.raw.r, ranef=ranef.r, SIMPLIFY=F)
  return(resids)
}

################################################################################
## Function that refits an lme model to the bootstrapped K function data
###############################################################################
refit.lmek <- function(mod, res.r){
 if(!is.null(mod)){
   ## sum up the random effects and residuals across all levels
    summed.res <- res.r$level1.resid.raw.r +
      Reduce('+', res.r$ranef.r)
    ## pull out the original data from the model
    k.data.frame <- getData(mod)
    ## bootstrap replicate of K function added to data
    k.data.frame$K.r <- fitted(mod) + summed.res
    ## refit model with new response
    mod.new <- try(update(mod, K.r~., data=k.data.frame,
                          correlation=mod$modelStruct$corStruct),
                   silent=TRUE)
    return(mod.new)
  }
  else  return(NULL)
}
#################################################################################
## Function to perform the bootstrapping on lme models and get out
## confidence intervals.
################################################################################
bootstrap.parallel.lme <- function(mods, resids, lin.comb.Ct, nboot,
                                   ncore=1, cltype='PSOCK', iseed=NULL)
{
  ## mod is the list of models
  ##lin.comb.Ct is the linear combintaion of the fixed parameters
  ##that is of interest
  ##nboot is the number of bootstrap simulations
  ## resids are the homogenised and variance corrected residuals and ranefs
  ###require(parallel) ## load the parallel package if needed
  
    cl <- makeCluster(ncore, type=cltype) ## make connections
    RNGkind("L'Ecuyer-CMRG")
    clusterSetRNGStream(cl = cl, iseed = iseed)


    on.exit({stopCluster(cl); print('clusters closed on exit')}) ## close connectons on
                                        # exit of function
  ## export all the functions and objects to the clusters
  clusterExport(cl, list('mods', 'resids', 'lin.comb.Ct', 'residual.randomise.lme',
                         'refit.lmek', 'getpars'), envir=environment())
  ##load nlme on all remote cores
  clusterEvalQ(cl, library(nlme))
  
  pars <- parSapply(cl, 1:nboot, function(i, mods, resids, lin.comb.Ct)
                      {

                        do.again <- TRUE ## to repeat if model doesn't converge

                        while(do.again){
                          ## get a back transformed bootstrap replicate of the
                                        #level 1 residuals and ranefs
                          res.rep <- residual.randomise.lme(mods=mods, resids=resids)
                          ## refit the models
                          newmods <- mapply(refit.lmek, mod=mods, res.r=res.rep,
                                            SIMPLIFY=FALSE)
                          ## set do.again to TRUE if the model didn't converge
                          do.again <- any(sapply(newmods, function(x) {
                            if(is.null(x))
                              return(FALSE)
                            else if(class(x)=='try-error')
                              return(TRUE)
                            else
                              return(any(is.na(fixef(x))))}))

                          if(do.again)
                            print('repeating sample')
                        }
                        ## pull out the parameters from the bootstrapped model
                        pars <- sapply(newmods, getpars, lin.comb.Ct=lin.comb.Ct,
                                       simplify=FALSE)
                        return(pars)
                      },  mods=mods, resids=resids, lin.comb.Ct=lin.comb.Ct,
                    simplify=FALSE)
  return(pars)
}
################################################################################
## Utility function to extract useful parameters from a model given a
## design matrix for the data to predict from.
################################################################################
getpars <- function(mod, lin.comb.Ct) {
  if(!is.null(mod)){
    beta.r <- fixef(mod)  ## fixed effecgts
    vcov.r <- vcov(mod) ## variacne covariance of fixed effects
    est.Kmean.r <- as.vector(lin.comb.Ct %*% fixef(mod)) ## predicted values
    est.Kse.r <- sqrt(diag(lin.comb.Ct %*% vcov(mod) %*%
                           t(lin.comb.Ct))) ## standard errorsof the predicted vals

  }
  ## if model failed, then just make everything NA
  else{
    beta.r <- rep(NA, ncol(lin.comb.Ct))
    vcov.r <- matrix(NA, ncol=length(beta.r),
                     nrow=length(beta.r))
    est.Kmean.r <-  rep(NA, nrow(lin.comb.Ct))
    est.Kse.r <- rep(NA, nrow(lin.comb.Ct))
  }
  return(list(pred.r=est.Kmean.r, se.pred.r=est.Kse.r,
              beta.r=beta.r, vcov.r=vcov.r))
}

################################################################################
## A function that wraps up the bootstrap function and returns the confidence
## intervals for the predictions and paramters
################################################################################
bootstrap.t.CI.lme <- function(mods, lin.comb.Ct, nboot, alpha, ncore=1, cltype='PSOCK',
                               transform=NULL)
{
  ##make sure we have the necessary packages
  ###require('abind')
  ###require('parallel')
    ##mods are the models
    ##lin.comb.Ct is the linear combination of the fixed parameters
                                        #that is of interest
    ##nboot is the number of bootstrap simulations
    ##(1-alpha) is the confidence level

    ## calculate the homogenised and corrected residuals and ranefs.
  resids <- lapply(mods, residual.homogenise.lme)
  
    ##calls estimator.bootstrap.distribution to supply bootstrap distribution of
    ##parameter estimates of interest and their standard errors
  boot.esti <- bootstrap.parallel.lme(mods=mods, resids=resids,
                                      lin.comb.Ct=lin.comb.Ct,
                                      nboot=nboot, ncore=ncore, cltype=cltype)
  
  ## Pull out the parameter estimates and predictions with SEs for fitted model
  sample.esti <- sapply(mods, getpars, lin.comb.Ct=lin.comb.Ct, simplify=FALSE)
  ## pull out the parameters and predictions with SEs for each bootstrap rep
  ## boot.pars <- sapply(boot.esti, function(r) # old version
  ##                       sapply(r, function(x) x$beta.r), simplify=FALSE)
  boot.pars <- sapply(boot.esti, function(r, est){
    mapply(function(x, est) {
      t.sim <- (x$beta.r- est$beta.r)/sqrt(diag(x$vcov.r))
      est$beta.r - t.sim * sqrt(diag(est$vcov.r))
    },  x=r, est=est, SIMPLIFY=TRUE)
  },est=sample.esti, simplify=FALSE)
  
  ## pull ou the cis of the parameter estimates
  boot.fix.cis <- apply(do.call('abind', args=list(what=boot.pars, along=3)),
                        c(2, 1), quantile, c(alpha/2, 1-alpha/2))

  sample.fix.cis <- aperm(sapply(sample.esti, function(x) x$beta.r), c(2,1))
  sample.fix.cis <- array(sample.fix.cis, dim=c(1, dim(sample.fix.cis)))
  
  modelpars <- abind(list('estimate'=sample.fix.cis, boot.fix.cis), along=1)

  ##construct the 't-distributions' for the predictions  
  t.score <- lapply(boot.esti, function(bootsamp, obssamp)
      {
          as.matrix(mapply(function(sim, obs){
                     t.r <- (sim$pred.r - obs$pred.r)/sim$se.pred.r
                     return(t.r)
                 },  bootsamp, obssamp))
      }, obssamp=sample.esti)
  
  ## turn these into a matrix
  t.score <- do.call('abind', args=list(what=t.score, along=3))
  uci <- apply(t.score, c(2, 1), quantile, alpha/2, na.rm=T) ## upper CI
  lci <- apply(t.score, c(2, 1), quantile, 1-alpha/2, na.rm=T) ##  lower CI

  ## turn into a list
  uci <- split(uci, row(uci))
  lci <- split(lci, row(lci))
  
  CIs <- mapply(function(obs, ucl, lcl){
    lower.CI <- obs$pred.r - lcl*obs$se.pred.r ## subtract t for lower limit
    upper.CI <- obs$pred.r- ucl*obs$se.pred.r ## subtract t for upper limit
    return(list(LCL=lower.CI, UCL=upper.CI))
  }, obs=sample.esti, ucl=uci, lcl=lci, SIMPLIFY=FALSE )
  
  # Extract the required parameters
  estimator <- lapply(sample.esti, function(x) x$pred.r)
  lower.CI <- lapply(CIs, function(x) x$LCL)
  upper.CI <- lapply(CIs, function(x) x$UCL)
  ## organise into return object
  bootstrap.CI <- list(estimator = estimator, lower = lower.CI, upper = upper.CI,
                       modelpars=modelpars )
  return(bootstrap.CI)
}

################################################################################
## Function that uses bootstrapping to compare two nested models
################################################################################
compare.mods.bootstrap <- function(modH0, modH1, res.r){
 if(!is.null(modH1)){
     summed.res <- res.r$level1.resid.raw.r +
         Reduce('+', res.r$ranef.r)

     k.data.frame <- getData(modH1) ## does not matter which model as data is
                                        #identical
    ## bootstrap replicate of K function
     ## formed by adding the expectation of the null model
     ## and the bootstrapped residuals

     k.data.frame$K.r <- fitted(modH0) + summed.res
     modH0.new <- try(update(modH0, K.r~., data=k.data.frame, method='ML',
                          correlation=modH0$modelStruct$corStruct),
                   silent=TRUE)
     modH1.new <- try(update(modH1, K.r~., data=k.data.frame, method='ML',
                              correlation=modH0$modelStruct$corStruct),
                      silent=TRUE)
     stat <- try((-2)*logLik(modH0.new) - (-2)*logLik(modH1.new), silent=TRUE)
     return(stat)
   }
 else  return(NULL)
}

#######################################################################
## Function that works as a wrapper to compare two models using
## bootstrapping and provides test statistics
#######################################################################
bootstrap.compare.lme <- function (mods, term, dists, nboot, ncore,
                                   cltype='PSOCK', iseed=NULL) 
{

    ## testdists <-  dist
    
    ## if(class(dists) == 'list'){
    ##     dists <- unique(do.call('c', dists))
    ##     dists[order[dists]]
    ## } ## this is just a segment that will eventually lead to more efficient code.
        
    if (!all(as.character(dists) %in% names(mods))) 
        stop("Some test distances have not been modelled")
    
    mods <- mods[as.character(dists)]
    
    modsH1 <- lapply(mods, function(mod) {
        if (!is.null(mod)) {
            k.data.frame <- getData(mod)
            mod <- try(update(mod, method = "ML", correlation = mod$modelStruct$corStruct), 
                       silent = TRUE)
            if (class(mod) == "try-error") {
                warning("ML full models  did not converge at some distances")
                mod <- NULL
            }
        }
        else mod <- NULL
        return(mod)
    })
    modsH0 <- lapply(modsH1, function(mod, term) {
        if (!is.null(mod)) {
            k.data.frame <- getData(mod)
            mod0 <- try(update(mod, paste("~.-", term), correlation = mod$modelStruct$corStruct), 
                        silent = TRUE)
            if (class(mod0) == "try-error") {
                warning("Null model fit did not converge at some distances")
                mod0 <- NULL
            }
        }
        else mod0 <- NULL
        return(mod0)
    }, term = term)
    
    ## remove models that didn't converge
    badmods <-  mapply(function(m0, m1) is.null(m0) | is.null(m1),
                       m0=modsH0, m1=modsH1)
    if(any(badmods))
        warning(paste('removed model for distances:', dists[badmods]))
    
    modsH0 <-  modsH0[!badmods]
    modsH1 <- modsH1[!badmods]
    
    obs.stat <- mapply(function(mod, mod0) {
        if (is.null(mod) | is.null(mod0)) 
            return(NULL)
        else (-2) * logLik(mod0) - (-2) * logLik(mod)
    }, mod = modsH1, mod0 = modsH0)
    
    cl <- makeCluster(ncore, type = cltype)
    RNGkind("L'Ecuyer-CMRG")
    clusterSetRNGStream(cl = cl, iseed = iseed)
    on.exit({
        stopCluster(cl)
        print("clusters closed on exit")
    })
    clusterEvalQ(cl, library(nlme))
    clusterEvalQ(cl, library(ReplicatedPointPatterns))
    
    clusterExport(cl, list("modsH1", "modsH0"), ## "residual.homogenise.lme", 
                  ##"residual.randomise.lme", "compare.mods.bootstrap"),           
                  envir = environment())
    
    boot.stat <- parSapply(cl, 1:ceiling(nboot * 1.2), function(i, 
                                                                modsH0, modsH1) {
        resids.H1 <- lapply(modsH1, residual.homogenise.lme)
        do.again <- TRUE
        while (do.again) {
            resids.resamp <- residual.randomise.lme(mods = modsH0, 
                                                    resids = resids.H1)
            bootstrap.stat <- mapply(compare.mods.bootstrap, 
                                     modsH0, modsH1, resids.resamp)
            do.again <- any(sapply(bootstrap.stat, function(x) {
                if (is.null(x)) return(FALSE)
                else return(class(x) == "try-error")
            }))
        }
        bootstrap.stat[sapply(bootstrap.stat, is.null)] <- NA
        return(bootstrap.stat)
    }, modsH0 = modsH0, modsH1 = modsH1, simplify = FALSE)
    
    goodsims <- sapply(boot.stat, function(x) all(is.numeric(x)))
    
    if (sum(goodsims) < nboot)
        warning(paste("Only ", sum(goodsims), "completed simulations - increase iterations?"))
    boot.stat <- boot.stat[goodsims][1:nboot]

    D.obs <- sum(obs.stat[as.character(dists[!badmods])])

    D.boot <- sapply(boot.stat[!sapply(boot.stat, is.null)], 
                     function(x, d) {
                         return(sum(x[as.character(d)]))
                     }, d = dists[!badmods], simplify = TRUE)
    
    p.val <- (sum(D.obs < D.boot) + 1)/(1 + length(D.boot))
    result <- list(D = D.obs, p = p.val, D.boot = D.boot, term = term, dists=dists[!badmods])
    class(result) <-  'kfunctionlmeanova'
    return(result)
}

## function to print kfunctionlmeanova object
print.kfunctionlmeanova <- function(obj){
    print(obj[c('D', 'p', 'term', 'dists')])
}
          

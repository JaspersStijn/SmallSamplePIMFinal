# Function as obtained from Amorim et al. (2018)
# Get confidence limits based on empirical likelihood

ci.EL_test = function(pim.fit, step=1, ELmethod="JEL", ProfileEL=FALSE, stat=qchisq(.95,length(coef.pim)), br.pim=FALSE, q=1, max_iter=1000, useMean=TRUE)
{
  available.ELmethods <- c("JEL", "AJEL", "AJELb")
  if (!(ELmethod %in% available.ELmethods))
    stop("Empirical Likelihood method not avaiable. Options are: AJEL or AJELb. JEL is the default choice")

  ## First get things that don't change
  link <- link(pim.fit)
  data <- as.data.frame(penv(pim.fit))
  form <- formula(formula(pim.fit))
  N  <- nrow(data)
  PO <- response(pim.fit)
  Z  <- model.matrix(pim.fit)
  coef.pim <- coef(pim.fit)

  if ( (length(coef.pim)>1 & !ProfileEL) | (ProfileEL & q>1) )
    stop("Only confidence intervals: confidence regions not available")

  if (br.pim) {
    NewScoreFunc <- function(Zbeta, Z, PO, link) BRScoreFunc(Zbeta, Z, PO, link)
    FEVenv <- new.pim.env(data)
    vcov.tmp <- NULL
    pim.fit.br <- pim.fit(Z,PO,construct = BRScore, link = link,
                          estim = 'estimator.nleqslv', start = coef(pim.fit), penv = FEVenv)
    coef.pim <- pim.fit.br$coefficients
  } else {
    NewScoreFunc <- function(Zbeta, Z, Y, link) pim:::U.sandwich(Zbeta, Z, Y, link)$U
  }


  ## Check if profile empirical is required. If so, user must provide q (the number of parameters of interest - first q). For example, if interested in the treatment effect only after adjusting for p covariates, the model must be $y ~ treatment + x_1 + ... + x_p$ and $q = 1$.
  if (ProfileEL) coef.pim <- coef.pim[1:q]

  find.root <- function(theta, pim.fit=pim.fit, ELmethod=ELmethod, ProfileEL=ProfileEL, br.pim=br.pim, stat=stat, useMean=useMean)
    emp.lik.func(pim.fit, theta, ELmethod, ProfileEL, br.pim, useMean) - stat

  coef.pim.aux <- coef.pim
  at.coef      <- emp.lik.func(pim.fit, coef.pim, ELmethod, ProfileEL, br.pim, useMean) - stat
  at.coef.aux  <- at.coef
  counter      <- 0

  while(sign(at.coef) == sign(at.coef.aux) & counter < max_iter) {
    sin.coef <- sign(coef.pim)
    if (coef.pim == 0) sin.coef = 1
    coef.pim.aux <- coef.pim.aux + sin.coef*step
    at.coef.aux  <- emp.lik.func(pim.fit, coef.pim.aux, ELmethod, ProfileEL, br.pim, useMean) - stat
    counter      <- counter + 1
  }
  if (counter >= max_iter) {
    LL <- sign(coef.pim.aux)*Inf
    #print("maximum iterations reached: no convergence or step size too small")
  } else {
    LL <- uniroot(find.root, interval=c(coef.pim.aux,coef.pim), pim.fit=pim.fit,
                  ELmethod=ELmethod, ProfileEL=ProfileEL, br.pim=br.pim, stat=stat, useMean=useMean)$root
  }

  ## Do the same to get an upper bound for the CI
  coef.pim.aux <- coef.pim
  at.coef.aux  <- at.coef
  counter <- 0

  while(sign(at.coef) == sign(at.coef.aux) & counter < max_iter) {
    sin.coef <- sign(coef.pim)
    if (coef.pim == 0) sin.coef = 1
    coef.pim.aux <- coef.pim.aux - sin.coef*step
    at.coef.aux  <- emp.lik.func(pim.fit, coef.pim.aux, ELmethod, ProfileEL, br.pim, useMean) - stat
    counter      <- counter + 1
  }
  if (counter >= max_iter) {
    UL <- sign(coef.pim.aux)*Inf
    #print("maximum iterations reached: no convergence or step size too small")
  } else {
    UL <- uniroot(find.root, interval=c(coef.pim.aux,coef.pim), pim.fit=pim.fit,
                  ELmethod=ELmethod, ProfileEL=ProfileEL, br.pim=br.pim, stat=stat, useMean=useMean)$root
  }

  ci.EL.res <- sort(c(LL,UL))

  return(list("CI" = ci.EL.res,"beta"=coef.pim))
}





gee_supressed = function (formula = formula(data), id = id, data = parent.frame(),
                          subset, na.action, R = NULL, b = NULL, tol = 0.001, maxiter = 25,
                          family = gaussian, corstr = "independence", Mv = 1, silent = TRUE,
                          contrasts = NULL, scale.fix = FALSE, scale.value = 1, v4.4compat = FALSE)
{
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$R <- m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <- m$corstr <- m$Mv <- m$silent <- m$contrasts <- m$family <- m$scale.fix <- m$scale.value <- m$v4.4compat <- NULL
  if (is.null(m$id))
    m$id <- as.name("id")
  if (!is.null(m$na.action) && m$na.action != "na.omit") {
    warning("Only 'na.omit' is implemented for gee\ncontinuing with 'na.action=na.omit'")
    m$na.action <- as.name("na.omit")
  }
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  y <- as.matrix(model.extract(m, "response"))
  x <- model.matrix(Terms, m, contrasts)
  Q <- qr(x)
  if (Q$rank < ncol(x))
    stop("rank-deficient model matrix")
  N <- rep(1, length(y))
  if (dim(y)[2] == 2) {
    N <- as.vector(y %*% c(1, 1))
    y <- y[, 1]
  }
  else {
    if (dim(y)[2] > 2)
      stop("Only binomial response matrices (2 columns)")
  }
  offset <- model.extract(m, offset)
  id <- model.extract(m, id)
  if (is.null(id)) {
    stop("Id variable not found")
  }
  nobs <- nrow(x)
  p <- ncol(x)
  xnames <- dimnames(x)[[2]]
  if (is.null(xnames)) {
    xnames <- paste("x", 1:p, sep = "")
    dimnames(x) <- list(NULL, xnames)
  }
  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()
  if (!is.null(b)) {
    beta <- matrix(as.double(b), ncol = 1)
    if (nrow(beta) != p) {
      stop("Dim beta != ncol(x)")
    }
    message("user's initial regression estimate")
  }
  else {
    mm <- match.call(expand.dots = FALSE)
    mm$R <- mm$b <- mm$tol <- mm$maxiter <- mm$link <- mm$varfun <- mm$corstr <- mm$Mv <- mm$silent <- mm$contrasts <- mm$scale.fix <- mm$scale.value <- mm$id <- NULL
    mm[[1]] <- as.name("glm")
    beta <- eval(mm, parent.frame())$coef
    beta <- as.numeric(beta)
  }
  if (length(id) != length(y))
    stop("Id and y not same length")
  maxclsz <- as.integer(max(unlist(lapply(split(id, id), "length"))))
  maxiter <- as.integer(maxiter)
  silent <- as.integer(silent)
  if (length(offset) <= 1)
    offset <- rep(0, length(y))
  if (length(offset) != length(y))
    stop("offset and y not same length")
  offset <- as.double(offset)
  if (!is.null(R)) {
    Rr <- nrow(R)
    if (Rr != ncol(R))
      stop("R is not square!")
    if (Rr < maxclsz)
      stop("R not big enough to accommodate some clusters.")
  }
  else {
    R <- matrix(as.double(rep(0, maxclsz * maxclsz)), nrow = maxclsz)
  }
  links <- c("identity", "log", "logit", "inverse", "probit",
             "cloglog")
  fams <- c("gaussian", "poisson", "binomial", "Gamma", "quasi")
  varfuns <- c("constant", "mu", "mu(1-mu)", "mu^2")
  corstrs <- c("independence", "fixed", "stat_M_dep", "non_stat_M_dep",
               "exchangeable", "AR-M", "unstructured")
  linkv <- as.integer(match(c(family$link), links, -1))
  famv <- match(family$family, fams, -1)
  if (famv < 1)
    stop("unknown family")
  if (famv <= 4)
    varfunv <- famv
  else varfunv <- match(family$varfun, varfuns, -1)
  varfunv <- as.integer(varfunv)
  corstrv <- as.integer(match(corstr, corstrs, -1))
  if (linkv < 1)
    stop("unknown link.")
  if (varfunv < 1)
    stop("unknown varfun.")
  if (corstrv < 1)
    stop("unknown corstr.")
  naivvar <- matrix(rep(0, p * p), nrow = p)
  robvar <- matrix(rep(0, p * p), nrow = p)
  phi <- as.double(scale.value)
  scale.fix <- as.integer(scale.fix)
  errstate <- as.integer(1)
  tol <- as.double(tol)
  Mv <- as.integer(Mv)
  maxcl <- as.integer(0)
  if (!(is.double(x)))
    x <- as.double(x)
  if (!(is.double(y)))
    y <- as.double(y)
  if (!(is.double(id)))
    id <- as.double(id)
  if (!(is.double(N)))
    N <- as.double(N)
  modvec <- as.integer(c(linkv, varfunv, corstrv))
  if (v4.4compat)
    compatflag <- 1
  else compatflag <- 0
  z <- .C(Cgee, x, y, id, N, offset, nobs, p, modvec, Mv, estb = beta,
          nv = naivvar, rv = robvar, sc = phi, wcor = R, tol, mc = maxcl,
          iter = maxiter, silent, err = errstate, scale.fix, as.integer(compatflag))
  if (z$err != 0)
    warning("Cgee had an error (code= ", z$err, ").  Results suspect.")
  if (min(eigen(z$wcor)$values) < 0) {
    warning("Working correlation estimate not positive definite")
    z$err <- z$err + 1000
  }
  fit <- list()
  attr(fit, "class") <- c("gee", "glm")
  fit$title <- "GEE:  GENERALIZED LINEAR MODELS FOR DEPENDENT DATA"
  fit$version <- "gee S-function, version 4.13 modified 98/01/27 (1998)"
  links <- c("Identity", "Logarithm", "Logit", "Reciprocal",
             "Probit", "Cloglog")
  varfuns <- c("Gaussian", "Poisson", "Binomial", "Gamma")
  corstrs <- c("Independent", "Fixed", "Stationary M-dependent",
               "Non-Stationary M-dependent", "Exchangeable", "AR-M",
               "Unstructured")
  fit$model <- list()
  fit$model$link <- links[linkv]
  fit$model$varfun <- varfuns[varfunv]
  fit$model$corstr <- corstrs[corstrv]
  if (!is.na(match(c(corstrv), c(3, 4, 6))))
    fit$model$M <- Mv
  fit$call <- call
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  fit$contrasts <- attr(x, "contrasts")
  fit$nobs <- nobs
  fit$iterations <- z$iter
  fit$coefficients <- as.vector(z$estb)
  fit$nas <- is.na(fit$coefficients)
  names(fit$coefficients) <- xnames
  eta <- as.vector(x %*% fit$coefficients)
  fit$linear.predictors <- eta
  mu <- as.vector(family$linkinv(eta))
  fit$fitted.values <- mu
  fit$residuals <- y - mu
  fit$family <- family
  fit$y <- as.vector(y)
  fit$id <- as.vector(id)
  fit$max.id <- z$mc
  z$wcor <- matrix(z$wcor, ncol = maxclsz)
  fit$working.correlation <- z$wcor
  fit$scale <- z$sc
  fit$robust.variance <- z$rv
  fit$naive.variance <- z$nv
  fit$xnames <- xnames
  fit$error <- z$err
  dimnames(fit$robust.variance) <- list(xnames, xnames)
  dimnames(fit$naive.variance) <- list(xnames, xnames)
  fit
}

require(gee)
environment(gee_supressed) <- asNamespace('gee')
assignInNamespace("gee", gee_supressed, ns = "gee")

GEE.var.gst_new = function (formula, id, family = binomial, data, corstr = "independence")
{
  if (is.null(data$id)) {
    index <- which(names(data) == id)
    data$id <- data[, index]
  }
  init <- model.frame(formula, data)
  init$num <- 1:length(init[, 1])
  if (any(is.na(init))) {
    index <- na.omit(init)$num
    data <- data[index, ]
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  else {
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  gee.fit <- gee_supressed(formula, data = data, id = id, family = family,
                           corstr = corstr)
  beta_est <- gee.fit$coefficient
  alpha <- 0
  len <- length(beta_est)
  len_vec <- len^2
  data$id <- gee.fit$id
  cluster <- cluster.size(data$id)
  ncluster <- max(cluster$n)
  size <- cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if (is.character(corstr)) {
    var <- switch(corstr, independence = cormax.ind(ncluster),
                  exchangeable = cormax.exch(ncluster, alpha), `AR-M` = cormax.ar1(ncluster,
                                                                                   alpha), unstructured = summary(gee.fit)$working.correlation,
    )
  }
  else {
    stop("'working correlation structure' not recognized")
  }
  if (is.character(family)) {
    family <- switch(family, gaussian = "gaussian", binomial = "binomial",
                     poisson = "poisson")
  }
  else {
    if (is.function(family)) {
      family <- family()[[1]]
    }
    else {
      stop("'family' not recognized")
    }
  }
  cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)
  step <- matrix(0, nrow = cluster$n[1], ncol = cluster$n[1])
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    if (family == "binomial") {
      resid <- (y - exp(covariate %*% beta_est)/(1 + exp(covariate %*%
                                                           beta_est))) %*% t(y - exp(covariate %*% beta_est)/(1 +
                                                                                                                exp(covariate %*% beta_est)))
      B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
      diag(B) <- 1/sqrt(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)
      step <- step + B %*% resid %*% B
    }
    else {
      stop("'family' not recognized")
    }
  }
  unstr <- step/(size - len)
  step11 <- matrix(0, nrow = len, ncol = len)
  step12 <- matrix(0, nrow = len, ncol = len)
  step13 <- matrix(0, nrow = len_vec, ncol = 1)
  step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
  p <- matrix(0, nrow = len_vec, ncol = size)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "binomial") {
      A <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
      diag(A) <- exp(covariate %*% beta_est)/(1 + exp(covariate %*%
                                                        beta_est))^2
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      xy <- t(D) %*% solve(Vi) %*% sqrt(A) %*% unstr %*%
        sqrt(A) %*% solve(Vi) %*% D
      xx <- t(D) %*% solve(Vi) %*% D
      step12 <- step12 + xy
      step11 <- step11 + xx
      step13 <- step13 + vec(xy)
      p[, i] <- vec(xy)
    }
  }
  for (i in 1:size) {
    dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
    step14 <- step14 + dif
  }
  cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
  cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
    step14 %*% kronecker(solve(step11), solve(step11))
  return(list(beta=beta_est,cov.beta = diag(cov.beta), cov.var = cov.var))
}

GEE.var.kc_new=function (formula, id, family = binomial, data, corstr = "independence")
{
  if (is.null(data$id)) {
    index <- which(names(data) == id)
    data$id <- data[, index]
  }
  init <- model.frame(formula, data)
  init$num <- 1:length(init[, 1])
  if (any(is.na(init))) {
    index <- na.omit(init)$num
    data <- data[index, ]
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  else {
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  gee.fit <- gee_supressed(formula, data = data, id = id, family = family,
                           corstr = corstr)
  beta_est <- gee.fit$coefficient
  alpha <- 0
  len <- length(beta_est)
  len_vec <- len^2
  data$id <- gee.fit$id
  cluster <- cluster.size(data$id)
  ncluster <- max(cluster$n)
  size <- cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if (is.character(corstr)) {
    var <- switch(corstr, independence = cormax.ind(ncluster),
                  exchangeable = cormax.exch(ncluster, alpha), `AR-M` = cormax.ar1(ncluster,
                                                                                   alpha), unstructured = summary(gee.fit)$working.correlation,
    )
  }
  else {
    stop("'working correlation structure' not recognized")
  }
  if (is.character(family)) {
    family <- switch(family, gaussian = "gaussian", binomial = "binomial",
                     poisson = "poisson")
  }
  else {
    if (is.function(family)) {
      family <- family()[[1]]
    }
    else {
      stop("'family' not recognized")
    }
  }
  cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)
  step11 <- matrix(0, nrow = len, ncol = len)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      xx <- t(D) %*% solve(Vi) %*% D
      step11 <- step11 + xx
    }
    else {
      stop("'family' not recognized")
    }
  }
  step12 <- matrix(0, nrow = len, ncol = len)
  step13 <- matrix(0, nrow = len_vec, ncol = 1)
  step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
  p <- matrix(0, nrow = len_vec, ncol = size)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      xy <- t(D) %*% solve(Vi) %*% mat.sqrt.inv(cormax.ind(cluster$n[i]) -
                                                  D %*% solve(step11) %*% t(D) %*% solve(Vi)) %*%
        (y - exp(covariate %*% beta_est)/(1 + exp(covariate %*%
                                                    beta_est)))

      step12 <- step12 + xy %*% t(xy)
      step13 <- step13 + vec(xy %*% t(xy))
      p[, i] <- vec(xy %*% t(xy))
    }
  }
  for (i in 1:size) {
    dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
    step14 <- step14 + dif
  }
  cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
  cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
    step14 %*% kronecker(solve(step11), solve(step11))
  return(list(beta=beta_est,cov.beta = diag(cov.beta), cov.var = cov.var))
}

GEE.var.pan_new = function(formula,id,family,data,corstr="independence"){
  if (is.null(data$id)) {
    index <- which(names(data) == id)
    data$id <- data[, index]
  }
  init <- model.frame(formula, data)
  init$num <- 1:length(init[, 1])
  if (any(is.na(init))) {
    index <- na.omit(init)$num
    data <- data[index, ]
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  else {
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  gee.fit <- gee_supressed(formula, data = data, id = id, family = family,
                           corstr = corstr)
  beta_est <- gee.fit$coefficient
  alpha <- 0
  len <- length(beta_est)
  len_vec <- len^2
  data$id <- gee.fit$id
  cluster <- cluster.size(data$id)
  ncluster <- max(cluster$n)
  size <- cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if (is.character(corstr)) {
    var <- switch(corstr, independence = cormax.ind(ncluster),
                  exchangeable = cormax.exch(ncluster, alpha), `AR-M` = cormax.ar1(ncluster,
                                                                                   alpha), unstructured = summary(gee.fit)$working.correlation,
    )
  }
  else {
    stop("'working correlation structure' not recognized")
  }
  if (is.character(family)) {
    family <- switch(family, gaussian = "gaussian", binomial = "binomial",
                     poisson = "poisson")
  }
  else {
    if (is.function(family)) {
      family <- family()[[1]]
    }
    else {
      stop("'family' not recognized")
    }
  }
  cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)
  step <- matrix(0, nrow = cluster$n[1], ncol = cluster$n[1])
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    if (family == "binomial") {
      resid <- (y - exp(covariate %*% beta_est)/(1 + exp(covariate %*%
                                                           beta_est))) %*% t(y - exp(covariate %*% beta_est)/(1 +
                                                                                                                exp(covariate %*% beta_est)))
      B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
      diag(B) <- 1/sqrt(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)
      step <- step + B %*% resid %*% B
    }
    else {
      stop("'family' not recognized")
    }
  }
  unstr <- step/size
  step11 <- matrix(0, nrow = len, ncol = len)
  step12 <- matrix(0, nrow = len, ncol = len)
  step13 <- matrix(0, nrow = len_vec, ncol = 1)
  step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
  p <- matrix(0, nrow = len_vec, ncol = size)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "binomial") {
      A <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
      diag(A) <- exp(covariate %*% beta_est)/(1 + exp(covariate %*%
                                                        beta_est))^2
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      xy <- t(D) %*% solve(Vi) %*% sqrt(A) %*% unstr %*%
        sqrt(A) %*% solve(Vi) %*% D
      xx <- t(D) %*% solve(Vi) %*% D
      step12 <- step12 + xy
      step11 <- step11 + xx
      step13 <- step13 + vec(xy)
      p[, i] <- vec(xy)
    }
  }
  for (i in 1:size) {
    dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
    step14 <- step14 + dif
  }
  cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
  cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
    step14 %*% kronecker(solve(step11), solve(step11))
  return(list(beta = beta_est,cov.beta = diag(cov.beta), cov.var = cov.var))
}

GEE.var.md_new=function (formula, id, family = binomial, data, corstr = "independence")
{
  if (is.null(data$id)) {
    index <- which(names(data) == id)
    data$id <- data[, index]
  }
  init <- model.frame(formula, data)
  init$num <- 1:length(init[, 1])
  if (any(is.na(init))) {
    index <- na.omit(init)$num
    data <- data[index, ]
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  else {
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  gee.fit <- gee_supressed(formula, data = data, id = id, family = family,
                           corstr = corstr)
  beta_est <- gee.fit$coefficient
  alpha <- 0
  len <- length(beta_est)
  len_vec <- len^2
  data$id <- gee.fit$id
  cluster <- cluster.size(data$id)
  ncluster <- max(cluster$n)
  size <- cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if (is.character(corstr)) {
    var <- switch(corstr, independence = cormax.ind(ncluster),
                  exchangeable = cormax.exch(ncluster, alpha), `AR-M` = cormax.ar1(ncluster,
                                                                                   alpha), unstructured = summary(gee.fit)$working.correlation,
    )
  }
  else {
    stop("'working correlation structure' not recognized")
  }
  if (is.character(family)) {
    family <- switch(family, gaussian = "gaussian", binomial = "binomial",
                     poisson = "poisson")
  }
  else {
    if (is.function(family)) {
      family <- family()[[1]]
    }
    else {
      stop("'family' not recognized")
    }
  }
  cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)
  step11 <- matrix(0, nrow = len, ncol = len)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      xx <- t(D) %*% solve(Vi) %*% D
      step11 <- step11 + xx
    }
    else {
      stop("'family' not recognized")
    }
  }
  step12 <- matrix(0, nrow = len, ncol = len)
  step13 <- matrix(0, nrow = len_vec, ncol = 1)
  step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
  p <- matrix(0, nrow = len_vec, ncol = size)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      xy <- t(D) %*% solve(Vi) %*% solve(cormax.ind(cluster$n[i]) -
                                           D %*% solve(step11) %*% t(D) %*% solve(Vi)) %*%
        (y - exp(covariate %*% beta_est)/(1 + exp(covariate %*%
                                                    beta_est)))
      step12 <- step12 + xy %*% t(xy)
      step13 <- step13 + vec(xy %*% t(xy))
      p[, i] <- vec(xy %*% t(xy))
    }
  }
  for (i in 1:size) {
    dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
    step14 <- step14 + dif
  }
  cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
  cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
    step14 %*% kronecker(solve(step11), solve(step11))
  return(list(beta=beta_est,cov.beta = diag(cov.beta), cov.var = cov.var))
}

GEE.var.mbn_new = function(formula, id, family = binomial, data, corstr = "independence",
                           d = 2, r = 1)
{
  if (is.null(data$id)) {
    index <- which(names(data) == id)
    data$id <- data[, index]
  }
  init <- model.frame(formula, data)
  init$num <- 1:length(init[, 1])
  if (any(is.na(init))) {
    index <- na.omit(init)$num
    data <- data[index, ]
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  else {
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }

  gee.fit <- gee_supressed(formula, data = data, id = id, family = binomial,
                           corstr = corstr)
  beta_est <- gee.fit$coefficient
  alpha <- 0
  scale <- summary(gee.fit)$scale
  len <- length(beta_est)
  len_vec <- len^2
  data$id <- gee.fit$id
  cluster <- cluster.size(data$id)
  ncluster <- max(cluster$n)
  size <- cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if (is.character(corstr)) {
    var <- switch(corstr, independence = cormax.ind(ncluster),
                  exchangeable = cormax.exch(ncluster, alpha), `AR-M` = cormax.ar1(ncluster,
                                                                                   alpha), unstructured = summary(gee.fit)$working.correlation,
    )
  }
  else {
    stop("'working correlation structure' not recognized")
  }
  if (is.character(family)) {
    family <- switch(family, gaussian = "gaussian", binomial = "binomial",
                     poisson = "poisson")
  }
  else {
    if (is.function(family)) {
      family <- family()[[1]]
    }
    else {
      stop("'family' not recognized")
    }
  }
  cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)
  step11 <- matrix(0, nrow = len, ncol = len)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    ncluster = cluster$n[i]
    var1 = var[1:ncluster, 1:ncluster]
    if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), ncluster) %*%
        var1 %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                            exp(covariate %*% beta_est))^2)), ncluster)
      xx <- t(D) %*% solve(Vi) %*% D
      step11 <- step11 + xx
    }
    else {
      stop("'family' not recognized")
    }
  }
  k <- (sum(cluster$n) - 1)/(sum(cluster$n) - len) * size/(size -
                                                             1)
  delta <- ifelse(size > ((d + 1) * len), len/(size - len),
                  1/d)
  step00 <- matrix(0, nrow = len, ncol = len)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    ncluster = cluster$n[i]
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    var1 = var[1:ncluster, 1:ncluster]
    if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), ncluster) %*%
        var1 %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                            exp(covariate %*% beta_est))^2)), ncluster)
      xy <- t(D) %*% solve(Vi) %*% (y - exp(covariate %*%
                                              beta_est)/(1 + exp(covariate %*% beta_est)))
      step00 <- step00 + xy %*% t(xy)
    }
  }
  xi <- pmax(r, sum(diag(solve(step11) %*% step00))/len)
  step12 <- matrix(0, nrow = len, ncol = len)
  step13 <- matrix(0, nrow = len_vec, ncol = 1)
  step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
  p <- matrix(0, nrow = len_vec, ncol = size)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    ncluster = cluster$n[i]
    var1 = var[1:ncluster, 1:ncluster]
    if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), ncluster) %*%
        var1 %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                            exp(covariate %*% beta_est))^2)), ncluster)
      xy <- t(D) %*% solve(Vi) %*% (k * (y - exp(covariate %*%
                                                   beta_est)/(1 + exp(covariate %*% beta_est))) %*%
                                      t(y - exp(covariate %*% beta_est)/(1 + exp(covariate %*%
                                                                                   beta_est))) + delta * xi * Vi) %*% solve(Vi) %*%
        D
      step12 <- step12 + xy
      step13 <- step13 + vec(xy)
      p[, i] <- vec(xy)
    }
  }
  for (i in 1:size) {
    dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
    step14 <- step14 + dif
  }
  cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
  cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
    step14 %*% kronecker(solve(step11), solve(step11))
  return(list(beta = beta_est,cov.beta = diag(cov.beta), cov.var = cov.var))
}

GEE.var.fg_new = function (formula, id, family = binomial, data, corstr = "independence",
                           b = 0.75)
{
  if (is.null(data$id)) {
    index <- which(names(data) == id)
    data$id <- data[, index]
  }
  init <- model.frame(formula, data)
  init$num <- 1:length(init[, 1])
  if (any(is.na(init))) {
    index <- na.omit(init)$num
    data <- data[index, ]
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  else {
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  gee.fit <- gee_supressed(formula, data = data, id = id, family = family,
                           corstr = corstr)
  beta_est <- gee.fit$coefficient
  alpha <- 0
  len <- length(beta_est)
  len_vec <- len^2
  data$id <- gee.fit$id
  cluster <- cluster.size(data$id)
  ncluster <- max(cluster$n)
  size <- cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if (is.character(corstr)) {
    var <- switch(corstr, independence = cormax.ind(ncluster),
                  exchangeable = cormax.exch(ncluster, alpha), `AR-M` = cormax.ar1(ncluster,
                                                                                   alpha), unstructured = summary(gee.fit)$working.correlation,
    )
  }
  else {
    stop("'working correlation structure' not recognized")
  }
  if (is.character(family)) {
    family <- switch(family, gaussian = "gaussian", binomial = "binomial",
                     poisson = "poisson")
  }
  else {
    if (is.function(family)) {
      family <- family()[[1]]
    }
    else {
      stop("'family' not recognized")
    }
  }
  cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)
  step11 <- matrix(0, nrow = len, ncol = len)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    ncluster = cluster$n[i]
    var1 = var[1:ncluster, 1:ncluster]
    if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), ncluster) %*%
        var1 %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                            exp(covariate %*% beta_est))^2)), ncluster)
      xx <- t(D) %*% solve(Vi) %*% D
      step11 <- step11 + xx
    }
    else {
      stop("'family' not recognized")
    }
  }
  step12 <- matrix(0, nrow = len, ncol = len)
  step13 <- matrix(0, nrow = len_vec, ncol = 1)
  step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
  p <- matrix(0, nrow = len_vec, ncol = size)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    ncluster = cluster$n[i]
    var1 = var[1:ncluster, 1:ncluster]
    if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), ncluster) %*%
        var1 %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                            exp(covariate %*% beta_est))^2)), ncluster)
      xx <- t(D) %*% solve(Vi) %*% D
      Qi <- xx %*% solve(step11)
      Ai <- diag((1 - pmin(b, diag(Qi)))^(-0.5))

      xy <- Ai %*% t(D) %*% solve(Vi) %*% (y - exp(covariate %*%
                                                     beta_est)/(1 + exp(covariate %*% beta_est)))
      step12 <- step12 + xy %*% t(xy)
      step13 <- step13 + vec(xy %*% t(xy))
      p[, i] <- vec(xy %*% t(xy))
    }
  }
  for (i in 1:size) {
    dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
    step14 <- step14 + dif
  }
  cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
  cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
    step14 %*% kronecker(solve(step11), solve(step11))
  return(list(beta=beta_est,cov.beta = diag(cov.beta), cov.var = cov.var))
}

GEE.var.mk_new = function (formula, id, family = binomial, data, corstr = "independence")
{
  if (is.null(data$id)) {
    index <- which(names(data) == id)
    data$id <- data[, index]
  }
  init <- model.frame(formula, data)
  init$num <- 1:length(init[, 1])
  if (any(is.na(init))) {
    index <- na.omit(init)$num
    data <- data[index, ]
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  else {
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  gee.fit <- gee_supressed(formula, data = data, id = id, family = family,
                           corstr = corstr, scale.value = 1)
  beta_est <- gee.fit$coefficient
  alpha <- 0
  len <- length(beta_est)
  len_vec <- len^2
  data$id <- gee.fit$id
  cluster <- cluster.size(data$id)
  ncluster <- max(cluster$n)
  size <- cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if (is.character(corstr)) {
    var <- switch(corstr, independence = cormax.ind(ncluster),
                  exchangeable = cormax.exch(ncluster, alpha), `AR-M` = cormax.ar1(ncluster,
                                                                                   alpha), unstructured = summary(gee.fit)$working.correlation,
    )
  }
  else {
    stop("'working correlation structure' not recognized")
  }
  if (is.character(family)) {
    family <- switch(family, gaussian = "gaussian", binomial = "binomial",
                     poisson = "poisson")
  }
  else {
    if (is.function(family)) {
      family <- family()[[1]]
    }
    else {
      stop("'family' not recognized")
    }
  }
  cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)
  step11 <- matrix(0, nrow = len, ncol = len)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      xx <- t(D) %*% solve(Vi) %*% D
      step11 <- step11 + xx
    }
    else {
      stop("'family' not recognized")
    }
  }
  step12 <- matrix(0, nrow = len, ncol = len)
  step13 <- matrix(0, nrow = len_vec, ncol = 1)
  step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
  p <- matrix(0, nrow = len_vec, ncol = size)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      xy <- t(D) %*% solve(Vi) %*% (y - exp(covariate %*%
                                              beta_est)/(1 + exp(covariate %*% beta_est)))
      step12 <- step12 + xy %*% t(xy)
      step13 <- step13 + vec(xy %*% t(xy))
      p[, i] <- vec(xy %*% t(xy))
    }
  }
  for (i in 1:size) {
    dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
    step14 <- step14 + dif
  }
  scale <- size/(size - len)
  cov.beta <- scale * solve(step11) %*% (step12) %*% solve(step11)
  cov.var <- scale^2 * size/(size - 1) * kronecker(solve(step11),
                                                   solve(step11)) %*% step14 %*% kronecker(solve(step11),
                                                                                           solve(step11))
  return(list(beta=beta_est,cov.beta = diag(cov.beta), cov.var = cov.var))
}

GEE.var.wl_new = function (formula, id, family = binomial, data, corstr = "independence")
{
  if (is.null(data$id)) {
    index <- which(names(data) == id)
    data$id <- data[, index]
  }
  init <- model.frame(formula, data)
  init$num <- 1:length(init[, 1])
  if (any(is.na(init))) {
    index <- na.omit(init)$num
    data <- data[index, ]
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  else {
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  gee.fit <- gee_supressed(formula, data = data, id = id, family = family,
                           corstr = corstr)
  beta_est <- gee.fit$coefficient
  alpha <- 0
  len <- length(beta_est)
  len_vec <- len^2
  data$id <- gee.fit$id
  cluster <- cluster.size(data$id)
  ncluster <- max(cluster$n)
  size <- cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if (is.character(corstr)) {
    var <- switch(corstr, independence = cormax.ind(ncluster),
                  exchangeable = cormax.exch(ncluster, alpha), `AR-M` = cormax.ar1(ncluster,
                                                                                   alpha), unstructured = summary(gee.fit)$working.correlation,
    )
  }
  else {
    stop("'working correlation structure' not recognized")
  }
  if (is.character(family)) {
    family <- switch(family, gaussian = "gaussian", binomial = "binomial",
                     poisson = "poisson")
  }
  else {
    if (is.function(family)) {
      family <- family()[[1]]
    }
    else {
      stop("'family' not recognized")
    }
  }
  m <- model.frame(formula, data)
  mat <- as.data.frame(model.matrix(formula, m))
  mat$subj <- rep(unique(data$id), cluster$n)
  cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)
  step01 <- matrix(0, nrow = len, ncol = len)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      xx <- t(D) %*% solve(Vi) %*% D
      step01 <- step01 + xx
    }
    else {
      stop("'family' not recognized")
    }
  }
  step <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "gaussian") {
      resid <- solve(cormax.ind(cluster$n[i]) - covariate %*%
                       solve(step01) %*% t(covariate) %*% solve(var_i)) %*%
        (y - covariate %*% beta_est)
      step <- step + resid %*% t(resid)
    }
    else if (family == "poisson") {
      B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
      diag(B) <- 1/sqrt(exp(covariate %*% beta_est))
      D <- mat.prod(covariate, exp(covariate %*% beta_est))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est))),
                 cluster$n[i]) %*% var_i %*% diag(sqrt(c(exp(covariate %*%
                                                               beta_est))), cluster$n[i])
      resid <- B %*% solve(cormax.ind(cluster$n[i]) - D %*%
                             solve(step01) %*% t(D) %*% solve(Vi)) %*% (y -
                                                                          exp(covariate %*% beta_est))
      step <- step + resid %*% t(resid)
    }
    else if (family == "binomial") {
      B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
      diag(B) <- 1/sqrt(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      resid <- B %*% solve(cormax.ind(cluster$n[i]) - D %*%
                             solve(step01) %*% t(D) %*% solve(Vi)) %*% (y -
                                                                          exp(covariate %*% beta_est)/(1 + exp(covariate %*%
                                                                                                                 beta_est)))
      step <- step + resid %*% t(resid)
    }
  }
  unstr <- step/size
  step11 <- matrix(0, nrow = len, ncol = len)
  step12 <- matrix(0, nrow = len, ncol = len)
  step13 <- matrix(0, nrow = len_vec, ncol = 1)
  step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
  p <- matrix(0, nrow = len_vec, ncol = size)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "binomial") {
      B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
      diag(B) <- exp(covariate %*% beta_est)/(1 + exp(covariate %*%
                                                        beta_est))^2
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 +
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 +
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      xy <- t(D) %*% solve(Vi) %*% sqrt(B) %*% unstr %*%
        sqrt(B) %*% solve(Vi) %*% D
      xx <- t(D) %*% solve(Vi) %*% D
      step12 <- step12 + xy
      step11 <- step11 + xx
      step13 <- step13 + vec(xy)
      p[, i] <- vec(xy)
    }
  }
  for (i in 1:size) {
    dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
    step14 <- step14 + dif
  }
  cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
  cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
    step14 %*% kronecker(solve(step11), solve(step11))
  return(list(beta=beta_est,cov.beta = diag(cov.beta), cov.var = cov.var))
}


test_firth = function(formula = formula(data), id = id, data = parent.frame(), corstr="independence",
                      est_dispersion = TRUE,iterlim=50){
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$R <- m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <-
    m$corstr <- m$Mv <- m$silent <- m$contrasts <-
    m$family <- m$scale.fix <- m$scale.value <- m$v4.4compat <- NULL
  m$est_dispersion <- NULL
  if(is.null(m$id)) m$id <- as.name("id")
  if(!is.null(m$na.action) && m$na.action != "na.omit") {
    warning("Only 'na.omit' is implemented for gee\ncontinuing with 'na.action=na.omit'")
    m$na.action <- as.name("na.omit")
  }
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  y <- as.matrix(model.extract(m, "response"))
  x <- model.matrix(Terms, m, contrasts)
  xx <- model.matrix(Terms, m, contrasts)
  xx <- as.data.frame(xx)
  id <- model.extract(m, id)
  if (dim(y)[2]==2) {
    N <- as.vector(y %*% c(1,1))
    y <- y[,1]
  }
  else {
    if (dim(y)[2]>2)
      stop("Only binomial response matrices (2 columns)")
  }

  if(is.null(id)) {
    stop("Id variable not found")
  }

  xnames <- dimnames(xx)[[2]]
  if(is.null(xnames)) {
    xnames <- paste("x", 1:p, sep = "")
    dimnames(xx) <- list(NULL, xnames)
  }

  k <- pp <- ncol(xx)
  ##########
  tx <- t(xx)
  tol<-2
  ####
  beta <- rep(0, k)
  while(tol > .0001){
    z <- as.matrix(xx)%*%beta
    pi <- as.vector(1/(1 + exp( - z )))
    W <- diag(pi*(1-pi))
    W05 <- W^.5
    XW2 <- (tx %*% W05)    #### X' (W ^ 1/2)
    XWXi <- solve(tx %*% W %*% as.matrix(xx))
    Hmat <- t(XW2) %*% XWXi %*% XW2
    UU <- y - pi + diag(diag(Hmat)) %*% (.5-pi)
    U.star <- tx %*% UU
    delta <- as.vector(XWXi %*% U.star)
    beta <- beta + delta
    tol <- max(abs(delta))
  }

  ############
  xx <- split(xx, id)
  y <- split(y, id)
  nc <- length(y)
  bet <- e <-  V <- DD <- part <- Z <- S1 <- S2 <- S3 <- FF <- sec.part <- list()
  p2 <- NULL
  del <- 10
  counter=0
  while((del > .0001) & (counter<iterlim)){
    counter = counter+1
    z <- lapply(xx, function(a) as.matrix(a) %*% beta)
    a <- lapply(z, function(a) exp(a))
    mu <- lapply(a, function(a) a/(1+a))
    var.mu <- lapply(mu, function(a) a * (1-a))
    W <- lapply(var.mu, function(x) diag(as.vector(x), nrow = length(x), ncol = length(x)))
    for(i in 1:nc){
      e[[i]] <- (y[[i]]-mu[[i]])/sqrt(var.mu[[i]])
    }

    if(corstr == "exchangeable"){
      cl.size <- sapply(e, function(x) length(x))
      frame <- lapply(e, function(x) data.frame(a=x, b=x))
      matr <- lapply(frame, function(x) with(x, sapply(a, function(x) x*b)))
      matr.sum <- lapply(matr, function(x){
        p1 <- as.matrix(x)
        ni <- nrow(p1)
        s <- sum(p1)-sum(diag(p1))
        mean.s <- ifelse(ni == 1, 0, s/(ni*(ni-1)))
        return(mean.s)

      })
      alpha <- mean(unlist(matr.sum))


      R <- lapply(cl.size, function(x) xch(x, alpha))
    }
    if(corstr == "ar1"){
      e1 <- lapply(e, function(x) x[-length(x)])
      e2 <- lapply(e, function(x) x[-1])
      ni <- sapply(e, function(x) length(x))
      temp <- list()
      for(i in 1:nc){
        temp[[i]] <- ifelse(ni[i] == 1, 0, (1/(ni[i]-1)) * sum(e1[[i]] * e2[[i]]))
      }
      alpha <- Reduce("+", temp)/nc

      R <- lapply(ni, function(x) ar1(x, alpha))
    }
    if(corstr == "unstr"){
      ni <- sapply(e, function(x) length(x))
      frame <- lapply(e, function(x) data.frame(a=x, b=x))
      matr <- lapply(frame, function(x) with(x, sapply(a, function(x) x*b)))
      new.mat <- sapply(ni, function(x) max(ni)-x)
      matr.s <- list()
      for(i in 1:nc){
        matr.s[[i]] <- rbind(cbind(as.matrix(matr[[i]]), matrix(0, nrow =  ni[i], ncol = new.mat[i])),
                             matrix(0, nrow =  new.mat[i], ncol = max(ni)))
      }
      alpha <- Reduce("+", matr.s)/nc
      R <- lapply(ni, function(x) alpha[1:x, 1:x])
    }
    if(corstr=="independence"){
      ni <- sapply(e, function(x) length(x))
      R <- lapply(ni, function(x) diag(x))
    }

    if(est_dispersion)
      phi <- mean(unlist(lapply(e, function(x) mean(x^2))))
    else
      phi <- 1
    for(i in 1:nc){
      V[[i]] <- phi * (W[[i]])^(1/2) %*% R[[i]] %*% (W[[i]])^(1/2)
    }
    for(i in 1:nc){
      DD[[i]] <- t(as.matrix(xx[[i]])) %*% W[[i]]

    }
    for(i in 1:nc){
      xi <- as.matrix(xx[[i]]); txi <- t(xi); W12 <- W[[i]]^(1/2)
      part[[i]] <- txi %*% W12 %*% ginv(R[[i]]) %*% W12 %*% xi
    }
    I <- Reduce("+", part)/phi
    Q <- lapply(mu, function(x) diag(0.5-as.vector(x), nrow = length(x), ncol = length(x)))
    for(i in 1:nc){
      Z[[i]] <- list()
      xi <- as.matrix(xx[[i]])
      pi <- ncol(xi)
      for(j in 1:pi){
        temp <- as.vector(xi[,j])
        Z[[i]][[j]] <- diag(temp, nrow = length(temp), ncol = length(temp))
      }
    }

    for(i in 1:nc){
      W12 <- W[[i]]^(1/2)
      xi <- as.matrix(xx[[i]])
      txi <- t(xi)
      Qi <- Q[[i]]
      S1[[i]] <- list()
      for(j in 1:pp){
        S1[[i]][[j]] <- txi %*% W12 %*% ginv(R[[i]]) %*% W12 %*% Qi %*% Z[[i]][[j]] %*% xi
      }

    }
    for(j in 1:pp){
      FF[[j]] <- matrix(0, pp, pp)
      for(i in 1:nc){
        FF[[j]] <- FF[[j]] + S1[[i]][[j]]
      }
      FF[[j]] <- 2 * FF[[j]]/phi
    }
    for(i in 1:pp){
      p2[i] <-  sum(diag(ginv(I) %*% FF[[i]]))
    }
    for(i in 1:nc){
      xi <- as.matrix(xx[[i]]);  txi <- t(xi);  Wi <- W[[i]];  yi <- y[[i]]; mui <- mu[[i]]
      sec.part[[i]] <- txi %*% Wi %*% ginv(V[[i]]) %*% (yi - mui)

    }
    Ustar <- Reduce("+", sec.part) + (0.5 * p2)
    tt <- ginv(I) %*% Ustar
    del <- max(abs(tt))
    beta <- beta + tt
  }

  ##########################################
  temp <- e <- V <- DD <- sec.part <- tempS <- list()
  ni <- NULL
  z <- lapply(xx, function(a) as.matrix(a) %*% beta)
  mu <- lapply(z, function(a) exp(a)/(1+exp(a)))
  var.mu <- lapply(mu, function(a) a * (1-a))
  W <- lapply(var.mu, function(x) diag(as.vector(x), nrow = length(x), ncol = length(x)))
  for(i in 1:nc){
    e[[i]] <- (y[[i]]-mu[[i]])/sqrt(var.mu[[i]])
  }
  if(corstr == "exchangeable"){
    cl.size <- sapply(e, function(x) length(x))
    frame <- lapply(e, function(x) data.frame(a=x, b=x))
    ni <-  cl.size
    matr <- lapply(frame, function(x) with(x, sapply(a, function(x) x*b)))
    matr.sum <- lapply(matr, function(x){
      p1 <- as.matrix(x)
      nni <- nrow(p1)
      s <- sum(p1)-sum(diag(p1))
      mean.s <- ifelse(nni == 1, 0, s/(nni*(nni-1)))
      return(mean.s)
    })
    alpha <- mean(unlist(matr.sum))
    R <- lapply(cl.size, function(x) xch(x, alpha))
  }
  if(corstr == "ar1"){
    e1 <- lapply(e, function(x) x[-length(x)])
    e2 <- lapply(e, function(x) x[-1])
    ni <- sapply(e, function(x) length(x))
    temp <- list()
    for(i in 1:nc){
      temp[[i]] <- ifelse(ni[i] == 1, 0, (1/(ni[i]-1)) * sum(e1[[i]] * e2[[i]]))
    }
    alpha <- Reduce("+", temp)/nc
    R <- lapply(ni, function(x) ar1(x, alpha))
  }
  if(corstr == "unstr"){
    ni <- sapply(e, function(x) length(x))
    frame <- lapply(e, function(x) data.frame(a=x, b=x))
    matr <- lapply(frame, function(x) with(x, sapply(a, function(x) x*b)))
    new.mat <- sapply(ni, function(x) max(ni)-x)
    matr.s <- list()
    for(i in 1:nc){
      matr.s[[i]] <- rbind(cbind(as.matrix(matr[[i]]), matrix(0, nrow =  ni[i], ncol = new.mat[i])),
                           matrix(0, nrow =  new.mat[i], ncol = max(ni)))
    }
    alpha <- Reduce("+", matr.s)/nc
    R <- lapply(ni, function(x) alpha[1:x, 1:x, drop = FALSE])
  }
  if(corstr=="independence"){
    ni <- sapply(e, function(x) length(x))
    R <- lapply(ni, function(x) diag(x))
  }

  for(i in 1:nc){
    V[[i]] <- phi * (W[[i]])^(1/2) %*% R[[i]] %*% (W[[i]])^(1/2)
  }
  for(i in 1:nc){
    DD[[i]] <- t(as.matrix(xx[[i]])) %*% W[[i]]

  }
  for(i in 1:nc){
    temp[[i]] <- (DD[[i]]) %*% ginv(V[[i]]) %*% t(DD[[i]])
  }
  Vm <- ginv(Reduce("+",temp))
  for(i in 1:nc){
    xi <- as.matrix(xx[[i]]); txi <- t(xi); Wi <- W[[i]]; yi <- y[[i]]
    mui <- mu[[i]]
    sec.part[[i]] <- txi %*% Wi %*% ginv(V[[i]]) %*% (yi - mui)

  }
  d.bar <- Reduce("+", sec.part)/nc
  for(i in 1:nc){
    t1 <- (sec.part[[i]] - d.bar)
    tempS[[i]] <- t1 %*% t(t1)
  }
  N <- sum(ni)
  B <- Reduce("+", tempS) * (nc/(nc-1)) * ((N-1)/(N-pp))
  fi <- max(1, sum(diag(Vm %*% B))/nrow(B))
  deln <- min(0.5, (nrow(B)/(nc - nrow(B))))
  VsM <- (Vm %*% B %*% Vm) + (fi * deln * Vm)
  diag(VsM) <- abs(diag(VsM))
  VsMod=sqrt(diag(VsM))
  dimen <- sapply(R, function(x) dim(x)[2])
  ##########################################
  t.swM=round(beta/VsMod, 4)
  p.value=round(pf(t.swM^2, 1, N-pp, lower.tail = F), 4)
  est.swM <- data.frame(coefficients=beta, std.err=VsMod, Wald=(t.swM)^2, p.val=p.value)
  row.names(est.swM) <- xnames

  fit <- list()
  #attr(fit, "class") <- c("gee")
  fit$call <- call
  fit$formula <- as.vector(attr(Terms, "formula"))
  fit$coefficient <- est.swM
  fit$correlation <- R[dimen == max(dimen)][[1]]
  fit$V = VsM
  fit

}















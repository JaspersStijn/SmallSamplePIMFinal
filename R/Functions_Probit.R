## Functions for small sample corrections
small_sample_GEE_WL = function (formula, id, family = binomial, data, corstr = "independence",link="probit",wl_orig=TRUE)
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
    for (i in 1:size) {
      y <- as.matrix(data$response[data$id == unique(data$id)[i]])
      covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                    mat$subj == unique(data$id)[i]))
      var_i = var[1:cluster$n[i], 1:cluster$n[i]]
      D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
      Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
      Vi = diag(1e-5,nrow(Vi))+Vi
      xx <- t(D) %*% solve(Vi) %*% D
      step01 <- step01 + xx
    }
    step01 = step01 + diag(1e-5,nrow(step01))

    step <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
    for (i in 1:size) {
      y <- as.matrix(data$response[data$id == unique(data$id)[i]])
      covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                    mat$subj == unique(data$id)[i]))
      var_i = var[1:cluster$n[i], 1:cluster$n[i]]

      B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
      diag(B) <- 1/sqrt(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est))+1e-5)
      D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
      Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
      Vi = diag(1e-5,nrow(Vi))+Vi

      resid <- B %*% solve(cormax.ind(cluster$n[i]) - D %*%
                             solve(step01) %*% t(D) %*% solve(Vi)) %*% (y - pnorm(covariate %*% beta_est))
      step <- step + resid %*% t(resid)
    }
    if(wl_orig){unstr <- step/size}
    if(!wl_orig){unstr <- step/(size-length(beta_est))}
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


      B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
      diag(B) <- pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est))
      D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
      Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
        var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
      Vi = diag(1e-5,nrow(Vi))+Vi

      xy <- t(D) %*% solve(Vi) %*% sqrt(B) %*% unstr %*%
        sqrt(B) %*% solve(Vi) %*% D
      xx <- t(D) %*% solve(Vi) %*% D
      step12 <- step12 + xy
      step11 <- step11 + xx
      step13 <- step13 + vec(xy)
      p[, i] <- vec(xy)

    }
    step11 = step11+diag(1e-5,nrow(step11))
    for (i in 1:size) {
      dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
      step14 <- step14 + dif
    }
    cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
    cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
      step14 %*% kronecker(solve(step11), solve(step11))
  return(list(beta = beta_est ,cov.beta = diag(cov.beta), cov.var = cov.var))
}

small_sample_GEE_GST = function (formula, id, family = binomial, data, corstr = "independence",link="probit",method="Gosh")
{
  if(link=="probit"){
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

    gee.fit <-  geese(formula, data = data, family = binomial(link="probit"), corstr = corstr, id = id)

    beta_est <- gee.fit$beta
    alpha <- gee.fit$alpha
    len <- length(beta_est)
    len_vec <- len^2
    data$id <- data[,id]
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
      print(corstr)
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
        print(family)
        stop("'family' not recognized")
      }
    }
    m <- model.frame(formula, data)
    mat <- as.data.frame(model.matrix(formula, m))
    mat$subj <- rep(unique(data$id), cluster$n)

    cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)
    step <- matrix(0, nrow = cluster$n[1], ncol = cluster$n[1])
    for (i in 1:size) {
      y <- as.matrix(data$response[data$id == unique(data$id)[i]])
      covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                    mat$subj == unique(data$id)[i]))
      if (family == "binomial") {
        resid <- (y - pnorm(covariate %*% beta_est)) %*% t(y - pnorm(covariate %*% beta_est))
        B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
        diag(B) <- 1/sqrt(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est))+1e-5)
        step <- step + B %*% resid %*% B
      }
      else {
        stop("'family' not recognized")
      }
    }
    step <- step+ diag(1e-5,nrow(step))
    if(method=="Gosh"){unstr <- step/(size - len)}
    if(method=="Pan"){unstr <- step/(size)}

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
        diag(A) <- pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est))

        D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
        Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
          var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
        Vi = diag(1e-5,nrow(Vi))+Vi

        xy <- t(D) %*% solve(Vi) %*% sqrt(A) %*% unstr %*%
          sqrt(A) %*% solve(Vi) %*% D
        xx <- t(D) %*% solve(Vi) %*% D
        step12 <- step12 + xy
        step11 <- step11 + xx
        step13 <- step13 + vec(xy)
        p[, i] <- vec(xy)
      }
    }
    step11=step11+ diag(1e-5,nrow(step11))
    step12=step12+ diag(1e-5,nrow(step12))
    step14=step14+ diag(1e-5,nrow(step14))

    for (i in 1:size) {
      dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
      step14 <- step14 + dif
    }
    cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
    cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
      step14 %*% kronecker(solve(step11), solve(step11))
  }
  return(list(beta = beta_est ,cov.beta = diag(cov.beta), cov.var = cov.var))
}


small_sample_GEE_MBN = function (formula, id, family = binomial, data, corstr = "independence",link="probit", d = 2, r = 1)
{
  if(link=="probit"){
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
    gee.fit <-    geese(formula, data = data, family = binomial(link="probit"), corstr = corstr, id = id)
    beta_est <- gee.fit$beta
    alpha <- gee.fit$alpha
    len <- length(beta_est)
    len_vec <- len^2
    data$id <- data[,id]
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
      print(corstr)
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
        print(family)
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
      var_i = var[1:ncluster, 1:ncluster]
      if (family == "binomial") {
        D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
        Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
          var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
        Vi = diag(1e-5,nrow(Vi))+Vi
        xx <- t(D) %*% solve(Vi) %*% D
        step11 <- step11 + xx
      }
      else {
        stop("'family' not recognized")
      }
    }
    step11 = step11+diag(1e-5,nrow(step11))
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
      var_i = var[1:ncluster, 1:ncluster]

      if (family == "binomial") {
        D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
        Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
          var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
        Vi = diag(1e-5,nrow(Vi))+Vi

        xy <- t(D) %*% solve(Vi) %*% (y - pnorm(covariate %*% beta_est))
        step00 <- step00 + xy %*% t(xy)
      }
    }
    step00 = step00+diag(1e-5,nrow(step00))

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
      var_i = var[1:ncluster, 1:ncluster]
      if (family == "binomial") {
        D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
        Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
          var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
        Vi = diag(1e-5,nrow(Vi))+Vi
        xy <- t(D) %*% solve(Vi) %*% (k * (y - pnorm(covariate %*% beta_est)) %*%
                                        t(y - pnorm(covariate %*% beta_est)) + delta * xi * Vi) %*% solve(Vi) %*%
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
    step12 = step12+diag(1e-5,nrow(step12))
    step14 = step14+diag(1e-5,nrow(step14))

    cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
    cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
      step14 %*% kronecker(solve(step11), solve(step11))
    return(list(beta = beta_est ,cov.beta = diag(cov.beta), cov.var = cov.var))
  }}

small_sample_GEE_MD = function (formula, id, family = binomial, data, corstr = "independence",link="probit")
{
  if(link=="probit"){
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

    gee.fit <-  geese(formula, data = data, family = binomial(link="probit"), corstr = corstr, id = id)
    beta_est <- gee.fit$beta
    alpha <- gee.fit$alpha
    len <- length(beta_est)
    len_vec <- len^2
    data$id <- data[,id]
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
      print(corstr)
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
        print(family)
        stop("'family' not recognized")
      }
    }
    m <- model.frame(formula, data)
    mat <- as.data.frame(model.matrix(formula, m))
    mat$subj <- rep(unique(data$id), cluster$n)

    cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)
    step11 <- matrix(0, nrow = len, ncol = len)
    for (i in 1:size) {
      y <- as.matrix(data$response[data$id == unique(data$id)[i]])
      covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                    mat$subj == unique(data$id)[i]))
      var_i = var[1:cluster$n[i], 1:cluster$n[i]]
      if (family == "binomial") {
        D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
        Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
          var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
        Vi = diag(1e-5,nrow(Vi))+Vi
        xx <- t(D) %*% solve(Vi) %*% D
        step11 <- step11 + xx
      }
      else {
        stop("'family' not recognized")
      }
    }

    step11 = step11+diag(1e-5,nrow(step11))

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
        D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
        Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
          var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
        Vi = diag(1e-5,nrow(Vi))+Vi
        xy <- t(D) %*% solve(Vi) %*% solve(cormax.ind(cluster$n[i]) -
                                             D %*% solve(step11) %*% t(D) %*% solve(Vi)) %*%
          (y - pnorm(covariate %*% beta_est))
        step12 <- step12 + xy %*% t(xy)
        step13 <- step13 + vec(xy %*% t(xy))
        p[, i] <- vec(xy %*% t(xy))
      }
    }

    for (i in 1:size) {
      dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
      step14 <- step14 + dif
    }
    step12 = step12+diag(1e-5,nrow(step12))
    step14 = step14+diag(1e-5,nrow(step14))
    cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
    cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
      step14 %*% kronecker(solve(step11), solve(step11))
  }
  return(list(beta = beta_est ,cov.beta = diag(cov.beta), cov.var = cov.var))
}

small_sample_GEE_KC = function (formula, id, family = binomial, data, corstr = "independence",link="probit")
{
  if(link=="probit"){
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

    gee.fit <-  geese(formula, data = data, family = binomial(link="probit"), corstr = corstr, id = id)
    beta_est <- gee.fit$beta
    alpha <- gee.fit$alpha
    len <- length(beta_est)
    len_vec <- len^2
    data$id <- data[,id]
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
      print(corstr)
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
        print(family)
        stop("'family' not recognized")
      }
    }
    m <- model.frame(formula, data)
    mat <- as.data.frame(model.matrix(formula, m))
    mat$subj <- rep(unique(data$id), cluster$n)

    cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)
    step11 <- matrix(0, nrow = len, ncol = len)
    for (i in 1:size) {
      y <- as.matrix(data$response[data$id == unique(data$id)[i]])
      covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                    mat$subj == unique(data$id)[i]))
      var_i = var[1:cluster$n[i], 1:cluster$n[i]]
      if (family == "binomial") {
        D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
        Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
          var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
        Vi = diag(1e-5,nrow(Vi))+Vi
        xx <- t(D) %*% solve(Vi) %*% D
        step11 <- step11 + xx
      }
      else {
        stop("'family' not recognized")
      }
    }
    step11 = step11+diag(1e-5,nrow(step11))

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
        D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
        Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
          var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
        Vi = diag(1e-5,nrow(Vi))+Vi

        xy <- t(D) %*% solve(Vi) %*% mat.sqrt.inv(cormax.ind(cluster$n[i]) -
                                                    D %*% solve(step11) %*% t(D) %*% solve(Vi)) %*%
          (y - pnorm(covariate %*% beta_est))
        step12 <- step12 + xy %*% t(xy)
        step13 <- step13 + vec(Re(xy %*% t(xy)))
        p[, i] <- vec(Re(xy %*% t(xy)))
      }
    }

    for (i in 1:size) {
      dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
      step14 <- step14 + dif
    }
    step12 = step12+diag(1e-5,nrow(step12))
    step14 = step14+diag(1e-5,nrow(step14))
    cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
    cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
      step14 %*% kronecker(solve(step11), solve(step11))
  }
  return(list(beta = beta_est ,cov.beta = diag(Re(cov.beta)), cov.var = cov.var))
}


small_sample_GEE_FG = function(formula, id, family = binomial, data, corstr = "independence",
                               b = 0.75,link="probit")
{
  if(link=="probit"){
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

    gee.fit <-  geese(formula, data = data, family = binomial(link="probit"), corstr = corstr, id = id)
    beta_est <- gee.fit$beta
    alpha <- gee.fit$alpha
    len <- length(beta_est)
    len_vec <- len^2
    data$id <- data[,id]
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
      print(corstr)
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
        print(family)
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
      var_i = var[1:ncluster, 1:ncluster]
      if (family == "binomial") {
        D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
        Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
          var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
        Vi = diag(1e-5,nrow(Vi))+Vi
        xx <- t(D) %*% solve(Vi) %*% D
        step11 <- step11 + xx
      }
      else {
        stop("'family' not recognized")
      }
    }
    step11 = step11+diag(1e-5,nrow(step11))
    step12 <- matrix(0, nrow = len, ncol = len)
    step13 <- matrix(0, nrow = len_vec, ncol = 1)
    step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
    p <- matrix(0, nrow = len_vec, ncol = size)
    for (i in 1:size) {
      y <- as.matrix(data$response[data$id == unique(data$id)[i]])
      covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                    mat$subj == unique(data$id)[i]))
      ncluster = cluster$n[i]
      var_i = var[1:ncluster, 1:ncluster]
      if (family == "binomial") {
        D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
        Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
          var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
        Vi = diag(1e-5,nrow(Vi))+Vi

        xx <- t(D) %*% solve(Vi) %*% D
        Qi <- xx %*% solve(step11)
        Ai <- diag((1 - pmin(b, diag(Qi)))^(-0.5))
        xy <- Ai %*% t(D) %*% solve(Vi) %*% (y - pnorm(covariate %*% beta_est))
        step12 <- step12 + xy %*% t(xy)
        step13 <- step13 + vec(xy %*% t(xy))
        p[, i] <- vec(xy %*% t(xy))
      }
    }
    for (i in 1:size) {
      dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
      step14 <- step14 + dif
    }
    step12 = step12+diag(1e-5,nrow(step12))
    step14 = step14+diag(1e-5,nrow(step14))
    cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
    cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
      step14 %*% kronecker(solve(step11), solve(step11))
  }
  return(list(beta = beta_est ,cov.beta = diag(cov.beta), cov.var = cov.var))
}

small_sample_GEE_MK = function (formula, id, family = binomial, data, corstr = "independence",link="probit")
{
  if(link=="probit"){
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

    gee.fit <-  geese(formula, data = data, family = binomial(link="probit"), corstr = corstr, id = id)

    beta_est <- gee.fit$beta
    alpha <- gee.fit$alpha
    len <- length(beta_est)
    len_vec <- len^2
    data$id <- data[,id]
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
      print(corstr)
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
        print(family)
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
        D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
        Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
          var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
        Vi = diag(1e-5,nrow(Vi))+Vi
        xx <- t(D) %*% solve(Vi) %*% D
        step11 <- step11 + xx
      }
      else {
        stop("'family' not recognized")
      }
    }
    step11 = step11+diag(1e-5,nrow(step11))

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
        D <- mat.prod(covariate, dnorm(covariate %*% beta_est))+1e-5
        Vi <- diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i]) %*%
          var_i %*% diag(sqrt(c(pnorm(covariate %*% beta_est)*(1-pnorm(covariate %*% beta_est)))), cluster$n[i])
        Vi = diag(1e-5,nrow(Vi))+Vi

        xy <- t(D) %*% solve(Vi) %*% (y - pnorm(covariate %*% beta_est))
        step12 <- step12 + xy %*% t(xy)
        step13 <- step13 + vec(xy %*% t(xy))
        p[, i] <- vec(xy %*% t(xy))
      }
    }
    for (i in 1:size) {
      dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
      step14 <- step14 + dif
    }
    step12 = step12+diag(1e-5,nrow(step12))
    step14 = step14+diag(1e-5,nrow(step14))
    scale <- size/(size - len)
    cov.beta <- scale * solve(step11) %*% (step12) %*% solve(step11)
    cov.var <- scale^2 * size/(size - 1) * kronecker(solve(step11),
                                                     solve(step11)) %*% step14 %*% kronecker(solve(step11),
                                                                                             solve(step11))
  }
  return(list(beta = beta_est ,cov.beta = diag(cov.beta), cov.var = cov.var))
}

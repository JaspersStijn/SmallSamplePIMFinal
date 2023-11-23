require(geesmv)
require(matrixcalc)
require(MASS)
source("./R/Functions_Logit.R")
source("./R/Functions_Probit.R")


GEE_MH_fit = function(data, response, treatment, control, correction = "MBN",link="logit"){

  ## Create the cluster indicators

  use_vals = unique(data[,treatment])
  if(length(use_vals) !=2) {print("error: we can only compare between two groups")}
  id.nonfac <- which(data[,treatment] == use_vals[1])
  id.fac <- which(data[,treatment] == use_vals[2])

  compare <- expand.grid(id.nonfac, id.fac)
  y1 <- data[id.nonfac,response]
  y2 <- data[id.fac,response]
  pseudo_y <- as.vector(t(outer(y1, y2, `<`)+0.5*(outer(y1, y2, `==`)))*1)
  pseudo_dat = data.frame("y"=pseudo_y)

  for(i in 1:length(control)){
    x_g1 <- data[id.nonfac,control[i]]
    x_g2 <- data[id.fac,control[i]]
    pseudo_x <- as.vector(t(outer(x_g1, x_g2, `-`)))
    pseudo_dat[,control[i]] = -pseudo_x
  }

  pseudo_dat$C1 = rep((1:length(id.nonfac)),rep(length(id.fac),length(id.nonfac)))
  pseudo_dat$C2 = rep((1:(length(id.fac))),length(id.nonfac))
  pseudo_dat$C1C2 = 1:nrow(pseudo_dat)

  form <- as.formula(paste("y~", paste(control, collapse="+")))

  if(link=="logit"){
  ## Run the three models

  if(correction=="MBN"){
    mod1 =  GEE.var.mbn_new(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

    mod3 = GEE.var.mbn_new(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

    pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
    mod2 = GEE.var.mbn_new(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

    V = mod1$cov.beta+
      mod2$cov.beta-
      mod3$cov.beta
  }

  if(correction=="Pan"){
    mod1 =  GEE.var.pan_new(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

    mod3 = GEE.var.pan_new(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

    pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
    mod2 = GEE.var.pan_new(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

    V = mod1$cov.beta+
      mod2$cov.beta-
      mod3$cov.beta
  }


  if(correction=="GST"){
    mod1 =  GEE.var.gst_new(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

    mod3 = GEE.var.gst_new(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

    pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
    mod2 = GEE.var.gst_new(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

    V = mod1$cov.beta+
      mod2$cov.beta-
      mod3$cov.beta
  }

  if(correction=="KC"){
    mod1 =  GEE.var.kc_new(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

    mod3 = GEE.var.kc_new(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

    pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
    mod2 = GEE.var.kc_new(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

    V = mod1$cov.beta+
      mod2$cov.beta-
      mod3$cov.beta
  }

  if(correction=="MD"){
    mod1 =  GEE.var.md_new(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

    mod3 = GEE.var.md_new(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

    pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
    mod2 = GEE.var.md_new(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

    V = mod1$cov.beta+
      mod2$cov.beta-
      mod3$cov.beta
  }


  if(correction=="FG"){
    mod1 =  GEE.var.fg_new(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

    mod3 = GEE.var.fg_new(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

    pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
    mod2 = GEE.var.fg_new(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

    V = mod1$cov.beta+
      mod2$cov.beta-
      mod3$cov.beta
  }

  if(correction=="MK"){
    mod1 =  GEE.var.mk_new(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

    mod3 = GEE.var.mk_new(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

    pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
    mod2 = GEE.var.mk_new(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

    V = mod1$cov.beta+
      mod2$cov.beta-
      mod3$cov.beta
  }

  if(correction=="WL"){
    mod1 =  GEE.var.wl_new(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

    mod3 = GEE.var.wl_new(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

    pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
    mod2 = GEE.var.wl_new(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

    V = mod1$cov.beta+
      mod2$cov.beta-
      mod3$cov.beta
  }

  if(correction=="Firth"){
    mod1 = test_firth(formula=form,data=pseudo_dat, id = C1,corstr = "independence")
    mod3 = test_firth(formula=form,data=pseudo_dat, id = C1C2,corstr = "independence")

    pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
    mod2 = test_firth(formula=form,data=pseudo_dat, id = C2,corstr = "independence")

    V = mod1$coefficient[,2]^2+
      mod2$coefficient[,2]^2-
      mod3$coefficient[,2]^2
    mod1$beta = mod1$coefficient[,1]
  }

  df = nrow(pseudo_dat) - length(V) - 1

  test_stat = mod1$beta[1]/sqrt(V[1])
  pval = 2*pt(abs(test_stat),lower.tail = FALSE,df=df)

  test_stat_all = mod1$beta/sqrt(V)
  pval_all = 2*pt(abs(test_stat_all),lower.tail = FALSE,df=df)


  ## Structured output

  fit <- list()
  fit$title <- "PIM with small sample correction"
  fit$model <- list()
  fit$model$link <- "logit"
  fit$formula <- form

  fit$coefficients <- as.vector(mod1$beta)
  fit$nas <- is.na(fit$coefficients)
  names(fit$coefficients) <- c("Intercept",control)
  fit$robust.variance <- V
  fit$test = test_stat_all
  fit$pval = pval_all
  new("pim_corrected", fit)
  return(list("beta_est" = mod1$beta[1],"beta_se" = sqrt(V[1]),"lower" = mod1$beta[1]-qt(0.975,df)*sqrt(V[1]),"upper" = mod1$beta[1]+qt(0.975,df)*sqrt(V[1]),pval=pval))
  }


  if(link=="probit"){
    ## Run the three models

    if(correction=="MBN"){
      mod1 =  small_sample_GEE_MBN(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

      mod3 = small_sample_GEE_MBN(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

      pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
      mod2 = small_sample_GEE_MBN(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

      V = mod1$cov.beta+
        mod2$cov.beta-
        mod3$cov.beta
    }

    if(correction=="Pan"){
      mod1 =  small_sample_GEE_GST(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence",method="Pan")

      mod3 = small_sample_GEE_GST(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence",method="Pan")

      pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
      mod2 = small_sample_GEE_GST(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence",method="Pan")

      V = mod1$cov.beta+
        mod2$cov.beta-
        mod3$cov.beta
    }


    if(correction=="GST"){
      mod1 =  small_sample_GEE_GST(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

      mod3 = small_sample_GEE_GST(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

      pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
      mod2 = small_sample_GEE_GST(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

      V = mod1$cov.beta+
        mod2$cov.beta-
        mod3$cov.beta
    }

    if(correction=="KC"){
      mod1 =  small_sample_GEE_KC(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

      mod3 = small_sample_GEE_KC(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

      pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
      mod2 = small_sample_GEE_KC(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

      V = mod1$cov.beta+
        mod2$cov.beta-
        mod3$cov.beta
    }

    if(correction=="MD"){
      mod1 =  small_sample_GEE_MD(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

      mod3 = small_sample_GEE_MD(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

      pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
      mod2 = small_sample_GEE_MD(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

      V = mod1$cov.beta+
        mod2$cov.beta-
        mod3$cov.beta
    }


    if(correction=="FG"){
      mod1 =  small_sample_GEE_FG(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

      mod3 = small_sample_GEE_FG(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

      pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
      mod2 = small_sample_GEE_FG(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

      V = mod1$cov.beta+
        mod2$cov.beta-
        mod3$cov.beta
    }

    if(correction=="MK"){
      mod1 =  small_sample_GEE_MK(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

      mod3 = small_sample_GEE_MK(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

      pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
      mod2 = small_sample_GEE_MK(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

      V = mod1$cov.beta+
        mod2$cov.beta-
        mod3$cov.beta
    }

    if(correction=="WL"){
      mod1 =  small_sample_GEE_WL(formula=form,id="C1",family=binomial,data=pseudo_dat,corstr="independence")

      mod3 = small_sample_GEE_WL(formula=form,id="C1C2",family=binomial,data=pseudo_dat,corstr="independence")

      pseudo_dat = pseudo_dat[order(pseudo_dat$C2),]
      mod2 = small_sample_GEE_WL(formula=form,id="C2",family=binomial,data=pseudo_dat,corstr="independence")

      V = mod1$cov.beta+
        mod2$cov.beta-
        mod3$cov.beta
    }

    if(correction=="Firth"){
      stop("We cannot apply Firth with probit link, please use option 'link=logit'")
      }

    df = nrow(pseudo_dat) - length(V) - 1

    test_stat = mod1$beta[1]/sqrt(V[1])
    pval = 2*pt(abs(test_stat),lower.tail = FALSE,df=df)

    test_stat_all = mod1$beta/sqrt(V)
    pval_all = 2*pt(abs(test_stat_all),lower.tail = FALSE,df=df)


    ## Structured output

    fit <- list()
    fit$title <- "PIM with small sample correction"
    fit$model <- list()
    fit$model$link <- "probit"
    fit$formula <- form

    fit$coefficients <- as.vector(mod1$beta)
    fit$nas <- is.na(fit$coefficients)
    names(fit$coefficients) <- c("Intercept",control)
    fit$robust.variance <- V
    fit$test = test_stat_all
    fit$pval = pval_all
    new("pim_corrected", fit)
    return(list("beta_est" = mod1$beta[1],"beta_se" = sqrt(V[1]),"lower" = mod1$beta[1]-qt(0.975,df)*sqrt(V[1]),"upper" = mod1$beta[1]+qt(0.975,df)*sqrt(V[1]),pval=pval))
    }

}


setClass( "pim_corrected", representation("list"))

# show method (here's how the output would be printed
# you can format to whatever you want... to show and how to show
setMethod("show", "pim_corrected", function(object) {
  out_f = data.frame(cbind(object$coefficients,sqrt(object$robust.variance),object$test,object$pval))
  names(out_f) = c("Estimate","SE","Test statistic","Pvalue")
  cat("------------------------\n")
  cat("Model Output: \n")
  cat("------------------------\n")
  print(out_f)
})

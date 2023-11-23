#------------------------------------------------------------------------
# Generate data for the simulation study in Jaspers et al. (2023)
#------------------------------------------------------------------------

# Scenario defines one of the two scenarios in the paper:
#      - Scenario 1 corresponds to the normal linear model
#      - Scenario 2 corresponds to the exponential model

# samplesize corresponds to the total sample size of the dataset

# alpha1 corresponds to the treatment effect (0, 0.5, 1 in the paper)

# p corresponds to the number of extra covariates (1 or 3 in the paper)

#------------------------------------------------------------------------

generate_data = function(scenario = 1, samplesize=20, alpha1=1, p=1){

  if(scenario == 1){

    x1 = rep(c(0,1),each = samplesize/2)
    x2 = runif(samplesize,0.1,2)

    if(p==3){
      x3 = runif(samplesize,0.1,2)
      x4 = runif(samplesize,0.1,2)
    }

    if(p==5){
      x3 = runif(samplesize,0.1,2)
      x4 = runif(samplesize,0.1,2)
      x5 = runif(samplesize,0.1,2)
      x6 = runif(samplesize,0.1,2)
    }

    if(p==1){
      y = rnorm(samplesize,alpha1*x1+1*x2,1)
      dat = data.frame(cbind(y,x1,x2))
    }

    if(p==3){
      y = rnorm(samplesize,alpha1*x1+1*x2+1*x3+1*x4,1)
      dat = data.frame(cbind(y,x1,x2,x3,x4))
    }

    if(p==5){
      y = rnorm(samplesize,alpha1*x1+1*x2+1*x3+1*x4+1*x5+1*x6,1)
      dat = data.frame(cbind(y,x1,x2,x3,x4,x5,x6))
    }


    id.nonfac <- which(dat$x1 == 0)
    id.fac <- which(dat$x1 == 1)
    compare <- expand.grid(id.nonfac, id.fac)
    y1 <- dat[id.nonfac,"y"]
    y2 <- dat[id.fac,"y"]
    pseudo_y <- as.vector(t(outer(y1, y2, `<`))*1)
    pseudo_dat = data.frame("y"=pseudo_y)

    for(i in 1:p){
      x_g1 <- dat[id.nonfac,paste0("x",i+1)]
      x_g2 <- dat[id.fac,paste0("x",i+1)]
      pseudo_x <- as.vector(t(outer(x_g1, x_g2, `-`)))
      pseudo_dat[,paste0("x",i+1)] = -pseudo_x
    }

    pseudo_dat$C1 = rep((1:(samplesize/2)),each=samplesize/2)
    pseudo_dat$C2 = rep((1:(samplesize/2)),samplesize/2)
    pseudo_dat$C1C2 = 1:nrow(pseudo_dat)

    pseudo_form <- as.formula(paste("y~", paste(rep("x",p), (1:p)+1, sep="", collapse="+")))

    return(list(data=dat,form=as.formula(paste("y~", paste(rep("x",p), (1:(p+1)), sep="", collapse="+"))),pseudo_dat=pseudo_dat,pseudo_form=pseudo_form,comp=compare))
  }

  if(scenario == 2){
    x1 = rep(c(0,1),each = samplesize/2)
    x2 = runif(samplesize,0.1,2)

    if(p==3){
      x3 = runif(samplesize,0.1,2)
      x4 = runif(samplesize,0.1,2)
    }

    if(p==5){
      x3 = runif(samplesize,0.1,2)
      x4 = runif(samplesize,0.1,2)
      x5 = runif(samplesize,0.1,2)
      x6 = runif(samplesize,0.1,2)
    }

    if(p==1){
      y = rexp(samplesize,exp(alpha1*x1+1*x2))
      dat = data.frame(cbind(y,x1,x2))
    }

    if(p==3){
      y = rexp(samplesize,exp(alpha1*x1+1*x2+1*x3+1*x4))
      dat = data.frame(cbind(y,x1,x2,x3,x4))
    }

    if(p==5){
      y = rexp(samplesize,exp(alpha1*x1+1*x2+1*x3+1*x4+1*x5+1*x6))
      dat = data.frame(cbind(y,x1,x2,x3,x4,x5,x6))
    }


    id.nonfac <- which(dat$x1 == 0)
    id.fac <- which(dat$x1 == 1)
    compare <- expand.grid(id.nonfac, id.fac)
    y1 <- dat[id.nonfac,"y"]
    y2 <- dat[id.fac,"y"]
    pseudo_y <- as.vector(t(outer(y1, y2, `<`))*1)
    pseudo_dat = data.frame("y"=pseudo_y)

    for(i in 1:p){
      x_g1 <- dat[id.nonfac,paste0("x",i+1)]
      x_g2 <- dat[id.fac,paste0("x",i+1)]
      pseudo_x <- as.vector(t(outer(x_g1, x_g2, `-`)))
      pseudo_dat[,paste0("x",i+1)] = -pseudo_x
    }

    pseudo_dat$C1 = rep((1:(samplesize/2)),each=samplesize/2)
    pseudo_dat$C2 = rep((1:(samplesize/2)),samplesize/2)
    pseudo_dat$C1C2 = 1:nrow(pseudo_dat)

    pseudo_form <- as.formula(paste("y~", paste(rep("x",p), (1:p)+1, sep="", collapse="+")))

    return(list(data=dat,form=as.formula(paste("y~", paste(rep("x",p), (1:(p+1)), sep="", collapse="+"))),pseudo_dat=pseudo_dat,pseudo_form=pseudo_form,comp=compare))
  }

}

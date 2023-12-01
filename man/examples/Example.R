require(devtools)
install_github("JaspersStijn/SmallSamplePIMFinal",force=TRUE);

library("SmallSamplePIM")

# Data example from Jaspers et al. (2023)

# Load packages
library(tidyverse)
library(MASS)
library(ggplot2)

## Plot data

ggplot(data=dat,aes(x=Dose,y=weight_use,fill=Sex))+geom_boxplot()+
  scale_fill_grey(start=0.4)+ylab("Weight difference after 90 days")+
  xlab("Pyridine Concentration (in ppm)")+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))


## Compare weight difference in specific subgroup, controlling for extra covariates
control = c("CreatT5","TPT5","ALBT5")
correc = "Firth"
response = "weight_use"
dat= dat %>% mutate_at(control, as.numeric)

subset_to_be_fitted = as.data.frame(subset(dat,Dose%in%c(0,50)&(Sex=="M")))

fit = GEE_MH_fit(data=subset_to_be_fitted,
                 response=response,
                 treatment = "Dose",
                 control = control,
                 correction = correc)


## get all estimates and confidence intervals for the data example
expit = function(x){return(exp(x)/(1+exp(x)))}

doses = unique(dat$Dose)
PI=c()
k=0
for(m in unique(dat$Sex)){
  for(i in 1:(length(doses)-1)){
    for(j in (i+1):length(doses)){
      k=k+1
      fit = GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(doses[i],doses[j])&(Sex==m))),response=response,treatment = "Dose", control = control,correction = correc)
      cat(round(expit(fit$beta_est)*100,2),
          paste(paste0(" (",round(expit(fit$lower)*100,2)),paste0(round(expit(fit$upper)*100,2),")"),sep=";"),"\n")
      PI[k] = paste0(round(expit(fit$beta_est)*100,2),
          paste(paste0(" (",round(expit(fit$lower)*100,2)),paste0(round(expit(fit$upper)*100,2),")"),sep=";"))

      }
  }
}


## get all adjusted p-values for the data example

doses = unique(dat$Dose)
pvals=c()
for(m in unique(dat$Sex)){
  for(i in 1:(length(doses)-1)){
    for(j in (i+1):length(doses)){
    pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(doses[i],doses[j])&(Sex==m))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
    }
  }
}
adjusted_p = p.adjust(pvals,method="BH")
round(adjusted_p,4)



for(i in 1:length(pvals)){
  cat(PI[i]," \\newline ",round(pvals[i],4)," \\newline ",round(adjusted_p[i],4),"\n")
}



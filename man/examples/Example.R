# Data example from Jaspers et al. (2023)

# Load packages
library(readxl)
library(tidyverse)
library(MASS)

source("./R/Main_Function.R")

## Read in data
dat = read_xlsx("./Data/Pyridine.xlsx")

## Data formatting

dat = dat %>% mutate(Dose = as.factor(Dose),
                     Sex = factor(Sex,labels=c("M","F")),
                     weight_use = as.numeric(weight_use))

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


## get all adjusted p-values for the data example

pvals=c()
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(0,50)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(0,100)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(0,250)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(0,500)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(0,1000)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(50,100)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(50,250)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(50,500)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(50,1000)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(100,250)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(100,500)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(100,1000)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(250,500)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(250,1000)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(500,1000)&(Sex=="M"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])

pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(0,50)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(0,100)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(0,250)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(0,500)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(0,1000)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(50,100)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(50,250)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(50,500)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(50,1000)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(100,250)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(100,500)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(100,1000)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(250,500)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(250,1000)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])
pvals=c(pvals,GEE_MH_fit(data=as.data.frame(subset(dat,Dose%in%c(500,1000)&(Sex=="F"))),response=response,treatment = "Dose", control = control,correction = correc)$pval[1])

adjusted_p = p.adjust(pvals,method="BH")
round(adjusted_p,4)



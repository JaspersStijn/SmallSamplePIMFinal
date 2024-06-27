# Introduction

R package related to 'Covariate-adjusted generalized pairwise comparisons in small samples', by Jaspers, Verbeeck and Thas (2024). Semiparametric Probabilistic Index Models allow for the comparison of two groups of observations, whilst adjusting for covariates, thereby fitting nicely within the framework of Generalized Pairwise Comparisons. In this paper, we
show that the parameters of the Probabilistic Index Model can be estimated using Generalized Estimating Equations, for which adjustments exist that lead to estimators of the
sandwich variance-covariance matrix with improved finite sample properties and that can deal with bias due to separa-
tion. In this way, appropriate inference.

# Usage of the package

Package can be installed by running

    require(devtools)
    install_github("JaspersStijn/SmallSamplePIMFinal",force=TRUE);
    library("SmallSamplePIM")

A self-explaining data example is provided in the R file 'Example.R' contained in the man->examples folder. The main function to be used is 

        GEE_MH_fit(data,response,treatment,control,correction,link)

The dataset contains a response variable, a treatment indicator (specifying which groups to be compared) and possible control covariates. The default correction used is MBN, corresponding to the Morel, Bokossa and Neerchal (2003) small sample correction. Other corrections are mentioned in the accompanying paper. The logit link is used as a deafult.

        

###########################################################################################
###########################################################################################
###########################################################################################
# NOTES:
#
# 1. This is a tutorial for fitting the lassoPL model described in Arbet et al (2017)
#
# 2. glmnet_lassoPL_functions.R contains the necessary functions (with helpful annotations)
#    in order to run the glmnet_lassoPL_tutorial code.
#
# 3. simdata.RData contains the simulated dataset used in this tutorial.
#    "Y" is the continuous response vector
#    "Xmat" is the 600 x 905 standardized SNP matrix.  Columns 1-5 represent the 5 causal SNPs (each with true regression coeffs equal to 0.15) 
#           columns 6-905 represent the 900 null SNPs (true reg. coeffs equal to 0)
#
# 4. glmnet_lassoPL_functions.R and glmnet_lassoPL_tutorial.R should only be used with 
#    small to moderate sized datasets (e.g. less than 100,000 predictors).  For larger 
#    datasets, use biglassoPL_functions.R and biglassoPL_tutorial.R instead
#
# REFERENCES:
#
# Arbet et al. "Resampling-based tests for Lasso in genome-wide association studies." BMC Genetics (Accepted June 2017, soon to be published)
#
# Wu, Tong Tong, et al. "Genome-wide association analysis by lasso penalized logistic regression." Bioinformatics 25.6 (2009): 714-721.
#
#
# QUESTIONS? email Jaron Arbet at: arbet003@umn.edu
###########################################################################################
###########################################################################################
###########################################################################################


source("glmnet_lassoPL_functions.R")

set.seed(1)
load("simdata.RData")

##################
#fit1: fits a lassoPL model using significant level 0.01
##################

fit1=lassoPL(Xmat=Xmat,Y=Y,repeats=50,alpha=0.01)
fit1$permLams #shows the 50 estimated values for the target value of lambda that controls the type-1-error rate at level alpha
fit1$lambda #final estimate used for lambda that controls the type-1-error rate at level alpha
fit1$DF #total number of significant predictors at level alpha (doesnt include intercept)
fit1$significant_predictors #all predictors whose regression coeffients are significantly different than 0, using significance testing level alpha

##################
# fit2: fits a lassoPL model using bonferroni correction: alpha= 0.05/ncol(Xmat)
##################
fit2=lassoPL(Xmat=Xmat,Y=Y,repeats=50,alpha=0.05/ncol(Xmat))
fit2$permLams #shows the 50 estimated values for the target value of lambda that controls the type-1-error rate at level alpha
fit2$lambda #final estimate used for lambda that controls the type-1-error rate at level alpha
fit2$DF #total number of significant predictors at level alpha (doesnt include intercept)
fit2$significant_predictors #all predictors whose regression coeffients are significantly different than 0, using significance testing level alpha

##################
#Diagnostics: 
#
# in order to ensure an accurate estimate for the target value of lambda,
# the left plot should converge to some constant, and the right plot should converge to 0
# as the number of 'repeats' used in the lassoPL model increases.
# see pg 12 and Figure 4 of Arbet et al (2017) for more info.  If these 2 conditions are not
# satisfied, then try increasing the number of "repeats" used in the lassoPL model
##################
fit1diag=PLdiag(fit1,Xmat,Y,plot=T)
fit2diag=PLdiag(fit2,Xmat,Y,plot=T)



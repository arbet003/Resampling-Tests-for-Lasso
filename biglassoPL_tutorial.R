###########################################################################################
###########################################################################################
###########################################################################################
# NOTES:
#
# 1. This is a tutorial for fitting the lassoPL model described in Arbet et al (2017)
#
# 2. biglassoPL_functions.R contains the necessary functions (with helpful annotations)
#    in order to run the biglassoPL_tutorial code.
#
# 3. simdata.RData contains the simulated dataset used in this tutorial.
#    "Y" is the continuous response vector
#    "Xmat" is the 600 x 905 standardized SNP matrix.  Columns 1-5 represent the 5 causal SNPs (each with true regression coeffs equal to 0.15) 
#           columns 6-905 represent the 900 null SNPs (true reg. coeffs equal to 0)
#
# 4. In real data applications, biglassoPL_tutorial.R and biglassoPL_functions.R should 
#    only be used for large datasets (e.g. more than 100,000 predictors).  For smaller
#    to moderate sized datasets (e.g. less than 100,000 predictors), the 
#    glmnet_lassoPL_functions.R and glmnet_lassoPL_tutorial.R should be used instead.
#    However, since I do not want to store a massive dataset on my github account,
#    for the purposes of this tutorial, a smaller dataset will be used to demonstrate
#    the biglassoPL() function.
#
# REFERENCES:
#
# Arbet et al. "Resampling-based tests for Lasso in genome-wide association studies." BMC Genetics (Accepted June 2017, soon to be published)
# Wu, Tong Tong, et al. "Genome-wide association analysis by lasso penalized logistic regression." Bioinformatics 25.6 (2009): 714-721.
# Yi, Hui, et al. "Penalized multimarker vs. single-marker regression methods for genome-wide association studies of quantitative traits." Genetics 199.1 (2015): 205-222.
#
#
# QUESTIONS? email Jaron Arbet at: arbet003@umn.edu
###########################################################################################
###########################################################################################
###########################################################################################


######################################
#Step 1: Load functions and Data
######################################
source("biglassoPL_functions.R") #requires packages: glmnet, bigmemory, biglasso, doSNOW
load("simdata.RData")
alpha=0.05/ncol(Xmat) #bonferroni correction
set.seed(1)

######################################
#Step 2: create big.matrix object
#
#this is a ONE TIME COMPUTATIONAL COST, after creating the big.matrix object,
#you can instantaneously load the big.matrix object in future R sessions
#at no computational cost
#
#if you get error "problem creating filebacked matrix" it is likely because you have
#already run this code and created the file-backed dmatrix.
#If this is the case, then ignore Step 2 and proceed to Steps 3-4
#
#After creating a file-backed big.matrix object with a particular backingfile 
#and descriptorfile name (e.g. bigXmat.bin and bigXmat.desc),
#you cannot create a NEW big.matrix object using these same file names.  If you wish to
#create a new big.matrix object, you must use new file names, or you could delete the old
#.bin and .desc files if you wish to use the same file names.
######################################

bigXmat=as.big.matrix(Xmat,backingfile="bigXmat.bin",descriptorfile="bigXmat.desc",shared=T)



######################################
#Step 3: find reasonable window for target value of lambda using glmnet
#
#This is a ONE TIME COMPUTATIONAL COST and requires loading the entire SNP matrix into RAM.
#However, if one does not want to load the entire SNP matrix into RAM, then:
#
#
#*******alternatively, one can use the following analytic formula (modified from Yi et al 2015
#       to find a reasonable window for the target value of lambda:
#
#       lambda.min= qnorm(high_alpha)*sd(Y)/(-sqrt(length(Y)))
#       lambda.max= qnorm(low_alpha)*sd(Y)/(-sqrt(length(Y)))
#
#       When alpha>=1/p, then set low_alpha and high_alpha to 1 or 2 order of
#       magnitudes below and above the target alpha. For example, suppose you want to
#       control the type-1-error rate at level alpha=0.01, then you could use 
#       low_alpha=0.001 and high_alpha=0.1, which for our simulated
#       dataset would result in: lambda.min= 0.05804 and lambda.max=0.13996
#      
#       However, when alpha<1/p, then set high_alpha = 10*(1/p) and low_alpha= (1/p)/10
#       which for our simulated dataset would result in:
#       lambda.min=0.10365 and lambda.max=0.16729
######################################
LW=find_lambda_window(Xmat,as.vector(Y),alpha=alpha)
LW$lambda.min
LW$lambda.max

#alternatively, if one does NOT want to load the entire SNP matrix into RAM, then 
#the following analytic formula can be used to get a reasonable window for lambda
#(see above notes for more info)

lambda.min.alt=qnorm(10/ncol(Xmat))*sd(Y)/(-sqrt(length(Y)))
lambda.max.alt=qnorm(1/(ncol(Xmat)*10))*sd(Y)/(-sqrt(length(Y)))

#######################################
#Step 4: Fit biglassoPL model that controls type-1-error rate at level 
#        alpha=0.05/ncol(Xmat), i.e. bonferroni correction
#
#Uses 4 cores to estimate the target value of lambda 20 times in parallel
#
#Note given that our sample simulated dataset is small (only 900 predictors),
#glmnet_lassoPL_functions.R will be significantly faster, however, for larger datasets
#(e.g. more than 100,000 predictors), biglassoPL_functions.R will likely be preferable
#because the entire SNP matrix does not need to be loaded into RAM
#
#######################################
bigXmat=attach.big.matrix("bigXmat.desc")
fit=biglassoPL(bigXmat,Y,cores=4,repeats=20,lambda.min=LW$lambda.min,lambda.max=LW$lambda.max,alpha=alpha)
fit$significant_predictors#all predictors whose regression coeffients are significantly different than 0, using significance testing level alpha
fit$permLams #shows the 20 estimated values for the target value of lambda that controls the type-1-error rate at level alpha.  The median of these 20 estimates is used as the final estimate for lambda.
fit$lambda.final #final value of lambda used that controls the type-1-error rate at level alpha


#######################################
#Step 5: Diagnostics: make sure you obtained an accurate estimate of the
#                     value of lambda that controls the type-1-error rate at level alpha
#
# in order to ensure an accurate estimate for the target value of lambda,
# the left plot should converge to some constant, and the right plot should converge to 0
# as the number of 'repeats' used in the biglassoPL model increases.
# see pg 12 and Figure 4 of Arbet et al (2017) for more info.  If these 2 conditions are not
# satisfied, then try increasing the number of "repeats" used in the lassoPL model
#######################################
fitdiag=PLdiag(fit,bigXmat,Y,cores=4)
plot_PLdiag(fitdiag)
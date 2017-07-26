###########################################################################################
###########################################################################################
###########################################################################################
# NOTES:
#
# 1. glmnet_lassoPL_functions.R contains the necessary functions in order to run the
#    glmnet_lassoPL_tutorial code.
#
# 2. glmnet_lassoPL_functions.R should only be used with small to moderate sized
#    datasets (e.g. less than 100,000 predictors).  For larger datasets,
#    use biglassoPL_functions.R instead
#
# 3. The user will primarily use only 2 functions within this file: lassoPL() to fit the lassoPL
# model described in Arbet et al (2017), and the PLdiag() diagnostic function to ensure
# that a good estimate for the value of lambda that controls the type-1-error rate at the
# desired level has been obtained.  
#
# 4. All functions within this file contain helpful annotations.
#
# 5. glmnet_lassoPL_functions.R does NOT use parallel computing, as it is unnecessary for
#    moderate sized datasets.  The biglassoPL files DO use parallel computing to estimate the target lambda
#
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




options(scipen=c(999))

##########################################################################
#load packages
##########################################################################
library("glmnet")

##########################################################################
# Misc. Functions
##########################################################################

is_TE = function(x){  #see if "object" has error or not
  inherits(x, "try-error")
}


##########################################################################################
# bisect(): Function to select lambda for Lasso that allows for an exact number of nonzero 
#           coeffs. Uses bisection algorithm of Wu 2009.
#
# Xmat: matrix of standardized SNPs
# Y: continuous response
# nonzero: the exact number of nonzero coeffs you want
# 
# returns the value of lambda that gives the desired number of nonzero coefficients for Lasso
##########################################################################################

bisect=function(Xmat,Y,nonzero=1){

    fit=glmnet(Xmat,Y,alpha=1,dfmax=nonzero+100)
    LDF=cbind(fit$df,fit$lambda)
    tmax=LDF[,2][max(which(LDF[,1]>=0& LDF[,1]<nonzero))]
    tmin= try(LDF[,2][min(which(LDF[,1]>=nonzero))])
    #if(is_TE(tmin)){
      #print("error: unable to calculate tmin, try increasing dfmax or nlam?")
    #}
    if(tmin<=0 |is.na(tmin)){
      tmin=0.00000000000001
    }
    if(tmax>100){
      count3=0
      dfmax2=nonzero+100
      while(count3<1){
        #print(paste("warning: bad starting lambda values (too high), will try again"))
        dfmax2=dfmax2+5
        fit=glmnet(Xmat,Y,alpha=1,dfmax=dfmax2)
        LDF=cbind(fit$df,fit$lambda)
        #print(LDF)
        tmax=LDF[,2][max(which(LDF[,1]>=0& LDF[,1]<nonzero))]
        tmin= try(LDF[,2][min(which(LDF[,1]>=nonzero))])
        #if(is_TE(tmin)){
          #print("error: unable to calculate tmin, try increasing dfmax")
        #}
        if(tmin<=0 |is.na(tmin)){
          tmin=0.00000000000001
        }
        if(tmax<10){
          count3=1
          #print("found good starting lambdas, will resume")
        }
      }
    }
    
    
    tL=tmin
    tU=tmax
    #print(c(tL,tU))
    tm=NA
    tm2=NA
    finalDF=NA
    lam=NA
    mdf=NA
    
    count=0
    count2=0
    while(count<1){
      if(!is.na(tm)){tm2=tm}
      fitL=glmnet(Xmat,Y,lambda=tL)
      fitU=glmnet(Xmat,Y,lambda=tU)
      if(nonzero>fitL$df){#print("problem, all possible lam values result in < desired nonzero coeffs")
        tL=0.00000000000001
      }
      
      if(nonzero<fitU$df){#print("problem, all possible lam values result in > desired nonzero coeffs")
        tU= tU+quantile(abs(Y),0.1)
      }
      if(fitL$df==nonzero){lambda=fitL$lambda
      count=1
      finalDF=nonzero
      }
      if(fitU$df==nonzero){lambda=fitU$lambda
      count=1
      finalDF=nonzero
      }
      if(fitU$df<nonzero&nonzero<fitL$df){tm= 0.5*(tL+tU)
      fitM=glmnet(Xmat,Y,lambda=tm)
      lam=fitM$lambda
      mdf=fitM$df
      if(fitM$df==nonzero){lambda=fitM$lambda
      count=count+1
      finalDF=fitM$df
      } else if(fitM$df<nonzero){
        tU=tm
      } else if(fitM$df>nonzero){
        tL=tm
      }
      }
      
      if(!is.na(tm)& !is.na(tm2) & tm==tm2){
        #count2= count2+1
        count2=5
      }
      else if(!is.na(tm)& !is.na(tm2) & tm!=tm2){
        count2=0
      }
      if(count2==5){
        count=1
        lambda=fitM$lambda
        finalDF=fitM$df
        print("tm has not changed for 2 iterations, will try again")
      }
      #print("current lambda and df:")
      #print(c(lam,mdf))
    }
    return(c(lambda,finalDF))
}

##################################################
# permlam(): Function that uses permutations to select lambda that controls overall 
#            type-1-error rate at level alpha
#
# Xmat: matrix of standardized SNPs
# Y: continuous response
# alpha: desired type-1-error rate. If alpha < 1/p, then this program can only 
#        handle alpha values of the following form: alpha= 1/(p*k) for some integer k,
#        and p=ncol(X)
#
# repeats: total number of times you want to estimate the target lambda
#
# returns the estimates of lambda that control the type-1-error at level alpha. For example,
#         if "repeats"= 3 then permlam() will estimate the target value of lambda 3 times,
#         then the lassoPL() function will fit the actual lassoPL model by using
#         the median of these 3 lambda estimates for the final estimate of lambda
##################################################

permlam=function(Xmat,Y,alpha=1/ncol(Xmat),repeats=3){
  ##### alpha >= 1/p
  if(alpha>=1/ncol(Xmat)){
    nonzero= round(alpha*ncol(Xmat))
    lams=vector()
    count=0
    start=Sys.time()
    while(count<repeats){
      ind=sample(1:length(Y),length(Y))
      y.perm=Y[ind]
      a=bisect(Xmat,y.perm,nonzero=nonzero)
      if(a[2]==nonzero){
        count=count+1
        end=Sys.time()
        print(paste("Done with repeat",count,", time=",round(end-start,3)))
        lams=rbind(lams,c(a,end-start))
        start=Sys.time()
      }
    }
    colnames(lams)=c("lambda","perm DF","time")
    return(lams)
  }
  
  #### alpha < 1/p
  if(alpha< 1/ncol(Xmat)){
    if(as.integer(round(1/alpha,0))%%ncol(Xmat)!=0) {
      stop("If alpha < 1/p, then this program can only handle alpha values of the following form: alpha= 1/(p*k) for some integer k, and p=ncol(X)")
    }
    nonzero=1
    reps2= 1/(alpha*ncol(Xmat))
    maxlams=vector()
    for(w in 1:repeats){
      lams=vector()
      count=0
      start=Sys.time()
      while(count<reps2){
        ind=sample(1:length(Y),length(Y))
        y.perm=Y[ind]
        a=bisect(Xmat,y.perm,nonzero=nonzero)
        if(a[2]==nonzero){
          count=count+1
          #print(paste("Done with repeat",count,", time=",round(end2-start2,3)))
          lams=rbind(lams,a)
        }
      }
      colnames(lams)=c("lambda","perm DF")
      end=Sys.time()
      maxlams=rbind(maxlams,c(max(lams[,1]),round(end-start,3),all(lams[,2]==1)))
      #print(maxlams[w,])
      print(paste("repeat",w,"finished, time=",round(end-start,3)))
    }
    colnames(maxlams)=c("lambda","time","DF Check")
    return(maxlams)
  }
}


##########################################################################
#  LassoPL(): fit a lassoPL model from Arbet et al (2017), i.e. use permutations to select
#             the value of lambda that controls the type-1-error rate at level alpha
#
# Xmat: matrix of standardized SNPs
# Y: continuous respons
# lambda: if one already knows the value of lambda they wish to use (e.g. from a previous lassoPL model)
#         then you can simply specify lambda here
# repeats: total number of times you want to estimate the target value of lambda, then
#          will use the median as the final estimate for the target lambda
#
# alpha: desired type-1-error rate. If alpha < 1/p, then this program can only 
#        handle alpha values of the following form: alpha= 1/(p*k) for some integer k,
#        and p=ncol(X)
#
#returns: betahat=estimated coefficients
#         significant_predictors: all predictors with regression coefficients that are significantly different than 0 using significance level alpha
#         DF=estimated number of degrees of freedom (nonzero coefficients), not counting intercept
#         lambda.final= final estimate of value of lambda that controls type-1-error rate at level alpha
#         permLams= output from permlam() function, contains all estimates for the target value of lambda,
#                   For example, if "repeats"= 3 then permlam() will estimate the target value of lambda 3 times,
#                   then the lassoPL() function will fit the actual lassoPL model by using
#                   the median of these 3 lambda estimates for the final estimate of lambda
#
#         time= total time taken to fit the lassoPL model
#########################################################################

lassoPL = function(Xmat, Y,lambda=NULL,repeats=10,alpha=0.05){
  start=Sys.time()
  Y=as.vector(Y)
  permLams=NULL
  
  if(!is.null(lambda)){
    model = glmnet(Xmat, Y,  alpha=1,lambda = lambda)
    lambda.final=lambda
  }
  
  if(is.null(lambda)){
      temp=permlam(Xmat,Y,repeats=repeats,alpha=alpha)
      lambda.final=median(temp[,1],na.rm=T)
      permLams=temp
      model = glmnet(Xmat, Y,alpha=1,lambda = lambda.final)
  }
  
  beta = coef(model)
  beta=as.vector(beta)
  names(beta)=c("Intercept",colnames(Xmat))
  end=Sys.time()
  time= end-start
  
  significant_predictors= beta[-1][which(beta[-1]!=0)]
  
  stats = list(betahat = beta,significant_predictors=significant_predictors ,DF=length(which(model$beta!=0)),
               lambda.final=lambda.final,permLams=permLams,time=time)
  
  return(stats)
}
#########################################################################################
# PLdiag(): diagnostic function for lassoPL, see pg 12 and Figure 4 of Arbet et al (2017).
#          In order to obtain a good estimate of the target value of lambda, the number of 
#          nonzero coefficients should converge to some constant (left plot),
#          and the number of discrepant SNPs should converge to 0 (right plot.  
#          If these 2 conditions are not met, try increasing the number of "repeats" used
#          in the lasso-PL model.
#
# model: a lassoPL() object
# X: matrix of standardized SNPs
# Y: continuous response
# plot: when = TRUE, will create the 2 plots described above
#
# returns object "stats":
#         the kth row of stats[,1] shows the total number of significant SNPs using k number of 'repeats' to estimate lambda in the lassoPL model...ideally should converge to some constant as the number of repeats increases
#         the kth row of stats[,2] shows the total number of discrepancies between models using k and (k-1) number of repeats to estimate lambda....ideally should converge to 0 as the number of repeats increases
#
#
#########################################################################################
PLdiag=function(model,X,Y,plot=T){
  stats=vector()
  lam=model$permLams[1,1]
  sig=list()
    lass=lassoPL(X,Y,lambda=lam)
    sig[[1]]=which(lass$beta[-1]!=0)
    miss=0

  stats=rbind(stats,c(lam,length(sig[[1]]),miss))
  
  for(i in 2:nrow(model$permLams)){
    lam=mean(model$permLams[1:i,1])
      lass=lassoPL(X,Y,lambda=lam)
      sig[[i]]=which(lass$beta[-1]!=0)
      miss= length(which(!sig[[i-1]]%in%sig[[i]]))+length(which(!sig[[i]]%in%sig[[i-1]]))
    stats=rbind(stats,c(lam,length(sig[[i]]),miss))
  }
  
  colnames(stats)=c("Lambda(B)","Num. Nonzero Coeffs","Num. Discrepencies")
  #stats2=list(stats=stats,sig=sig)
  stats2=list(stats=stats)
  
  ####plots
  par(mfrow=c(2,2))
  plot(1:nrow(stats),stats[,2],main="Model Size",xlab="Num. of 'repeats' used in lassoPL",ylab="Num. of Significant SNPs")
  plot(1:nrow(stats),stats[,3],main="Model Discrepancies",xlab="Num. of 'repeats' used in lassoPL",ylab="Num. of Discrepancies")
  return(stats2)
}

###########################################################################################
###########################################################################################
###########################################################################################
# NOTES:
#
# 1. biglassoPL_functions.R contains the necessary functions in order to run the
#    biglassoPL_tutorial code.
#
# 2. In real data applications, biglassoPL_tutorial.R and biglassoPL_functions.R should 
#    only be used for large datasets (e.g. more than 100,000 predictors).  For smaller
#    to moderate sized datasets (e.g. less than 100,000 predictors), the 
#    glmnet_lassoPL_functions.R and glmnet_lassoPL_tutorial.R should be used instead.
#    However, since I do not want to store a massive dataset on my github account,
#    for the purposes of this tutorial, a smaller dataset will be used to demonstrate
#    the biglassoPL() function.
#
# 3. The user will primarily use only 4 functions within this file: the find_lambda_window()
#    function to find a reasonable window for the target value of lambda the gives the desired
#    type-1-error rate; biglassoPL() to fit the lassoPL model described in Arbet et al (2017);
#    and the PLdiag() and plot_PLdiag() diagnostic functions to ensure that a good estimate
#    for the value of lambda that controls the type-1-error rate at the desired level has
#    been obtained.  
#
# 4. All functions within this file contain helpful annotations.
#
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



library("glmnet")
library("biglasso")
library("doSNOW")


##############################################
#Misc. functions
##############################################
DF=function(i,mat){
  return(sum(mat[,i]!=0))
}

is_TE = function(x){  #see if "object" has error or not
  inherits(x, "try-error")
}
##############################################

##########################################################################################
# find_lambda_window(): Function to find a reasonable window for the target value of lambda
#                       that controls the type-1-error rate at level alpha.  This window
#                       is supplied to the biglassoPL function via the lambda.min and
#                       lambda.max arguments.
#
#This function uses glmnet, therefore the entire SNP matrix must be loaded in RAM.
#This is a one-time computational cost.  If it is not possible to load the entire 
#matrix in RAM, then one can ignore the find_lambda_window() function and simply
#supply your own guess at a reasonable window for the target lambda via the 
#lambda.min and lambda.max arguments of the biglassoPL() function.
#
# Xmat: matrix of standardized SNPs
# Y: continuous response
# alpha: the desired type-1-error rate
# nlambda: the number of candidate lambda values used
# 
# returns a reasonable window for the target value of lambda (lambda.min and lambda.max),
#     which is supplied to the biglassoPL function via the lambda.min and lambda.max arguments
##########################################################################################
find_lambda_window=function(Xmat,Y,dfmax=1000,alpha=1/ncol(Xmat),nlambda=100){
  if(alpha< 1/ncol(Xmat)){
    if(as.integer(round(1/alpha,0))%%ncol(Xmat)!=0) {
      stop("If alpha < 1/p, then this program can only handle alpha values of the following form: alpha= 1/(p*k) for some integer k, and p=ncol(X)")
    }
  nonzero=1
  } else {
    nonzero=round(alpha*ncol(Xmat))
  }

    ind=sample(1:length(Y),length(Y))
    Y2=Y[ind]
    Y2=Y
  start2=Sys.time()
  start=Sys.time()
  print("starting initial fit:")
  fit=glmnet(Xmat,Y2,nlambda=nlambda,dfmax=dfmax)
  end=Sys.time()
  print(paste("finished initial fit, time=",end-start))
  df=unlist(lapply(1:ncol(fit$beta),DF,mat=fit$beta))
  LDF=cbind(df,fit$lambda)
  #print(LDF)
  tmin= try(LDF[,2][min(which(LDF[,1]>= nonzero+30))])
  print(tmin)
  if(is_TE(tmin)|is.na(tmin)){
    stop("all lambdas fit to model result in  < desired number of nonzero coeffs.  Try again but increase nlambda")
  }
  else {
    print(paste("found lambda.min=",tmin))
  }
  tmax=LDF[,2][max(which(LDF[,1]>=0& LDF[,1]<nonzero))]
  if(tmax>1000){
    count=0
    while(count<1){
      start=Sys.time()
      print("bad initial lambda values, try again:")
      fit=glmnet(Xmat,Y2,lambda= tmin*2)
      end=Sys.time()
      print(paste("time=",end-start))
      df=unlist(lapply(1:ncol(fit$beta),DF,mat=fit$beta))
      #LDF=cbind(df,fit$lambda)
      print(LDF)
      tmax=LDF[,2][max(which(LDF[,1]==0))]
      print(tmax)
      if(!is.na(tmax)){
        count=1
      }
    }
  }
  
  print(paste("found lambda.max=",tmax))
  lambda.max=tmax
  end2=Sys.time()
  return(list(lambda.min=tmin,lambda.max=tmax,lambda_table=LDF,time= end2-start2))
}


##########################################################################################
# bisect3(): Function to select lambda for biglasso that allows for an exact number of nonzero 
#           coeffs. Uses bisection algorithm of Wu 2009.  This function is called automatically
#           within biglassoPL() and thus one should not need to directly use this function.
#
# Xmat: big.matrix object of standardized SNPs
# Y: continuous response
# nonzero: the exact number of nonzero coeffs you want
# lambda.min and lambda.max: use find_lambda_window() function then input its results here
# 
# returns the value of lambda that gives the desired number of nonzero coefficients for biglasso
##########################################################################################
bisect3=function(Xmat,Y,nonzero=1,dfmax= nonzero+20,lambda.min=NULL,lambda.max=NULL){
  start=Sys.time()
  descriptorfile=describe(Xmat)
  lambda.final=NA
  tmin=lambda.min
  tmax=lambda.max
  tL=tmin
  tU=tmax

  tm=NA
  tm2=NA
  finalDF=NA
  lam=NA
  mdf=NA
  
  count=0
  count2=0
  while(count<1){
    if(!is.na(tm)){tm2=tm}
    fitL=biglasso(Xmat,Y,lambda=tL)      
    fitU=biglasso(Xmat,Y,lambda=tU)
    if(nonzero>sum(fitL$beta[-1]!=0)){print("problem, all possible lam values result in < desired nonzero coeffs")
                                      tL=0.00000000000001
    }
    
    if(nonzero<sum(fitU$beta[-1]!=0)){print("problem, all possible lam values result in > desired nonzero coeffs")
                                      tU= tU+quantile(abs(Y),0.1)
    }
    
    if(sum(fitL$beta[-1]!=0)==nonzero){lambda=fitL$lambda
                                       count=1
                                       finalDF=nonzero
    }
    if(sum(fitU$beta[-1]!=0)==nonzero){lambda=fitU$lambda
                                       count=1
                                       finalDF=nonzero
    }
    
    if(!is.null(lambda.min)& !is.null(lambda.max) & sum(fitL$beta[-1]!=0)<nonzero){
      stop("initial window created by user specified lambda.min and lambda.max does not contain the target lambda value:lambda.min is too large... use a smaller lambda.min")
    }
    
    if(!is.null(lambda.min)& !is.null(lambda.max) & sum(fitU$beta[-1]!=0)>nonzero){
      print("lambda.max too small, will increase")
      tU= lambda.max+quantile(abs(Y),0.1)
      
    }
    
    
    if(sum(fitU$beta[-1]!=0)<nonzero&nonzero<sum(fitL$beta[-1]!=0)){tm= 0.5*(tL+tU)
                                                                    fitM=biglasso(Xmat,Y,lambda=tm)
                                                                    lam=fitM$lambda
                                                                    mdf=sum(fitM$beta[-1]!=0)
                                                                    if(mdf==nonzero){lambda=fitM$lambda
                                                                                     count=count+1
                                                                                     finalDF=sum(fitM$beta[-1]!=0)
                                                                    } else if(mdf<nonzero){
                                                                      tU=tm
                                                                    } else if(mdf>nonzero){
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
      finalDF=mdf
      print("tm has not changed for 2 iterations, will stop")
    }
    #print("current lambda and df:")
    #print(c(lam,mdf))
  }
  end=Sys.time()
  return(c(lambda,finalDF,end-start))
  
}


##################################################
# permlam3(): Function that uses permutations to select lambda that controls overall 
#            type-1-error rate at level alpha.  This function is called automatically
#            within biglassoPL() and thus the user should not need to directly use this function.
#
# Xmat: big.matrix object of standardized SNPs
# Y: continuous response
# alpha: desired type-1-error rate. If alpha < 1/p, then this program can only 
#        handle alpha values of the following form: alpha= 1/(p*k) for some integer k,
#        and p=ncol(X)
# dfmax = maximum number of nonzero coefficients allowed (can be used to simplify computation)
# lambda.min and lambda.max: obtained from the find_lambda_window() function
#
# repeats: total number of times you want to estimate the target lambda
# cores: number of cores used in parallel
#
# returns the estimates of lambda that control the type-1-error at level alpha. For example,
#         if "repeats"= 3 then permlam3() will estimate the target value of lambda 3 times,
#         then the biglassoPL() function will fit the actual lassoPL model by using
#         the median of these 3 lambda estimates for the final estimate of lambda
##################################################
permlam3=function(Xmat,Y,alpha=1/ncol(Xmat),dfmax=100,lambda.min=NULL,lambda.max=NULL,repeats=4,cores=4){
  progress <- function(iter) {print(paste("Finished Repeat Num. ",iter,sep=""))}
  opts <- list(progress = progress)
  
  #uses bisect3
  descriptorfile=describe(Xmat)
  if(is.null(lambda.min)|is.null(lambda.max)){
    stop("must specify reasonable starting window (lambda.min and lambda.max) for target lambda.  Use find_lambda_window() function to find starting window")
  }
  if(alpha>=1/ncol(Xmat)){
    nonzero= round(alpha*ncol(Xmat))
    lams=vector()
    lams=foreach(i=1:repeats,.combine=rbind,.options.snow = opts)%dopar%{
      Xmat <- attach.big.matrix(descriptorfile)
      count=0
      start=Sys.time()
      while(count<1){
        ind=sample(1:length(Y),length(Y))
        y.perm=Y[ind]
        a=bisect3(Xmat,y.perm,nonzero=nonzero,dfmax=dfmax,lambda.min=lambda.min,lambda.max=lambda.max)
        if(a[2]==nonzero){
          count=count+1
          #print(paste("Done with repeat",count,", time=",round(end-start,3)))
          #lams=rbind(lams,c(a,end-start))
        }
      }
      end=Sys.time()
      print(paste("done with repeat",i,", time=",end-start))
      return(c(a[1:2],end-start))
    }
    #colnames(lams)=c("lambda","perm DF","time")
    return(lams)
  }
  
  #### alpha < 1/p
  if(alpha< 1/ncol(Xmat)){
    if((1/alpha)%%ncol(Xmat)!=0) {
      stop("If alpha < 1/p, then this program can only handle alpha values of the following form: alpha= 1/(p*k) for some integer k, and p=ncol(X)")
    }    
    nonzero=1
    reps2= 1/(alpha*ncol(Xmat)) 
    maxlams=foreach(w=1:repeats,.combine=rbind,.options.snow=opts)%dopar%{
      Xmat <- attach.big.matrix(descriptorfile)
      lams=vector()
      count=0
      start=Sys.time()
      while(count<reps2){
        ind=sample(1:length(Y),length(Y))
        y.perm=Y[ind]
        a=bisect3(Xmat,y.perm,nonzero=nonzero,dfmax=dfmax,lambda.min=lambda.min,lambda.max=lambda.max)
        if(a[2]==nonzero){
          count=count+1
          lams=rbind(lams,a)
        }
      }
      #colnames(lams)=c("lambda","perm DF")
      end=Sys.time()
      print(paste("done with repeat",w,"time=",end-start))
      return(c(max(lams[,1]),round(end-start,3),all(lams[,2]==1)))
    }
    colnames(maxlams)=c("lambda","time","DF Check")
    return(maxlams)
  }
}

##########################################################################
#  LassoPL(): fit a lassoPL model from Arbet et al (2017), i.e. use permutations to select
#             the value of lambda that controls the type-1-error rate at level alpha
#
# Xmat: big.matrix of standardized SNPs
# Y: continuous respons
# lambda: if one already knows the value of lambda they wish to use (e.g. from a previous biglassoPL model)
#         then you can simply specify lambda here
# repeats: total number of times you want to estimate the target value of lambda, then
#          will use the median as the final estimate for the target lambda
#
# alpha: desired type-1-error rate. If alpha < 1/p, then this program can only 
#        handle alpha values of the following form: alpha= 1/(p*k) for some integer k,
#        and p=ncol(X)
#
#cores: number of cores to be used in parallel
#
#returns: betahat=estimated coefficients
#         #significant_predictors: all predictors with regression coefficients that are significantly different than 0 using significance level alpha
#         DF=estimated number of degrees of freedom (nonzero coefficients), not counting intercept
#         lambda.final= final estimate of value of lambda that controls type-1-error rate at level alpha
#         permLams= output from permlam3() function, contains all estimates for the target value of lambda,
#                   For example, if "repeats"= 3 then permlam() will estimate the target value of lambda 3 times,
#                   then the biglassoPL() function will fit the actual lassoPL model by using
#                   the median of these 3 lambda estimates for the final estimate of lambda
#
#         time= total time taken to fit the lassoPL model
#########################################################################
biglassoPL = function(Xmat, Y,lambda=NULL,lambda.min=NULL,lambda.max=NULL,dfmax=1000,repeats=4,alpha=1/ncol(Xmat),cores=4){
  #uses bisect3 - Wu's bisection Alg
  #can estimate 20 lambda_alphas in parallel
  start=Sys.time()
  permLams=NULL
  
  descriptorfile=describe(Xmat)
  
  if(!is.null(lambda)){
    model = biglasso(Xmat, Y,lambda=lambda)
    lambda.final=lambda
  }  
  
  if(is.null(lambda)){
    print(paste("setting up cluster with",cores,"cores"))
    progress <- function(iter) {print(paste("Finished Repeat Num. ",iter,sep=""))}
    opts <- list(progress = progress)
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    clusterExport(cl, c("%dopar%","foreach","bisect3","is_TE","biglasso","describe","DF","attach.big.matrix"))
    print("finished setting up cluster")

      #if(parallel==F){
      #  temp=permlam(Xmat,Y,alpha=alpha,repeats=repeats,intercept=intercept,fit_type="lasso",standardize=standardize)
      #  lambda.final=mean(temp[,1],na.rm=T)
      #  permLams=temp
      #  model = glmnet(Xmat, Y,penalty.factor=penalty.factor,alpha=1,lambda = lambda.final,intercept=intercept,standardize=standardize)
      #}
      temp=permlam3(Xmat,Y,lambda.min=lambda.min,dfmax=dfmax,lambda.max=lambda.max,alpha=alpha,repeats=repeats)
      if(repeats==1){
        lambda.final=temp[1]
      }
      else {
        lambda.final=median(temp[,1])
      }
      permLams=temp
      colnames(permLams)=c("lambda_alpha_hat","time","DF")
      model = biglasso(Xmat, Y,lambda = lambda.final)
      stopCluster(cl)
    }
  
  #start2=Sys.time()
  #fitted.vals= Xmat%*%model$beta[-1] + model$beta[1]
  #end2=Sys.time()
  #print(paste("calculated fitted values: time=",end2-start2))
  beta=as.vector(model$beta)
  names(beta)=c("Intercept",colnames(Xmat))
  
  significant_predictors= beta[-1][which(beta[-1]!=0)]
  
  #residuals=Y-fitted.vals
  end=Sys.time()
  time= end-start

  stats = list(betahat = beta,significant_predictors=significant_predictors, DF=length(which(model$beta[-1]!=0)),
               lambda.final=lambda.final,permLams=permLams,time=time,alpha=alpha
               
  )
  return(stats)
}

#########################################################################################
# PLdiag(): diagnostic function for biglassoPL, see pg 12 and Figure 4 of Arbet et al (2017).
#          In order to obtain a good estimate of the target value of lambda, the number of 
#          nonzero coefficients should converge to some constant (left plot),
#          and the number of discrepant SNPs should converge to 0 (right plot.  
#          If these 2 conditions are not met, try increasing the number of "repeats" used
#          in the lasso-PL model.
#
# model: a biglassoPL() object
# X: a big.matrix object of standardized SNPs
# Y: continuous response
#
# returns object "stats":
#         the kth row of stats[,1] shows the total number of significant SNPs using k number of 'repeats' to estimate lambda in the biglassoPL model...ideally should converge to some constant as the number of repeats increases
#         the kth row of stats[,2] shows the total number of discrepancies between models using k and (k-1) number of repeats to estimate lambda....ideally should converge to 0 as the number of repeats increases
#
#
#########################################################################################
PLdiag=function(model,X,Y,cores=4){
  progress <- function(iter) {print(paste("Finished Repeat Num. ",iter,sep=""))}
  opts <- list(progress = progress)
  cl2 <- makeCluster(cores)
  registerDoSNOW(cl2)
  clusterExport(cl2, c("%dopar%","foreach","bisect3","is_TE","biglasso","describe","biglassoPL","DF","attach.big.matrix"))
  descriptorfile=describe(X)  
  stats=foreach(i=1:nrow(model$permLams),.options.snow = opts)%dopar%{
    X <- attach.big.matrix(descriptorfile)
    lam=mean(model$permLams[1:i,1])
    fit=biglassoPL(X,Y,lambda=lam)
    sig=which(fit$beta[-1]!=0)
    #miss= length(which(!sig[[i-1]]%in%sig[[i]]))+length(which(!sig[[i]]%in%sig[[i-1]]))
    return(list(stat=c(lam,as.integer(length(sig))),sig=sig))
  }
  stat1=vector()
  sig=list()
  for(i in 1:length(stats)){
    stat1=rbind(stat1,unlist(stats[[i]][1]))
    sig[[i]]=as.vector(stats[[i]][2])
  }
  miss=rep(0,length(stats))
  for(i in 2:length(stats)){
    miss[i]= length(which(!sig[[i-1]]%in%sig[[i]]))+length(which(!sig[[i]]%in%sig[[i-1]]))
  }
  stats2=data.frame(stat1,miss)
  colnames(stats2)=c("Lambda(B)","Num. Nonzero Coeffs","Num. Discrepencies")
  stopCluster(cl2)
  
  return(list(stats=stats2))
}

#########################################################################################
# plot_PLdiag(): plots results from PLdiag function. In order to obtain a good estimate of 
#          the target value of lambda, the number of nonzero coefficients should converge 
#          to some constant (left plot), and the number of discrepant SNPs should converge
#          to 0 (right plot. If these 2 conditions are not met, try increasing the number 
#          of "repeats" used in the lasso-PL model.
#
# PLdiag_obj: an object from PLdiag()
#
#########################################################################################
plot_PLdiag=function(PLdiag_obj){
  par(mfrow=c(1,2))
  plot(1:nrow(PLdiag_obj[[1]]),PLdiag_obj[[1]][,2],main="Model Size",xlab="Num. of 'repeats' used in biglassoPL",ylab="Num. of Significant SNPs")
  plot(1:nrow(PLdiag_obj[[1]]),PLdiag_obj[[1]][,3],main="Model Discrepancies",xlab="Num. of 'repeats' used in biglassoPL",ylab="Num. of Discrepancies")
}


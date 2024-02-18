simdata<-function(n=n,m=m,a=a,b=bb,r=r,ppp=ppp,c=100000000,ll=10000000,seed=1,cov=NULL,g=NULL,lll=NULL){

  if(is.null(cov)==FALSE){
    return(simdata_cov(n=n,m=m,a=a,r=r,b=b,g=g,lll=lll,ppp=ppp,c=c,ll=ll,seed=seed))
  }else {
    ##to specify covariance matrix among all covariates
    R <- matrix(c(rep(a^2, m^2)), ncol = m)
    xcol=rep(a,m)
    R =cbind(xcol,R)
    xrow=rep(a,m+1)
    R =rbind(xrow,R)
    diag(R) <- a^2+1
    R[1,1] = 1
    S=R
    beta = c(r, b)
    id.iter = NA
    id.study  = NA

    ## Scale parameter (the smaller lambda, the greater becomes T)
    lambda <- exp(-6)#0.000001#1.7#eta.t <- -6

    ## Shape parameter
    nue <- 2

    ## Sample size
    n <- n

    ## Number of predictors
    p <- length(beta)


    ## Generate column vector of coefficients
    beta <- matrix(beta, ncol = 1)

    ## Generate correlated covariate vectors using a multivariate normal
    mu <- rep(0, p)
    X <- scale(mvrnorm(n, mu, S))


    ## Calculate survival times
    TT <- (-log(runif(n)) / (lambda * exp(X %*% beta)))^(1/nue)
    CT<- runif(n,0,c)
    event<-ifelse(TT<CT,1,0)
    T<-ifelse(event==1,TT,CT)
    cr<-1-sum(event)/n

    event[T>=ll]=1
    T[T>=ll]=ll
    llt=sum(T>=ll)/n


    M <- scale(matrix(rnorm(n*(ppp-m)), nrow=n, ncol=ppp-m))
    M=cbind(X[,2:(m+1)],M)
    MT=X[,2:(m+1)]
    colnames(M)=paste0('V', 1:ncol(M))
    X=X[,1]
    dat <- data.frame(T = T, X, M,event = event)

    ## cen =30% (0.30) of all marriages are getting divorced, i.e. 70% of all
    ## observations are censored ("event = rbinom(n, 1, 0.30)")
    ##dat <- data.frame(T = T, X, event = rbinom(n, 1, cen))
    ##dat$event <- ifelse(dat$T >= cen*100, 0, dat$event)
    ##dat$T <- ifelse(dat$T >= cen*100, cen*100, dat$T)

    dat$id.iter <- id.iter <- seq(1,n)
    dat$id.study  <- id.study
    ## Returning a matrix speeds-up things a lot... lesson learned.
    dat <- as.matrix(dat)
    return(list(dat=dat,MT=MT,llt=llt,cr=cr))
  }
}

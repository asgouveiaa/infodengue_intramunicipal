Rt<-function(obj, count = "casos", gtdist, meangt, sdgt, CI = "beta", alpha = .95, a0 = 2 , b0 = 3){
  
  
  if(!any(c(count,"SE") %in% names(obj))) stop("obj must be a data.frame with variables 
                                              SE and var, at least. Consider using getCases")
  y <- obj[,count]
  le <- length(y)
  if (le < 2*meangt) warning("you need a time series                           
                             with size at least 2 generation intervals to estimate Rt")
  
  if (gtdist == "normal") ga <- rev(dnorm(x = 1:le, mean = meangt, sd = sdgt))
  if (gtdist == "delta")  {
    ga <- rep(0, le)
    ga[le - meangt] <- 1
  }
  
  obj$Rt <- NA
  obj$lwr <- NA
  obj$upr <- NA
  obj$p1 <- NA
  
  for (t in ceiling(2*meangt):le){
    num = y[t]
    deno = sum(y[1:t] * ga[(le-t+1):le]) # equation 4.1 in Wallinga and Lipsitch 2007
    obj$Rt[t]<-num/deno
    if (CI == "beta"){
      obj$p1[t] <- 1 - pbeta(.5, shape1 = num, shape2 = deno)
      obj[t, c("lwr","upr")] <- ll(betaconf(alpha = alpha, x = num, 
                                            n = num + deno, a = a0, b = b0 ))
    }
  }
  
  obj
}


## obtain 100\lapha confidence/credibility intervals for the success probability \theta
betaconf <- function(alpha = .95, x, n, a = 1, b = 1, CP = "FALSE"){
  if(CP=="TRUE"){  
    lower <- 1 - qbeta((1-alpha)/2, n + x - 1, x)
    upper <- 1 - qbeta((1+alpha)/2, n - x, x + 1)
  }else{
    lower <- qbeta( (1-alpha)/2, a + x, b + n - x)
    upper <- qbeta(( 1+alpha)/2, a + x, b + n - x)  
  } 
  return(c(lower, upper))
  #CP stands for Clopper-Pearson
  #Default is 'Bayesian' with an uniform prior over p (a=b=1)
}


# R = theta/(1-theta)
ll <- function(x) {x/(1-x)}
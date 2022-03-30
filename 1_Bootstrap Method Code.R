# This is the simulation of bootstrap method.
# ci.pb is the confidence interval function

#zero-inflated gamma distribution data generator
rzigamma=function(n,delta,shape,scale){
  #n is the sample size.
  #delta is the proportion of 0 in binomial distribution.
  ni0=rbinom(1,n,prob=delta) #the prob of zeros which are generated using this way.
  ni1=n-ni0
  bdata=array(0,dim = ni0)
  gdata=rgamma(n=ni1,shape=shape,scale=scale)
  zigdata=c(bdata,gdata)
  return(list(ni0=ni0,ni1=ni1,zigdata=zigdata,zigdata1=gdata))
}

gam.mles= function(x){
  n1= length(x)
  xb = mean(x)
  xd = mean(log(x))
  s = log(xb)-xd
  a0 = (3.0-s+sqrt((s-3.0)^2+24.0*s))/12.0/s
  l = 1
  repeat{
    ans = (log(a0)-digamma(a0)-s)
    a1 = a0-ans/(1.0/a0-trigamma(a0))
    if(abs(ans) <= 1.0e-7 | l >= 30){break}
    a0 = a1
    l = l+1}
  ah = a1; bh = xb/a1
  return(c(ah,bh))
}
a.pivot = function(nr, x, cl){
  pivot = seq(1:nr)
  al = (1-cl)/2; ll = floor(nr*al);
  n1 = length(x)
  ml =  gam.mles(x)
  ah = ml[1]; bh = ml[2]
  for(i in 1:nr){
    x = rgamma(n1, shape = ah, scale = bh)
    xbs = mean(x); xds = mean(log(x))
    ml = gam.mles(x)
    ahs = ml[1]; bhs = ml[2]
    pivot[i] = ahs
  }
  return(pivot)
}

ci_.pb<-function(nr, data, cl){
  zigamma1 = data
  n0 = zigamma1$ni0; n1 = zigamma1$ni1
  n = n0+n1
  pb=a.pivot (nr, zigamma1$zigdata1, cl)
  zw=rnorm(nr)
  Tw<-(n0+zw^2/2)/(n+zw^2)-zw*sqrt(n0*(1-n0/n)+zw^2/4)/(n+zw^2)
  M=1/sqrt((1-Tw)*pb)
  return(M)
}
# the difference of CVs
ci.pb = function(zigamma1, zigamma2, nr, alpha){
  cl = 1-alpha
  CI1 = ci_.pb(nr, zigamma1, cl)
  CI2 = ci_.pb(nr, zigamma2, cl)
  M = sort(CI1-CI2)
  l = M[nr*(1-cl)/2]; u = M[nr*(1-(1-cl)/2)]
  return(list(lower = l, upper = u))
}



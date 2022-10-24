## extended generalised Pareto negative log-likelihood functions

## model 1 ##

.zidegpd1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  zidegpd1d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]],likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.zidegpd1.d12 <- function(pars, likdata) {
  zidegpd1d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]],likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}


.iG1 <- function(v, kappa) v^(1/kappa) # Quantile func of G(v)= v^Kappa

##Quantile func of DEGDP1 
#.iG1<- function(p, sigma, xi, kappa) {
#  xi <- sign(xi) * pmax(abs(xi), 1e-06)
#  out <- 1 - (p)^(1/kappa)
#  q<-ceiling(sigma * (out^(-xi) - 1)/xi) -1
#  q[q < 0] <- 0
#  return(q)
#}


.zidegpd1fns <- list(d0=.zidegpd1.d0, d120=.zidegpd1.d12, d340=NULL, m=1, iG=.iG1)


## model 2 ##

.G2 <- function(v, kappa1, kappa2, p) {
  p * v^kappa1 + (1 - p) * v^kappa2
}

.iG2i <- function(v, kappa1, kappa2, p) {
  vv <- range(c(v^1/kappa1, v^1/kappa2))
  d <- diff(vv)
  lo <- vv[1]
  while(.G2(lo, kappa1, kappa2, p) - v > 0) lo <- max(0, lo - d)
  hi <- vv[2]
  while(.G2(hi, kappa1, kappa2, p) - v > 0) hi <- min(1, hi + d)
  uniroot(function(x) .G2(x, kappa1, kappa2, p) - v, c(lo, hi))$root
}

.iG2 <- function(v, kappa1, kappa2, p) {
  n <- length(v)
  vapply(seq_len(n), function(i) .iG2i(v[i], kappa1[i], kappa2[i], p[i]), double(1))
}

.zidegpd2.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  zidegpd2d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$X[[6]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.zidegpd2.d12 <- function(pars, likdata) {
  zidegpd2d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$X[[6]] ,likdata$y[,1], likdata$dupid, likdata$duplicate)
}


.zidegpd2fns <- list(d0=.zidegpd2.d0, d120=.zidegpd2.d12, d340=NULL, m=2, iG=.iG2)


## model 3 ##

.zidegpd3.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  zidegpd3d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]],likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.zidegpd3.d12 <- function(pars, likdata) {
  zidegpd3d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]],likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.iG3 <- function(v, delta) 1 - qbeta(1 -v, 1/delta, 2)^(1/delta)

##Quantile func of DEGDP2
#.iG2<- function(p, delta,  sigma, xi){
#  q<- ceiling(qgpd(.G2(p, delta), sigma=sigma, xi=xi))-1
#  q[q < 0] <- 0
#  return(q)
#  
#}

####################################
.zidegpd3fns <- list(d0=.zidegpd3.d0, d120=.zidegpd3.d12, d340=NULL, m=3, iG=.iG3)

## model 4 ##

.zidegpd4.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  zidegpd4d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]],likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.zidegpd4.d12 <- function(pars, likdata) {
  zidegpd4d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]],likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}


.iG4 <- function(v, delta, kappa) 1 - qbeta(1 - v^(2/kappa), 1/delta, 2)^(1/delta)

##Quantile func of DEGDP3
#.iG3<- function(p, delta,kappa,  sigma, xi){
#  q<- ceiling(qgpd(.G3(p, delta, kappa), sigma=sigma, xi=xi))-1
#  q[q < 0] <- 0
#  return(q)
#  
#}


.zidegpd4fns <- list(d0=.zidegpd4.d0, d120=.zidegpd4.d12, d340=NULL, m=4, iG=.iG4)


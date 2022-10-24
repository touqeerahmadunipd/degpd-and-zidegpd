## Discrete extended generalised Pareto negative log-likelihood functions

## model 1 ##

.degpd1.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  degpd1d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.degpd1.d12 <- function(pars, likdata) {
  degpd1d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
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


.degpd1fns <- list(d0=.degpd1.d0, d120=.degpd1.d12, d340=NULL, m=1, iG=.iG1)


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

.degpd2.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  degpd2d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.degpd2.d12 <- function(pars, likdata) {
  degpd2d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}


.degpd2fns <- list(d0=.degpd2.d0, d120=.degpd2.d12, d340=NULL, m=2, iG=.iG2)


## model 3 ##

.degpd3.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  degpd3d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.degpd3.d12 <- function(pars, likdata) {
  degpd3d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}


.iG3 <- function(v, delta) 1 - qbeta(1 -v, 1/delta, 2)^(1/delta)

.degpd3fns <- list(d0=.degpd3.d0, d120=.degpd3.d12, d340=NULL, m=3, iG=.iG3)

##Quantile func of DEGDP2
#.iG2<- function(p, delta,  sigma, xi){
#  q<- ceiling(qgpd(.G2(p, delta), sigma=sigma, xi=xi))-1
#  q[q < 0] <- 0
#  return(q)
#  
#}

####################################
.degpd3fns <- list(d0=.degpd3.d0, d120=.degpd3.d12, d340=NULL, m=3, iG=.iG3)

## model 4 ##

.degpd4.d0 <- function(pars, likdata) {
  if (likdata$censored)
    stop("Censored likelihoods not currently available for extended GPDs.")
  degpd4d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.degpd4.d12 <- function(pars, likdata) {
  degpd4d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}


.iG4 <- function(v, delta, kappa) 1 - qbeta(1 - v^(2/kappa), 1/delta, 2)^(1/delta)

##Quantile func of DEGDP3
#.iG3<- function(p, delta,kappa,  sigma, xi){
#  q<- ceiling(qgpd(.G3(p, delta, kappa), sigma=sigma, xi=xi))-1
#  q[q < 0] <- 0
#  return(q)
#  
#}


.degpd4fns <- list(d0=.degpd4.d0, d120=.degpd4.d12, d340=NULL, m=4, iG=.iG4)


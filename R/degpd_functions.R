
#' Fit degpd and zidegpd

p.G <- function(u, type = 1, prob, kappa, delta) {
  if (!type %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  if (type %in% c(1,2, 3, 5,6) && missing(kappa)) {
    stop("Argument `kappa' missing.")
  }
  if (type %in% c(4, 5, 6) && missing(delta)) {
    stop("Argument `delta' missing.")
  }
  if (type == 6 && missing(prob)) {
    stop("Argument `prob' missing.")
  }
  if (type == 0) {
    return(u)
  } else if (type == 1) {
    return(u^kappa)
  } else if (type == 2) {
    F.min <- 1-pnorm(1, mean = 0, sd = 1/sqrt(kappa)) #Phi(sqrt(kappa)*(a-u)--->Phi(sqrt(kappa)*(0-1)= 1-Phi(sqrt(kappa)) cdf of standard normal dist, denonirator Eq.7 A flexible extended generalized Pareto distribution for tail estimation
    F.max <- pnorm(1, mean = 1, sd = 1/sqrt(kappa))#Phi(sqrt(kappa)*(b-u)--->Phi(sqrt(kappa)*(1-1)=1/2
    p <- (pnorm(u, mean = 1, sd = 1/sqrt(kappa)) - F.min)/(F.max -F.min)
    #p <- (pnorm(q, mean = 1, sd = 1/sqrt(kappa)) - F.min)/(F.max -F.min)
    return(p)
  }else if (type == 3) {
    lower=1/32
    upper=1/2
    stopifnot(lower <= upper,kappa > 0)
    aa <- rep(lower, length(u))
    bb <- rep(upper, length(u))
    normalize.factor <- pbeta(bb, kappa, kappa) - pbeta(aa, kappa, kappa)
    #tt <- pbeta(apply(cbind(apply(cbind((upper-lower)*u +lower , bb), 1, min), aa), 1, max), kappa, kappa)
    tt <- pbeta((upper-lower)*u +lower, kappa, kappa)
    tt <- tt / normalize.factor
    return(tt)
  }else if (type == 4) {
    return(1 - pbeta((1 - u)^delta, 1/delta, 2))
  } else if (type == 5) {
    return((1 - pbeta((1 - u)^delta, 1/delta, 2))^(kappa/2))
  } else if (type == 6) {
    return(prob * u^kappa + (1 - prob) * u^delta)
  }
}

#p.G(runif(10), type = 5, delta = 2, kappa = 1, prob = 0.1)

d.G <- function(u, type = 1, prob = NA, kappa = NA, delta = NA, log = FALSE) {
  if (!type %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  if (type %in% c(1,2, 3, 5,6) && missing(kappa)) {
    stop("Argument `kappa' missing.")
  }
  if (type %in% c(4, 5, 6) && missing(delta)) {
    stop("Argument `delta' missing.")
  }
  if (type == 6 && missing(prob)) {
    stop("Argument `prob' missing.")
  }
  if (log == FALSE) {
    if (type == 0) {
      return(1)
    } else if (type == 1) {
      return(kappa * u^(kappa - 1))
    } else if (type == 2) {
      F.min <-1- pnorm(1, mean = 0, sd = 1/sqrt(kappa))
      F.max <- pnorm(1, mean = 1, sd = 1/sqrt(kappa))
      d<- dnorm(u,mean = 1, sd = 1/sqrt(kappa))
      den <- (sqrt(kappa)*d)/(F.max -F.min)
      return(den)
    }else if (type == 3) {
      lower=1/32
      upper=1/2
      stopifnot(lower <= upper,kappa > 0)
      tt <- rep(0, length(u))
      normalize.factor <- pbeta(upper, kappa, kappa) - pbeta(lower, kappa, kappa)
      tt[u >= lower & u <= upper] <- dbeta(u[u >= lower & u <= upper],
                                           kappa, kappa) / normalize.factor
      return(tt)
    }
    else if (type == 4) {
      return(dbeta((1 - u)^delta, 1/delta, 2) * delta * (1 - u)^(delta - 1))
    } else if (type == 5) {
      return((kappa/2) * (1 - pbeta((1 - u)^delta, 1/delta, 2))^(kappa/2 - 1) * dbeta((1 - u)^delta, 1/delta, 2) * delta * (1 -
                                                                                                                              u)^(delta - 1))
    } else if (type == 6) {
      return(prob * kappa * u^(kappa - 1) + (1 - prob) * delta * u^(delta - 1))
    }
  } else {
    if (type == 0) {
      return(0)
    } else if (type == 1) {
      return(log(kappa) + (kappa - 1) * log(u))
    }else if (type == 2) {
      F.min <-1- pnorm(1, mean = 0, sd = 1/sqrt(kappa), log.p = TRUE)
      F.max <- pnorm(1, mean = 1, sd = 1/sqrt(kappa), log.p = TRUE)
      d<- dnorm(u,mean = 1, sd = 1/sqrt(kappa), log.p = TRUE)
      den <- (sqrt(kappa)*d)/(F.max -F.min)
      return(den)
    }else if (type == 3) {
      lower=1/32
      upper=1/2
      stopifnot(lower <= upper,kappa > 0)
      tt <- rep(0, length(u))
      normalize.factor <- pbeta(upper, kappa, kappa, log = TRUE) - pbeta(lower, kappa, kappa, log = TRUE)
      tt[u >= lower & u <= upper] <- dbeta(u[u >= lower & u <= upper],
                                           kappa, kappa,log = TRUE) / normalize.factor
      return(tt)
    }else if (type == 4) {
      return(dbeta((1 - u)^delta, 1/delta, 2, log = TRUE) + log(delta) + (delta - 1) * log(1 - u))
    } else if (type == 5) {
      return(log(kappa/2) + (kappa/2 - 1) * log(1 - pbeta((1 - u)^delta, 1/delta, 2)) + dbeta((1 - u)^delta, 1/delta, 2, log = TRUE) +
               log(delta) + (delta - 1) * log(1 - u))
    } else if (type == 6) {
      return(log(prob * kappa * u^(kappa - 1) + (1 - prob) * delta * u^(delta - 1)))
    }
  }
}
# x<-runif(10)
# d.G(x,  type=2,kappa=2)
##
q.G <- function(u, type = 1, prob = NA, kappa = NA, delta = NA) {
  if (!type %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  if (type %in% c(1,2, 3, 5,6) && missing(kappa)) {
    stop("Argument `kappa' missing.")
  }
  if (type %in% c(4, 5, 6) && missing(delta)) {
    stop("Argument `delta' missing.")
  }
  if (type == 6 && missing(prob)) {
    stop("Argument `prob' missing.")
  }
  if (type == 0) {
    return(u)
  } else if (type == 1) {
    return(u^(1/kappa))
  } else if (type == 2) {
    F.min <- 1-pnorm(1, mean = 0, sd = 1/sqrt(kappa))
    F.max <- pnorm(1, mean = 1, sd = 1/sqrt(kappa))
    q <- (qnorm(u* (F.max - F.min) + (F.min), mean = 0, 
                sd = 1))/sqrt(kappa) +1
    return(q)
  }else if (type == 3) {
    lower=1/32
    upper=1/2
    tt <- u
    pin <- pbeta(lower, kappa, kappa) + u * (pbeta(upper, kappa, kappa) - pbeta(lower, kappa, kappa))
    tt <- (qbeta(pin, kappa, kappa)-lower)/(upper-lower)
    return(tt)
  }else if (type == 4) {
    return(1 - qbeta(1 - u, 1/delta, 2)^(1/delta))
  } else if (type == 5) {
    return(1 - qbeta(1 - u^(2/kappa), 1/delta, 2)^(1/delta))
  } else if (type == 6) {
    dummy.func <- function(u, p, prob = NA, kappa = NA, delta = NA) {
      return(p.G(u = u, prob = prob, kappa = kappa, delta = delta, type = 4) - p)
    }
    find.root <- function(u, prob = NA, kappa = NA, delta = NA) {
      return(uniroot(dummy.func, interval = c(0, 1), p = u, prob = prob, kappa = kappa, delta = delta)$root)
    }
    return(sapply(u, FUN = find.root, prob = prob, kappa = kappa, delta = delta))
  }
}
#q.G(runif(10), kappa = 2, type = 4)
##
r.G <- function(n, prob = NA, kappa = NA, delta = NA, type = 1, unifsamp = NULL, direct = FALSE) {
  if (!type %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  if (type %in% c(1,2, 3, 5,6) && missing(kappa)) {
    stop("Argument `kappa' missing.")
  }
  if (type %in% c(4, 5, 6) && missing(delta)) {
    stop("Argument `delta' missing.")
  }
  if (type == 6 && missing(prob)) {
    stop("Argument `prob' missing.")
  }
  if (is.null(unifsamp)) {
    unifsamp <- runif(n)
  }
  if (type != 6 | (type == 6 & direct)) {
    return(q.G(unifsamp, prob = prob, kappa = kappa, delta = delta, type = type))
  } else if (type == 6 & !direct) {
    components <- sample(x = c(1, 2), size = n, replace = TRUE, prob = c(prob, 1 - prob))
    res <- c()
    res[components == 1] <- q.G(unifsamp[components == 1], prob = NA, kappa = kappa, delta = NA, type = 1)
    res[components == 2] <- q.G(unifsamp[components == 2], prob = NA, kappa = delta, delta = NA, type = 1)
    return(res)
  }
}

#r.G(10, type = 6, delta = 1, kappa = 1, prob = 0.1)
############################################
##    EXTENDED GPD TYPE 1 to 4           ##
############################################

pegpd <- function(q, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  return(p.G(extraDistr::pgpd(q, sigma = sigma, xi = xi), prob = prob, kappa = kappa, delta = delta, type = type))
}

pegpd(runif(10), type = 3, kappa = 1, sigma = 1, xi=0.1)
##
degpd <- function(x, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, log = FALSE) {
  if (log == FALSE) {
    return(d.G(extraDistr::pgpd(x, sigma = sigma, xi = xi), prob = prob, kappa = kappa, delta = delta, type = type) * extraDistr::pgpd(x, sigma = sigma,
                                                                                                                                       xi = xi))
  } else {
    return(d.G(extraDistr::pgpd(x, sigma = sigma, xi = xi), prob = prob, kappa = kappa, delta = delta, type = type, log = TRUE) + extraDistr::dgpd(x,
                                                                                                                                                   sigma= sigma, xi = xi, log = TRUE))
  }
}

degpd(runif(10), type = 2, kappa = 1, sigma = 1, xi=0.3)
#
qegpd <- function(p, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  return(extraDistr::qgpd(q.G(p, prob = prob, kappa = kappa, delta = delta, type = type), sigma = sigma, xi = xi))
}

qegpd(runif(10), type = 3, kappa = 1, sigma = 1, xi=0.3)


regpd <- function(n, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL, censoring = c(0, Inf)) {
  return(extraDistr::qgpd(r.G(n, prob = prob, kappa = kappa, delta = delta, type = type, unifsamp), sigma = sigma,
                          xi = xi))
}

regpd(10, type = 2, kappa = 1, sigma = 1, xi=0.3)

############################################
##   Discrete EXTENDED GPD TYPE 1 to 4    ##
############################################
##################Density function##########

ddiscegpd <- function(x, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  return(pegpd(x+1,prob = prob , kappa =kappa, delta = delta, sigma = sigma, xi = xi, type = type)- pegpd(x,prob = prob , kappa =kappa, delta = delta, sigma = sigma, xi = xi, type = type))
}
ddiscegpd(runif(10),type = 3, kappa = 1, sigma = 1, xi=0.3)
##################CDF################################
pdiscegpd <- function(q, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  return(pegpd(q+1,prob = prob , kappa =kappa, delta = delta, sigma = sigma, xi = xi, type = type))
}
################Quantile function###############################
qdiscegpd <- function(p, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  q<-ceiling(qegpd(p, prob = prob, kappa = kappa, delta = delta, type = type, sigma = sigma, xi=xi))-1
  q[q < 0] <- 0
  if(is.matrix(p)) {
    q <- matrix(q, ncol = ncol(p), nrow = nrow(p))
    colnames(q) <- colnames(p)
    rownames(q) <- rownames(p)
  }
  return(q)
}


rdiscegpd <- function(n, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL) {
  return(floor(regpd(n, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi,type = type, unifsamp)))
}

# x<-rdiscegpd(10000, type = 3, kappa = 10, sigma = 1, xi=0.3)
#################################################################################################
################Zero Inflated Discrete EGPD  ###################################################
################################################################################################
#--------------------------------------------------------------------------------------------------
dzidiscegpd<-function(x,  pi = NA,prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1)
{
  if (type==1) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==2) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==3) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==4) {
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
  }
  if (type==5) {
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==6) {
    if (any(prob <= 0) | any(prob >= 1) )  stop(paste("prob must be between 0 and 1", "\n", ""))
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (any(pi <= 0) | any(pi >= 1) )  stop(paste("pi must be between 0 and 1", "\n", ""))
  #if (any(prob <= 0) | any(prob >= 1) )  stop(paste("prob must be between 0 and 1", "\n", ""))
  #if (any(kappa <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  #if (any(delta <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", ""))
  if (any(xi <= 0) )  stop(paste("xi must be greater than 0", "\n", ""))
  if (any(x < 0) )  stop(paste("x must be 0 or greater than 0", "\n", ""))
  ly <- max(length(x),length(pi),length(prob),length(kappa),length(delta),length(sigma),length(xi))
  x <- rep(x, length = ly)
  pi <- rep(pi, length = ly)
  prob <- rep(prob, length = ly)
  kappa <- rep(kappa, length = ly)
  delta <- rep(delta, length = ly)
  sigma <- rep(sigma, length = ly)
  xi <- rep(xi, length = ly)
  fy <- rep(0, ly)
  fy <- ifelse((x==0), pi+(1-pi)*ddiscegpd(0, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type), (1-pi)*ddiscegpd(x, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type))
  fy
}

dzidiscegpd(x, type = 3, pi=0.1, kappa = 1, sigma = 1, xi=0.3)
##################################
# CDF ---------------------------------------------------------------------------------------------
pzidiscegpd<-function(q, pi = NA,prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1)
{
  if (type==1) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==2) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==3) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==4) {
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
  }
  if (type==5) {
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==6) {
    if (any(prob <= 0) | any(prob >= 1) )  stop(paste("prob must be between 0 and 1", "\n", ""))
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (any(pi <= 0) | any(pi >= 1) )  stop(paste("pi must be between 0 and 1", "\n", ""))
  #if (any(prob <= 0) | any(prob >= 1) )  stop(paste("prob must be between 0 and 1", "\n", ""))
  #if (any(kappa <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  #if (any(delta <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", ""))
  if (any(xi <= 0) )  stop(paste("xi must be greater than 0", "\n", ""))
  
  if (any(q < -1) )  stop(paste("y must be 0 or greater than 0", "\n", ""))
  ly <- max(length(q),length(pi),length(prob),length(kappa),length(delta),length(sigma),length(xi))
  q <- rep(q, length = ly)
  pi <- rep(pi, length = ly)
  prob <- rep(prob, length = ly)
  kappa <- rep(kappa, length = ly)
  delta <- rep(delta, length = ly)
  sigma <- rep(sigma, length = ly)
  xi <- rep(xi, length = ly)
  cdf <- rep(0,ly)
  cdf <- pdiscegpd(q, prob = prob, kappa = kappa, delta =delta, sigma = sigma, xi = xi, type = type)
  cdf <- ifelse(q < 0, 0, pi + (1 - pi) * cdf)
  cdf
}

################################################################################################
#Quantile function-----------------------------------------------------------------------------------------

qzidiscegpd <- function(p,  pi = NA,prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1)
{
  
  if (type==1) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==2) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==3) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==4) {
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
  }
  if (type==5) {
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==6) {
    if (any(prob <= 0) | any(prob >= 1) )  stop(paste("prob must be between 0 and 1", "\n", ""))
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (any(pi <= 0) | any(pi >= 1) )  stop(paste("pi must be between 0 and 1", "\n", ""))
  #if (any(prob <= 0) | any(prob >= 1) )  stop(paste("prob must be between 0 and 1", "\n", ""))
  #if (any(kappa <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  #if (any(delta <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  if (any(sigma <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  if (any(xi <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  if (any(p <= 0) | any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))
  # if (log.p == TRUE) p <- exp(p)   else p <- p
  #if (lower.tail == TRUE)  p <- p  else p <- 1 - p
  ly <- max(length(p),length(pi),length(prob),length(kappa),length(delta),length(sigma),length(xi))
  p <- rep(p, length = ly)
  pi <- rep(pi, length = ly)
  prob <- rep(prob, length = ly)
  kappa <- rep(kappa, length = ly)
  delta <- rep(delta, length = ly)
  sigma <- rep(sigma, length = ly)
  xi <- rep(xi, length = ly)
  pnew <- ((p-pi)/(1-pi))-(1e-7)
  pnew <- ifelse(pnew > 0, pnew,0 )
  q<- qdiscegpd(pnew, prob = prob, kappa = kappa, delta =delta, sigma = sigma, xi = xi, type = type)
  q
}

###Zero Inflated discrete degpd ###
rzidiscegpd <- function(n,pi=NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL ) {
  z <- rbinom(n,size=1,prob=pi)
  y <- (1-z)*rdiscegpd(n, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type, unifsamp = NULL)
  return(y)
}

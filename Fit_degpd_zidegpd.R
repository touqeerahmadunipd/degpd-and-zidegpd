setwd("D:/DEVGAM/evgam-master/") #set working directory 

dest <- "./R/"      # this function all R function 
files = list.files(dest, full.names = T)
for (i in 1:length(files)) {
  source(files[i])
}

dest <- "./src/"  # This function will call all C++ files 
files = list.files(dest, full.names = T)
for (i in 1:length(files)) {
  Rcpp::sourceCpp(files[i])
}
################DEGPD1 
library(mev)
library(evd)
n1 <-20 # number of locations
n2 <- 500# number of realizations at each location
test_df <- data.frame(x = rep(rnorm(n2), n1))

test_df$y <- floor(rextgp(n1*n2,kappa = 10, sigma =1, xi=0.2, type =1 ))
plot(table(test_df$y))

######Polynomials------------------------------------------
#test_pars <- poly(test_df$x, 6) %*% sapply(c(15, 10, 10), function(x) runif(6, -x, x))
#test_df$y <- rdegpd(n1*n2,kappa = 10*exp(test_pars[, 1]), sigma =exp(test_pars[, 2]), xi=exp(test_pars[, 3])/5, type =1 )

inits1<- c(2.0896946, 0.1762738, 0.95)
fmla_degpd1 <- list(lsigma = y ~ s(x), lxi = ~1, lkappa = ~1)
degpd1<- devgam(fmla_degpd1, data = test_df, family = "degpd",degpd.args = list(m=1), trace = 2, inits = inits1) # To fit degpd 2,3 and 4, put m=2, 3 and 4. 
summary(degpd1)# To fit zidegpd, change family as "zidegpd" and change put m=1, 3,4 for zidegpd 1, 3, 4.  zidegpd 2 is not developed yet. 
est_par =predict(degpd1, type = 'response')
head(est_par)
plot(degpd1)
q<-predict(degpd1,  prob = c(.50,0.75,0.99))

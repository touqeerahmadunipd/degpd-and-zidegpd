**Fitting of degpd-and-zidegpd**

The "Fit_degpd_zidegpd.R" function will work to fit simply discrete extended generalized Pareto distribution (degpd) and zero-inflated discrete extended generalized Pareto distribution (zidegpd) with thier GAM forms as well. These models are proposed in the paper "Ahmad, T., Gaetan, C., & Naveau, P. (2022). Modelling of discrete extremes through extended versions of discrete generalized Pareto distribution. arXiv. https://doi.org/10.48550/arXiv.2210.15253." We are using the functions of evgam R package (Youngman, 2020): An R package for Generalized Additive Extreme Value Models. 
https://doi.org/10.48550/arXiv.2003.04067 behind to run our own developed R code.

The example with fitting of degpd 1 model is shown "Fit_degpd_zidegpd.R" code. The other degpd 2, 3 and 4 models can be fitted by changing "***m=1 is corresponding to degpd 1 in paper", "m=2 is corresponding*** $$G\left(u\right)={u}^{\kappa}$$", "m=3 is corresponding to degpd 3 in paper" and "m=4 is corresponding to degpd 3 in paper"***. The zidegpd models can also be fitted by changing family " **degpd" to "zidegpd**" and putting "**m=1 is corresponding to zidegpd 1 in paper", "m=3 is corresponding to zidegpd 2 in paper" and "m=4 is corresponding to degpd 3 in paper**". In addition zidegpd 4 is not developed yet.
$$\left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right)$$
**Note** 
1. Ignore the error
"Error in Rcpp::sourceCpp(files[i]) : 
  The filename 'Makevars' does not have an extension of .cc or .cpp so cannot be compiled." The code is working correctly, this error appearing from evgam package code, we are using some function behind from evgam. Surely, we will introduce our proposed models as new families in evgam.
  
2. We developed and executed this code on window operating system, you may face problem when you excute on MAC, we will try to check the code on MAC also, and provide the piece of code as soon as posible.    

**For further discussion feel free to write me on (touqeer.ahmad@phd.unipd.it)**
 


**_To fit discrete extended generalized Pareto distribution (degpd) and zero-inflated discrete extended generalized Pareto distribution (zidegpd), we developed a code for new families and run it using evgam package fucntion_**


**Fitting of degpd-and-zidegpd**

The "Fit_degpd_zidegpd.R" function will work to fit simply discrete extended generalized Pareto distribution (degpd) and zero-inflated discrete extended generalized Pareto distribution (zidegpd) with thier GAM forms as well. These models are proposed in the paper "Ahmad, T., Gaetan, C., & Naveau, P. (2023). A regression model for count data with extreme observations." We are using the functions of evgam R package (Youngman, 2020): An R package for Generalized Additive Extreme Value Models. 
https://doi.org/10.48550/arXiv.2003.04067 behind to run our own developed R code.

The example with fitting of degpd 1 model is shown "Fit_degpd_zidegpd.R" code. The other degpd 2, 3 and 4 models can be fitted by changing m. In the code, the m=1 is corresponding to model $$G\left(u; \psi\right)={u}^{\kappa},$$
m=2 is corresponding to model
$$G\left(u;\psi\right)= p{u}^{\kappa_1} + \left(1-p\right){u}^{\kappa_2},$$
m=3 is corresponding to model
$$G\left(u;\psi\right)=1-D_{\delta}\{\left(1-u\right)^{\delta}\},$$
and m=4 is corresponding to model
$$G\left(u;\psi\right)=\left[1-D_{\delta}\{(1-u)^{\delta}\}\right]^{\kappa/2}$$
The zidegpd models can also be fitted by changing family " **degpd**" to "**zidegpd**" and by changing **m**. The m=1 is corresponding to model $$G\left(u; \psi\right)={u}^{\kappa},$$
m=2 is corresponding to model (**not developed yet**)
$$G\left(u;\psi\right)= p{u}^{\kappa_1} + \left(1-p\right){u}^{\kappa_2},$$
m=3 is corresponding to model
$$G\left(u;\psi\right)=1-D_{\delta}\{\left(1-u\right)^{\delta}\},$$
and m=4 is corresponding to model
$$G\left(u;\psi\right)=\left[1-D_{\delta}\{(1-u)^{\delta}\}\right]^{\kappa/2}$$
**Note** 
1. Ignore the error
"Error in Rcpp::sourceCpp(files[i]) : 
  The filename 'Makevars' does not have an extension of .cc or .cpp so cannot be compiled." The code is working correctly, this error appearing from evgam package code, we are using some function behind from evgam. Surely, we will introduce our proposed models as new families in evgam.
  
2. We developed and executed this code on window operating system, you may face problem when you excute on MAC, we will try to check the code on MAC also, and provide the piece of code as soon as posible.    

**For further discussion feel free to write me on (touqeer.ahmad@ensai.fr)**
 


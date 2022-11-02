**Fitting of degpd-and-zidegpd**

The "Fit_degpd_zidegpd.R" function will work to fit simply discrete extended generalized Pareto distribution (degpd) and zero-inflated discrete extended generalized Pareto distribution (zidegpd) with thier GAM forms as well. These models are proposed in the paper "Ahmad, T., Gaetan, C., & Naveau, P. (2022). Modelling of discrete extremes through extended versions of discrete generalized Pareto distribution. arXiv. https://doi.org/10.48550/arXiv.2210.15253." We are using the functions of evgam R package (Youngman, 2020) behind to run our own developed R code.

The example with fitting of degpd 1 model is shown "Fit_degpd_zidegpd.R" code. The other degpd 2, 3 and 4 models can be fitted by changing "m=1" to "m=2", "m=3" and "m=4". The zidegpd models can also be fitted by changing family " degpd" to "zidegpd" and putting "m=1", "m=3" and "m=4". In addition zidegpd 2 is not developed yet.

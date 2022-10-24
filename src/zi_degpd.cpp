// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0;
// //' Zero Inflated Discrete Extended generalized Pareto distribution of type 1 (zideGPD1) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each zideGPD parameter
// //' @param X1 a design matrix for the zideGPD log scale parameter
// //' @param X2 a design matrix for the zideGPD shape parameter
// //' @param X3 a design matrix for the zideGPD log kappa
// //' @param X3 a design matrix for the zideGPD logit pi
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return zidegpd1d0 a scalar, the negative log-liklihood
// //' @return zidegpd1d12 a matrix, first then second derivatives w.r.t. eeGPD1 parameters
// //' @return zidegpd1d34 a matrix, third then fourth derivatives w.r.t. eeGPD1 parameters (Not given)
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]

double zidegpd1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, arma::vec yvec, const arma::uvec& dupid, int dcate)
{  
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec logitpivec = X4 * Rcpp::as<arma::vec>(pars[3]);
  int nobs = yvec.size();
   
    if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
    lkappavec = lkappavec.elem(dupid);
    logitpivec = logitpivec.elem(dupid);
  }
  
  double y, lsigma, lxi, lkappa, logitpi;
  double e1,e2,e3, e4, e5;
  double ee1,ee2,ee3, ee4, ee5;
  double hi; 
  double lo;
  double nllh=0.0;
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    lkappa = lkappavec[j];
    logitpi = logitpivec[j];
    if(y>0){ 
	e1=1/exp(lxi);
    e2 =  exp(lxi) / exp(lsigma);
	e3= 1+ (y+1)*e2;
	e4 = 1+  y*e2;
	e5= exp(lkappa);
    hi=  R_pow(1- R_pow(1/e3, e1) , e5);
    lo= R_pow(1-R_pow(1/e4, e1) , e5);
	nllh += -log((1/(1+exp(logitpi)))*(hi - lo)); 
	 
}else{
	 ee1=1/exp(lxi);
    ee2 =  exp(lxi) / exp(lsigma);
    ee3 = 1+  ee2;
    ee4= exp(lkappa);
    ee5= R_pow(1- R_pow(1/ee3, ee1) , ee4);
	nllh += -log(exp(logitpi)/(1+exp(logitpi))+ (1/(1+exp(logitpi)))*ee5);
}
   //if (!ISNA(nllh)){
    //nllh = 1e20;
    // break;
//}
 //Rprintf("hi %f ",hi); 
   // Rprintf("lo %f \n",lo);
  }
  return(nllh);
}


// //' @rdname ziegpd1d0
// [[Rcpp::export]]
arma::mat zidegpd1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3,  arma::mat X4,arma::vec yvec, const arma::uvec& dupid, int dcate )
{
  
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec logitpivec = X4 * Rcpp::as<arma::vec>(pars[3]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs,  14);
  
  if (dcate == 1) {
   lsigmavec = lsigmavec.elem(dupid);
   lxivec = lxivec.elem(dupid);
   lkappavec = lkappavec.elem(dupid);
   logitpivec = logitpivec.elem(dupid);
  }
  
  double y, lsigma, lxi, lkappa, logitpi;
  double ee1, ee2, ee3, ee4, ee5,ee7, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21,ee22, ee23, ee24,  ee28, ee29,ee30;
  double ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39, ee40;
  double ee41, ee42,ee43,  ee44, ee45, ee46, ee47;
  double ee49, ee51, ee53, ee55, ee56, ee57;
  double ee58, ee59, ee60, ee61, ee62;
  double ee63, ee64, ee65, ee67;
  double ee69, ee72,  ee74, ee75, ee76;
  double ee77, ee78,ee79, ee80;
  double ee84, ee85, ee86, ee91;
  double ee92,ee93, ee94, ee95, ee98,ee99;
  double ee102, ee103, ee104, ee105, ee108, ee109;
  double ee113, ee114,ee116, ee117, ee119,ee122, ee124, ee126;
  double eee1, eee2, eee3,eee4, eee5, eee6, eee7, eee8;
  double eee9, eee10;
  double eee11, eee12, eee13, eee14, eee15, eee16;
  double eee17, eee18,eee19, eee20,eee21,  eee22,eee23, eee24;
  double eee25,eee27,eee28, eee30, eee31,eee32;
  double eee33, eee35,eee36,eee37,eee38, eee40,  eee41;
  double eee42, eee43,  eee45, eee46, eee47,  eee49, eee50, eee51,  eee52, eee54, eee55; 
  
  for (int j=0; j < nobs; j++) {
    
    y  = yvec[j];
      lsigma  = lsigmavec[j];
      lxi  = lxivec[j];
      lkappa  = lkappavec[j];
      logitpi  = logitpivec[j];
   
    if(y>0){   
    ee1 = exp(lxi);
    ee2 = exp(lsigma);
    ee3 = 1/ee1;
    ee4 = 1 + y;
    ee5 = exp(lkappa);
    ee7 = ee4 * ee1/ee2;
    ee9 = y * ee1/ee2;
    ee10 = ee7 + 1;
    ee11 = 1 + ee9;
    ee12 = R_pow(ee10,ee3);
    ee13 = R_pow(ee11,ee3);
    ee14 = exp(logitpi);
    ee15 = 1 - 1/ee12;
    ee16 = 1 - 1/ee13;
    ee17 = R_pow(ee15,ee5);
    ee18 = R_pow(ee16,ee5);
    ee19 = ee5 - 1;
    ee20 = 1 + ee14;
    ee21 = 1 + ee3;
    ee22 = ee17 - ee18;
    ee23 = R_pow(ee15,ee19);
    ee24 = R_pow(ee16,ee19);
    ee28 = ee17 - (ee22 * ee14/ee20 + ee18);
    ee29 = R_pow(ee10,ee21);
    ee30 = R_pow(ee11,ee21);
    ee31 = ee3 - 1;
    ee32 = log1p(ee7);
    ee33 = log1p(ee9);
    ee34 = 2/ee1;
    ee35 = 1 - ee14/ee20;
    ee36 = log(ee15);
    ee37 = log(ee16);
    ee38 = ee23 * ee4;
    ee39 = y * ee24;
    ee40 = ee38/ee29;
    ee41 = ee17 * ee36;
    ee42 = ee18 * ee37;
    ee43 = ee39/ee30;
    ee44 = ee29 * ee2;
    ee45 = R_pow(ee10,ee31);
    ee46 = ee30 * ee2;
    ee47 = R_pow(ee11,ee31);
    ee49 = ee45 * ee4/ee2;
    ee51 = ee12 * ee32/ee1;
    ee53 = ee13 * ee33/ee1;
    ee55 = y * ee47/ee2;
    ee56 = ee49 - ee51;
    ee57 = R_pow(ee10,ee34);
    ee58 = R_pow(ee11,ee34);
    ee59 = ee55 - ee53;
    ee60 = ee28 * ee2;
    ee61 = ee12 * ee1;
    ee62 = ee13 * ee1;
    ee63 = ee28 * ee20;
    ee64 = ee41 - ee42;
    ee65 = ee43 - ee40;
    ee67 = ee56 * ee23/ee57;
    ee69 = ee24 * ee59/ee58;
    ee72 = ee32/ee61 - ee4/ee44;
    ee74 = ee33/ee62 - y/ee46;
    ee75 = ee5 - 2;
    ee76 = ee23 * ee72;
    ee77 = ee24 * ee74;
    ee78 = 2 * ee21;
    ee79 = ee67 - ee69;
    ee80 = ee77 - ee76;
    ee84 = R_pow(ee60,2);
    ee85 = ee67 - (ee79 * ee14/ee20 + ee69);
    ee86 = ee22 * ee35;
    ee91 = ee40 + ee14 * ee65/ee20 - ee43;
    ee92 = R_pow(ee15,ee75);
    ee93 = ee41 - (ee64 * ee14/ee20 + ee42);
    ee94 = R_pow(ee16,ee75);
    ee95 = R_pow(ee63,2);
    ee98 = ee12 * ee21 * ee4 * ee1;
    ee99 = ee94 * ee19;
    ee102 = y * ee21 * ee13 * ee1;
    ee103 = R_pow(ee44,2);
    ee104 = ee56 * ee92;
    ee105 = R_pow(ee61,2);
   // ee107 = ee86/ee28 - 1;
    ee108 = R_pow(ee46,2);
    ee109 = R_pow(ee62,2);
    ee113 = ee98/ee2 - ee29 * ee32/ee1;
    ee114 = R_pow(ee10,ee78);
    ee116 = ee23 * ee5 * ee36;
    ee117 = ee92 * ee19;
    ee119 = ee24 * ee5 * ee37;
    ee122 = R_pow(ee11,ee78);
   // ee123 = 1 + 3/ee1;
    ee124 = ee3 - ee78;
    ee126 = ee102/ee2 - ee30 * ee33/ee1;


out(j, 0) =  -(ee35 * ee5 * ee65/ee60); //# w.r.t. lsigma
    out(j, 1) = -(ee80 * 
                       ee35 * ee5/ee28); // # w.r.t. lxi 
    out(j, 2) = -(ee64 * ee35 * ee5/ee28);// # w.r.t. lkappa
    out(j, 3) = ee86 * ee14/ee63; // # logitpi
    
    out(j, 4) = -(((R_pow(y,2) * 
                         (ee24 * ee21 * R_pow(ee11,ee124) * ee1 - ee99/ee122) - (R_pow(ee10,ee124) * 
                                                                            ee23 * ee21 * ee1 - ee117/ee114) * R_pow(ee4,2))/(ee28 * R_pow(ee2,2)) - 
                        (ee60 - ee91 * ee5) * ee65/ee84) * ee35 * ee5); //# w.r.t (lsigma, lsigma)
    out(j, 5) = -((ee91 * 
                        ee80 * ee5/ee60 + y * (((ee46 - ee102)/ee108 + (ee47 * 
                                                                          ee1 * ee33/ee109 - 1/ee30)/ee2) * ee24 - ee99 * ee74/ee46) - 
                        (((ee44 - ee98)/ee103 + (ee45 * ee1 * ee32/ee105 - 1/ee29)/ee2) * 
                           ee23 - ee117 * ee72/ee44) * ee4) * ee35 * ee5/ee28); // # w.r.t (lsigma, xi)
    out(j, 6) =  -((ee91 * ee64 * ee5/ee28 + y * (ee119/ee30 + 
                                                     ee24/ee30) - (ee116/ee29 + ee23/ee29) * ee4) * ee35 * 
                      ee5/ee60);  // # w.r.t (lsigma, lkappa)
    out(j, 7) =  (ee91 * ee22 * ee20/ee95 + ee65/ee63) *ee35 * ee5 * ee14/ee2; // # w.r.t (l, logitpi)
    out(j, 8) = -(((ee94 * 
                        ee74 * ee59/ee58 - ee104 * ee72/ee57) * ee19 + ee24 * 
                       (y * (1/ee46 + ee2 * ee126/ee108) - (ee13 + ee55 - ee53) * 
                          ee1 * ee33/ee109) - (((ee113 * ee2/ee103 + 1/ee44) * 
                                                  ee4 - (ee49 + ee12 - ee51) * ee1 * ee32/ee105) * ee23 + 
                                                 ee85 * ee80 * ee5/ee28)) * ee35 * ee5/ee28); // # w.r.t (lxi, lxi)
    out(j, 9) = -((ee56 * 
                        (ee116/ee57 + ee23/ee57) - (ee85 * ee64 * ee5/ee28 + 
                                                      (ee119/ee58 + ee24/ee58) * ee59)) * ee35 * ee5/ee28); // # wer.t (lxi, lkappa)
      
    out(j, 10) = (ee79/ee63 - ee85 * ee22 * ee20/ee95) * ee35 * ee5 * ee14; //  # w.r.t (lxi, logitpi) 
    out(j, 11) = -(((ee17 * 
                          R_pow(ee36,2) - (ee93 * ee64/ee28 + ee18 *  R_pow(ee37,2))) * ee5 + 
                         ee41 - ee42) * ee35 * ee5/ee28); //  # w.r.t (lkappa, lkappa)
    out(j, 12) =  (ee64/ee63 - 
                       ee22 * ee93 * ee20/ee95) * ee35 * ee5 * ee14; //  # w.r.t (lkappa, logitpi)
    out(j, 13) =  ee22 * R_pow(ee35,2) * ee14/ee63; // # weret (logitpi, logitpi)
  }

else {
     	   
     eee1 = exp(lxi);
    eee2 = exp(lsigma);
    eee3 = eee1/eee2;
    eee4 = 1 + eee3;
    eee5 = 1/eee1;
    eee6 = exp(lkappa);
    eee7 = R_pow(eee4,eee5);
    eee8 = 1 - 1/eee7;
    eee9 = exp(logitpi);
    eee10 = R_pow(eee8,eee6);
    eee11 = eee10 + eee9;
    eee12 = 1 + eee5;
    eee13 = R_pow(eee4,eee12);
    eee14 = eee6 - 1;
    eee15 = R_pow(eee8,eee14);
    eee16 = log1p(eee3);
    eee17 = eee13 * eee2;
    eee18 = log(eee8);
    eee19 = eee11 * eee13;
    eee20 = 1 + eee9;
    eee21 = R_pow(eee4,(eee5 - 1));
    eee22 = eee7 * eee1;
    eee23 = R_pow(eee4,(2/eee1));
    eee24 = eee19 * eee2;
    eee25 = eee21/eee2;
    eee27 = eee7 * eee16/eee1;
    eee28 = 1/eee17;
    eee30 = eee15 * eee6;
    eee31 = eee25 - eee27;
    eee32 = 1 - eee11/eee20;
    eee33 = 2 * eee6;
    eee35 = eee16/eee22 - eee28;
    eee36 = R_pow(eee24,2);
    eee37 = R_pow(eee8,(eee33 - 1));
    eee38 = R_pow(eee8,(eee6 - 2));
    eee40 = eee12 * eee7 * eee1;
    eee41 = eee11 * eee23;
    eee42 = eee38 * eee14;
    eee43 = eee11 * eee2;
   // eee44 = R_pow(eee11,2);
    eee45 = R_pow(eee17,2);
    eee46 = R_pow(eee22,2);
    eee47 = eee32 * eee15;
    eee49 = R_pow(eee8,(2 * eee14)) * eee6;
    eee50 = eee15/eee13;
    eee51 = eee15/eee23;
    eee52 = eee10 * eee6;
    eee54 = eee40/eee2 - eee13 * eee16/eee1;
    eee55 = eee17 - eee40;

     
    //# first derivatives
    out(j, 0) =  eee30/eee24;  //  # weereetee lsigma
    out(j, 1) = eee30 * eee35/eee11;  //# w.r.t. lxi
    out(j, 2)  = -(eee52 * eee18/eee11); //# w.r.t. lkappa
    out(j, 3)  =   -(eee32 * eee9/eee11);// # w.r.t. logitpi
    
    //# second derivatives
    out(j, 4)  = -(((eee11 * eee55 - eee30) * 
                     eee15/eee36 + eee42/(eee11 * R_pow(eee4,(2 * eee12)) * R_pow(eee2,2))) * 
                    eee6);//  # weer.t (lsigma, lsigma)
    out(j, 5)  =  ((eee55/eee45 + (eee21 * eee1 * eee16/eee46 - 
                                  1/eee13)/eee2) * eee15 + (eee49/eee19 - eee42/eee13) * eee35/eee2) * eee6/eee11; //# w.r.t (lsigma, xi)
    out(j, 6)  = -(((eee37/eee19 - eee50) * eee6 * 
                     eee18 - eee50) * eee6/eee43);// # weereet (lsigma, lkappa)
      out(j, 7) = -((eee47/eee19 + 
                       eee15/(eee20 * eee13)) * eee6 * eee9/eee43);// # w.r.t (lsigma, logitpi)
     out(j, 8)  =  ((eee54 * eee2/eee45 + eee28 - 
                       (eee25 + eee7 - eee27) * eee1 * eee16/eee46) * eee15 + (eee42/eee23 - 
                                                                          eee49/eee41) * eee31 * eee35) * eee6/eee11; //# w.r.t (lxi, lxi)
      out(j, 9)  =  -(((eee51 - 
                          eee37/eee41) * eee6 * eee18 + eee51) * eee31 * eee6/eee11);//# w.r.t (lxi, lkappa)
      
    out(j, 10)  = (eee47/eee41 + eee15/(eee20 * eee23)) * eee31 * eee6 * eee9/eee11;//# w.r.t (lxi, logitpi)
    out(j, 11)  = -(((eee10 - R_pow(eee8,eee33)/eee11) * 
                      eee6 * eee18 + eee10) * eee6 * eee18/eee11);//   # w.r.t (lkappa, lkappa)
    out(j, 12)  =  (eee32 * 
                     eee10/eee11 + eee10/eee20) * eee6 * eee9 * eee18/eee11;// # w.r.t ( lkappa, logitpi)
    out(j, 13)  = -(eee32 * (1 - (1/eee11 + 
                                  1/eee20) * eee9) * eee9/eee11);//# w.r.t (logitpi, logitpi)                                                                          
  }
}
     
   return out;
}
    


// //' Discrete Extended generalized Pareto distribution of type 2 (zideGPD1) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each zideGPD parameter
// //' @param X1 a design matrix for the zideGPD log scale parameter
// //' @param X2 a design matrix for the zideGPD log shape parameter
// //' @param X3 a design matrix for the zideGPD log kappa1
// //' @param X4 a design matrix for the deGPD log kappa2 parameter
// //' @param X5 a design matrix for the zideGPD logit p parameter
// //' @param X6 a design matrix for the zideGPD logit pi kappa
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return zidegpd1d0 a scalar, the negative log-liklihood
// //' @return zidegpd1d12 a matrix, first then second derivatives w.r.t. zideGPD2 parameters
// //' @return zidegpd1d34 a matrix, third then fourth derivatives w.r.t. zideGPD1 parameters (Not given)
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double zidegpd2d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4 , const arma::mat& X5, const arma::mat& X6, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
  
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappa1vec = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec lkappa2vec = X4 * Rcpp::as<arma::vec>(pars[3]);
  arma::vec logitpvec = X5 * Rcpp::as<arma::vec>(pars[4]);
  arma::vec logitpivec = X5 * Rcpp::as<arma::vec>(pars[5]);
  int nobs = yvec.size();
  
   if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
     lkappa1vec = lkappa1vec.elem(dupid);
    lkappa2vec = lkappa2vec.elem(dupid);
    logitpvec = logitpvec.elem(dupid);
    logitpivec = logitpivec.elem(dupid);
  }
  
  
  double y, lsigma, lxi,lkappa1, lkappa2,logitp, logitpi;
  double e1,e2,e3, e4, e5,e6,e7,e8,e9,e10,e11,e12;
  double ee1,ee2,ee3, ee4, ee5, ee6,ee7, ee8;
  double hi; 
  double lo;
  double nllh=0.0;
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    lkappa1 = lkappa1vec[j];
    lkappa2 = lkappa2vec[j];
    logitp = logitpvec[j];
    logitpi = logitpivec[j];
    if(y>0){ 
	e1=1.0/exp(lxi);
    e2 =  exp(lxi) / exp(lsigma);
    e3= 1+ (y+1)*e2;
    e4= R_pow(1/e3, e1);
    e5= R_pow(1-e4, exp(lkappa1)); //(H(y+1))^lkappa1
    e6= 1+ y*e2;
    e7= R_pow(1/e6, e1);
    e8= R_pow(1-e7, exp(lkappa1)); //(H(y))^lkappa1
    e9= exp(logitp)/ (1+exp(logitp));
	hi= e9*(e5-e8);
	e10= R_pow(1-e4, exp(lkappa2)); //(H(y+1))^lkappa2
	e11= R_pow(1-e7, exp(lkappa2)); //(H(y))^lkappa2
	e12= 1/(1+exp(logitp));
	lo= e12*(e10-e11);
    nllh += -log((1/(1+exp(logitpi)))*(hi + lo));
}else{
	ee1=1.0/exp(lxi);
    ee2 =  exp(lxi) / exp(lsigma);
    ee3= 1+ ee2;
    ee4= R_pow(1/ee3, ee1);
    ee5= R_pow(1-ee4, exp(lkappa1)); //(H(1))^lkappa1
    ee6= R_pow(1-ee4, exp(lkappa2)); //(H(1))^lkappa2
    ee7= exp(logitp)/ (1+exp(logitp));
	hi= ee7*ee5;
	ee8= 1/(1+exp(logitp));
	lo= ee8*ee6;
	nllh += -log(exp(logitpi)/(1+exp(logitpi))+ (1/(1+exp(logitpi)))*(hi+lo));
	}
   //if (!ISNA(nllh)){
    //nllh = 1e20;
    // break;
//}
 //Rprintf("hi %f ",hi); 
   // Rprintf("lo %f \n",lo);
  }
  return(nllh);
}


// //' @rdname zidegpd1d0
// [[Rcpp::export]]
arma::mat zidegpd2d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3,  arma::mat X4, arma::mat X5,arma::mat X6, arma::vec yvec, const arma::uvec dupid, int dcate )
{
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappa1vec = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec lkappa2vec = X4 * Rcpp::as<arma::vec>(pars[3]);
  arma::vec logitpvec = X5 * Rcpp::as<arma::vec>(pars[4]);
  arma::vec logitpivec = X5 * Rcpp::as<arma::vec>(pars[5]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs,  27);
  
    if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);   
	lkappa1vec = lkappa1vec.elem(dupid);
    lkappa2vec = lkappa2vec.elem(dupid);
    logitpvec = logitpvec.elem(dupid);
    logitpivec = logitpivec.elem(dupid);
  }
  
  
  double y, lsigma, lxi,lkappa1, lkappa2,logitp, logitpi;
  double ee1, ee2, ee3, ee4, ee6, ee8,  ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21,ee22, ee23, ee24, ee25,ee26, ee27,    ee28,ee29, ee30;
  double ee31, ee32, ee33, ee34,  ee136, ee37, ee39, ee40;
  double ee41,  ee42,ee43,ee44,ee45,ee46,ee47,ee48;
  double ee49, ee50,ee51,ee52,ee53,ee55,   ee57;
  double  ee59,  ee61,ee62,  ee63, ee64 ;
  double  ee65, ee66,  ee67, ee68 ;
  double ee70, ee72 , ee74 , ee75, ee76;
  double ee77, ee78,ee79,ee80, ee81,ee82, ee83;
  double ee84, ee85, ee86, ee87,ee88, ee90;
  double ee92, ee93 , ee94,ee97,ee98, ee99,ee100,ee101;
  double  ee102, ee103 , ee104,ee107, ee108,ee109,ee110, ee111,  ee112, ee113 ,ee115;
  double ee116, ee117,  ee119,ee121,ee123,ee124,ee125;
  double ee129, ee130, ee137 ,ee138,   ee141;
  double ee142, ee143, ee144,ee145, ee146,ee150;
  double ee151, ee152, ee153, ee154,ee155,ee158, ee160, ee162,ee163;
  double ee166,ee167, ee168, ee169,  ee171, ee175,ee176,ee177,  ee179;
  double  ee184;
  double ee186, ee187,ee189,ee190, ee192,ee194,ee195, ee197,ee198,ee200, ee203,ee204,ee205,ee207,ee208;

  double eee1, eee2, eee3,eee4, eee5, eee6, eee7, eee8;
  double eee9, eee10;
  double eee11, eee12, eee13, eee14,  eee15, eee16;
  double eee17, eee18, eee19,eee20,eee21,  eee22,eee23, eee24;
  double eee25,eee26,eee27,eee28,eee29, eee30, eee32;
  double eee33,eee34, eee35,eee36,eee37,eee38, eee39,eee40;
  double eee42, eee43, eee44, eee45,  eee47,  eee48, eee50,eee51,  eee52,eee54, eee55, eee56; 
  double eee57,eee58,  eee59, eee60, eee62, eee63, eee64,  eee68, eee70;
  double eee71, eee72,  eee73,eee74, eee76, eee80, eee86;
  double eee89,eee90,eee92,  eee94;
  
  
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    lkappa1 = lkappa1vec[j];
    lkappa2 = lkappa2vec[j];
    logitp = logitpvec[j];
    logitpi = logitpivec[j];
   
    if(y>0){   
    ee1 = exp(lxi);
    ee2 = exp(lsigma);
    ee3 = 1/ee1;
    ee4 = 1 + y;
    ee6 = ee4 * ee1/ee2;
    ee8 = y * ee1/ee2;
    ee9 = ee6 + 1;
    ee10 = 1 + ee8;
    ee11 = exp(lkappa2);
    ee12 = exp(logitpi);
    ee13 = R_pow(ee9,ee3);
    ee14 = R_pow(ee10,ee3);
    ee15 = 1 - 1/ee13;
    ee16 = 1 - 1/ee14;
    ee17 = exp(lkappa1);
    ee18 = 1 + ee12;
    ee19 = R_pow(ee15,ee11);
    ee20 = R_pow(ee16,ee11);
    ee21 = 1 + ee3;
    ee22 = exp(logitp);
    ee23 = ee12/ee18;
    ee24 = 1 - ee23;
    ee25 = R_pow(ee15,ee17);
    ee26 = R_pow(ee16,ee17);
    ee27 = ee11 - 1;
    ee28 = ee25 - ee26;
    ee29 = ee17 - 1;
    ee30 = ee19 - ee20;
    ee31 = R_pow(ee9, ee21);
    ee32 = R_pow(ee10,ee21);
    ee33 = ee28 * ee24;
    ee34 = ee33 * ee22;
    ee37 = ee30 * ee12/ee18 + ee20;
    ee39 = ee34 + ee19 - ee37;
    ee40 = ee3 - 1;
    ee41 = log1p(ee6);
    ee42 = log1p(ee8);
    ee43 = 2/ee1;
    ee44 = R_pow(ee15,ee27);
    ee45 = R_pow(ee16,ee27);
    ee46 = log(ee15);
    ee47 = log(ee16);
    ee48 = R_pow(ee15,ee29);
    ee49 = R_pow(ee16,ee29);
    ee50 = ee31 * ee2;
    ee51 = ee32 * ee2;
    ee52 = R_pow(ee9,ee40);
    ee53 = R_pow(ee10,ee40);
    ee55 = ee52 * ee4/ee2;
    ee57 = ee13 * ee41/ee1;
    ee59 = ee14 * ee42/ee1;
    ee61 = y * ee53/ee2;
    ee62 = ee55 - ee57;
    ee63 = R_pow(ee9,ee43);
    ee64 = R_pow(ee10,ee43);
    ee65 = ee61 - ee59;
    ee66 = ee13 * ee1;
    ee67 = ee14 * ee1;
    ee68 = ee44 * ee4;
    ee70 = 1 + ee22;
    ee72 = ee41/ee66 - ee4/ee50;
    ee74 = ee42/ee67 - y/ee51;
    ee75 = y * ee45;
    ee76 = ee68/ee31;
    ee77 = ee19 * ee46;
    ee78 = ee20 * ee47;
    ee79 = ee75/ee32;
    ee80 = ee48 * ee4;
    ee81 = y * ee49;
    ee82 = ee39 * ee18;
    ee83 = ee39 * ee2;
    ee84 = ee80/ee31;
    ee85 = ee81/ee32;
    ee86 = ee85 - ee84;
    ee87 = ee77 - ee78;
    ee88 = ee79 - ee76;
    ee90 = ee62 * ee44/ee63;
    ee92 = ee45 * ee65/ee64;
    ee93 = ee25 * ee46;
    ee94 = ee26 * ee47;
    ee97 = ee28 * ee22 + ee19 - ee20;
    ee98 = ee93 - ee94;
    ee99 = 2 * ee21;
    ee100 = ee17 * ee22;
    ee101 = 1 - ee22/ee70;
    ee102 = ee23 - 1;
    ee103 = ee48 * ee72;
    ee104 = ee49 * ee74;
    ee107 = ee62 * ee48/ee63 - ee49 * ee65/ee64;
    ee108 = ee90 - ee92;
    ee109 = ee28 * ee101;
    ee110 = (ee104 - ee103) * ee17;
    ee111 = ee44 * ee72;
    ee112 = ee45 * ee74;
    ee113 = ee17 - 2;
    ee115 = ee100 * ee86 + ee11 * ee88;
    ee116 = ee11 - 2;
    ee117 = ee97 * ee24;
    ee119 = ee109 * ee24 + ee30 * ee102/ee70;
    ee121 = ee110 * ee22 + (ee112 - ee111) * ee11;
    ee123 = ee24 * ee17 * ee22;
    ee124 = R_pow(ee82,2);
    ee125 = R_pow(ee83,2);
    ee129 = ee107 * ee24 * ee17 * ee22 + (ee90 - (ee108 * ee12/ee18 + 
                                                     ee92)) * ee11;
    ee130 = (ee76 + ee12 * ee88/ee18 - ee79) * ee11;
    ee136 = ee13 * ee21 * ee4 * ee1;
    ee137 = ee77 - (ee87 * ee12/ee18 + ee78);
    ee138 = ee123 * ee86;
    ee141 = y * ee21 * ee14 * ee1;
    ee142 = ee138 - ee130;
    ee143 = R_pow(ee50,2);
    ee144 = R_pow(ee66,2);
    ee145 = R_pow(ee51,2);
    ee146 = R_pow(ee67,2);
    ee150 = ee136/ee2 - ee31 * ee41/ee1;
    ee151 = R_pow(ee9,ee99);
    ee152 = R_pow(ee15,ee113);
    ee153 = R_pow(ee15,ee116);
    ee154 = R_pow(ee16,ee113);
    ee155 = R_pow(ee16,ee116);
    ee158 = R_pow(ee10,ee99);
   // ee159 = 1 + 3/ee1;
    ee160 = ee3 - ee99;
    ee162 = ee141/ee2 - ee32 * ee42/ee1;
    ee163 = R_pow(ee39,2);
  //  ee165 = ee117/ee39 - 1;
    ee166 = ee98 * ee24;
    ee167 = ee154 * ee29;
    ee168 = ee155 * ee27;
    ee169 = ee108 * ee11;
    ee171 = (ee150 * ee2/ee143 + 1/ee50) * ee4 - (ee55 + ee13 - 
                                                     ee57) * ee1 * ee41/ee144;
    ee175 = (ee50 - ee136)/ee143 + (ee52 * ee1 * ee41/ee144 - 
                                       1/ee31)/ee2;
    ee176 = ee62 * ee152;
    ee177 = ee62 * ee153;
    ee179 = (ee51 - ee141)/ee145 + (ee53 * ee1 * ee42/ee146 - 
                                       1/ee32)/ee2;
    //ee183 = R_pow(ee9,ee159);
    ee184 = R_pow(ee9,ee160);
    ee186 = ee48 * ee17 * ee46;
    ee187 = ee152 * ee29;
    ee189 = ee44 * ee11 * ee46;
    ee190 = ee153 * ee27;
    ee192 = ee49 * ee17 * ee47;
    ee194 = ee45 * ee11 * ee47;
    ee195 = R_pow(ee24,2);
   // ee196 = R_pow(ee10,ee159);
    ee197 = R_pow(ee10,ee160);
    ee198 = R_pow(ee4,2);
    ee200 = 1/ee82 - ee117 * ee18/ee124;
    ee203 = 1/ee70;
    ee204 = R_pow(ee46,2);
    ee205 = R_pow(ee47,2);
    ee207 = y * (1/ee51 + ee2 * ee162/ee145) - (ee14 + ee61 - 
                                                   ee59) * ee1 * ee42/ee146;
    ee208 = R_pow(y,2);
    
    //# First derivatives
    out(j, 0)= -(ee24 * ee115/ee83);// #w.r.t lsigma
    out(j, 1)=  -(ee121 * 
                    ee24/ee39); //  #w.r.t lxi
    out(j, 2)=  -(ee166 * ee17 * ee22/ee39); // #w.r.t lkappa1
    out(j, 3)=   -(ee87 * 
                     ee24 * ee11/ee39); //  #weret lkappa2
    out(j, 4)=  -(ee119 * ee22/ee39); //  #w.r.t logitp
    out(j, 5)= ee117 * ee12/ee82; // ##w.r.t logitpi
    
    //##Second derivatives 
    
    out(j, 6)=  -(((ee100 * (ee208 * 
                               (ee49 * ee21 * ee197 * ee1 - ee167/ee158) - (ee184 * 
                                                                              ee48 * ee21 * ee1 - ee187/ee151) * ee198) + ee11 * (ee208 * 
                                                                                                                                    (ee45 * ee21 * ee197 * ee1 - ee168/ee158) - (ee184 * 
                                                                                                                                                                                   ee44 * ee21 * ee1 - ee190/ee151) * ee198))/(ee39 * R_pow(ee2,2)) - 
                     (ee83 + ee138 - ee130) * ee115/ee125) * ee24); // #weret lsigma, lsigma
    out(j, 7)=  -(ee24 * 
                    (ee100 * (y * (ee179 * ee49 - ee167 * ee74/ee51) - (ee175 * 
                                                                          ee48 - ee187 * ee72/ee50) * ee4) + ee11 * (y * (ee179 * 
                                                                                                                            ee45 - ee168 * ee74/ee51) - (ee175 * ee44 - ee190 * 
                                                                                                                                                           ee72/ee50) * ee4) - ee121 * ee142/ee83)/ee39); // #weret lsigma, lxi
    out(j, 8)= -(ee123 * 
                   (y * (ee192/ee32 + ee49/ee32) - ((ee186/ee31 + ee48/ee31) * 
                                                      ee4 + ee98 * ee142/ee39))/ee83); // #weret lsigma, lkappa1
    out(j, 9)=  -(ee24 * 
                     ee11 * (y * (ee194/ee32 + ee45/ee32) - ((ee189/ee31 + 
                                                                ee44/ee31) * ee4 + ee87 * ee142/ee39))/ee83); // #weret lsigma, lkappa2
    out(j, 10)=   -((ee101 * 
                       ee24 * ee17 * ee86 + ee11 * ee102 * ee88/ee70 - ee119 * 
                       ee142/ee39) * ee22/ee83); // #weret lsigma, logitp
    out(j, 11)=  (ee115/ee82 - ee97 * 
                    ee142 * ee18/ee124) * ee24 * ee12/ee2; // #weret logitpi
    
    out(j, 12)= -((((ee154 * ee74 * 
                       ee65/ee64 - ee176 * ee72/ee63) * ee29 + ee49 * ee207 - 
                      ee171 * ee48) * ee17 * ee22 + ((ee155 * ee74 * ee65/ee64 - 
                                                        ee177 * ee72/ee63) * ee27 + ee45 * ee207 - ee171 * ee44) * 
                     ee11 - ee129 * ee121/ee39) * ee24/ee39); // #weret lxi, lxi
    out(j, 13)=  -((ee62 * 
                      (ee186/ee63 + ee48/ee63) - (ee129 * ee98/ee39 + (ee192/ee64 + 
                                                                         ee49/ee64) * ee65)) * ee24 * ee17 * ee22/ee39); // #weret lxi, lkappa1
    out(j, 14)= -((ee62 * 
                     (ee189/ee63 + ee44/ee63) - (ee129 * ee87/ee39 + (ee194/ee64 + 
                                                                        ee45/ee64) * ee65)) * ee24 * ee11/ee39); //#weret lxi, lkappa2
    out(j, 15)=  -((ee107 * 
                      ee101 * ee24 * ee17 + ee169 * ee102/ee70 - ee129 * ee119/ee39) * 
                     ee22/ee39); // #weret lxi, logitp
    out(j, 16)=   ((ee107 * ee17 * ee22 + ee169)/ee82 - 
                     ee129 * ee97 * ee18/ee124) * ee24 * ee12;  //#weret lxi, logitpi
    out(j, 17)=  -(((ee25 * ee204 - (R_pow(ee98,2)* ee24 * ee22/ee39 + 
                                       ee26 * ee205)) * ee17 + ee93 - ee94) * ee24 * ee17 * 
                     ee22/ee39); // #weret lkappa1, lkappa1
    
    out(j, 18)=  ee98 * ee87 * ee195 * ee17 * ee11 * ee22/ee163; // #weret lkappa1, lkappa2
    out(j, 19)=   -(ee98 * (1 - (ee119/ee39 + 
                                   ee203) * ee22) * ee24 * ee17 * ee22/ee39);// #weret lkappa1, logitp
    out(j, 20)= ee166 *ee200 * ee17 * ee22 * ee12; // #weret lkappa1, logitpi
    out(j, 21)= -(((ee19 * 
                      ee204 - (ee137 * ee87/ee39 + ee20 * ee205)) * ee11 + 
                     ee77 - ee78) * ee24 * ee11/ee39); // #weret lkappa2, lkappa2
    out(j, 22)=  -((ee87 * 
                      ee102/ee70 - ee119 * ee137/ee39) * ee11 * ee22/ee39); // #weret lkappa2, logitp 
    out(j, 23)= (ee87/ee82 - ee97 * ee137 * ee18/ee124) * ee24 *ee11 * ee12; // #weret lkappa2, logitpi
    out(j, 24)= -(ee119 * 
                    (1 - (ee33/ee39 + ee203) * ee22) * ee22/ee39); // #weret logitp, logitp
    out(j, 25)= -((ee119 * 
                     ee97/ee39 + ee30/ee70 - ee109) * ee24 * ee22 * ee12/ee82); // #weret logitp, logitpi
    out(j, 26)=  ee97 * (ee24/ee82 - ((ee20 - ee19) * ee24 + 
                                        ee19 - ee37) * ee12/ee124) * ee24 * ee12; // #wer.t logitpi, logitpi

}
else {  
    
    eee1 = exp(lxi);
    eee2 = exp(lsigma);
    eee3 = eee1/eee2;
    eee4 = 1 + eee3;
    eee5 = 1/eee1;
    eee6 = R_pow(eee4,eee5);
    eee7 = 1 - 1/eee6;
    eee8 = exp(logitp);
    eee9 = exp(lkappa1);
    eee10 = exp(lkappa2);
    eee11 = 1 + eee8;
    eee12 = exp(logitpi);
    eee13 = R_pow(eee7,eee9);
    eee14 = R_pow(eee7,eee10);
    eee15 = (eee13 * eee8 + eee14)/eee11;
    eee16 = eee15 + eee12;
    eee17 = 1 + eee5;
    eee18 = eee16 * eee11;
    eee19 = R_pow(eee4,eee17);
    eee20 = eee9 - 1;
    eee21 = eee10 - 1;
    eee22 = R_pow(eee7,eee20);
    eee23 = R_pow(eee7,eee21);
    eee24 = eee22 * eee9;
    eee25 = eee23 * eee10;
    eee26 = log(eee7);
    eee27 = eee24 * eee8;
    eee28 = R_pow(eee18,2);
    eee29 = log1p(eee3);
    eee30 = R_pow(eee4,(2/eee1));
    eee32 = eee27/eee19 + eee25/eee19;
    eee33 = eee19 * eee2;
    eee34 = 1 + eee12;
    eee35 = 1 - eee8/eee11;
    eee36 = R_pow(eee4,(eee5 - 1));
    eee37 = eee6 * eee1;
    eee38 = eee18 * eee2;
    eee39 = eee13 * eee35;
    eee40 = eee36/eee2;
    eee42 = eee6 * eee29/eee1;
    eee43 = 1/eee33;
    eee44 = eee27 + eee25;
    eee45 = eee40 - eee42;
    eee47 = eee29/eee37 - eee43;
    eee48 = eee14/eee11;
    eee50 = eee39 - eee48;
    eee51 = 1 - eee16/eee34;
    eee52 = R_pow(eee38,2);
    eee54 = eee27/eee30 + eee25/eee30;
    eee55 = 2 * eee17;
    eee56 = eee13 + eee12;
    eee57 = eee24 * eee26;
    eee58 = R_pow(eee7,(eee9 - 2));
    eee59 = eee25 * eee26;
    eee60 = R_pow(eee7,(eee10 - 2));
    eee62 = eee17 * eee6 * eee1;
    eee63 = R_pow(eee4,eee55);
    eee64 = 1/eee34;
   // eee66 = eee18 * eee19 * eee2;
    eee68 = eee51/eee16 + eee64;
    eee70 = eee62/eee2 - eee19 * eee29/eee1;
    eee71 = eee32 * eee13;
    eee72 = eee32 * eee14;
    eee73 = R_pow(eee33,2);
    eee74 = R_pow(eee37,2);
    eee76 = eee22 * eee35 * eee9;
  //  eee77 = eee22 + eee57;
    eee80 = eee58 * eee9 * eee20 * eee8;
    eee86 = R_pow(eee7,(eee9 + eee10)) * eee9 * eee10 * eee8 * R_pow(eee26,2)/eee28;
  //  eee87 = eee23 + eee59;
    eee89 = eee60 * eee10 * eee21;
    eee90 = eee13 * eee9;
    eee92 = eee14 * eee10 * eee26;
  //  eee93 = R_pow(eee4,(1 + 3/eee1));
    eee94 = R_pow(eee4,(eee5 - eee55));
    
    //#First derivatives 
    out(j, 0)=eee32/eee38; // # weereet lsigma
    out(j, 1)= eee44 * eee47/eee18; // # weereet lxi
    out(j, 2)= -(eee90 * eee8 * eee26/eee18); // # weereet lkappa1
    out(j, 3)= -(eee92/eee18); // # weereet lkappa2
    out(j, 4)= -(eee50 * eee8/eee18); // # weereet logitp
    out(j, 5)= -(eee51 * eee12/eee16); //  # weereet logitpi
    ////#second derivatives 
    out(j, 6)=((eee22 * eee17 * eee94 * eee1 - 
                 eee58 * eee20/eee63) * eee9 * eee8 + (eee23 * eee17 * eee94 * 
                                                    eee1 - eee60 * eee21/eee63) * eee10)/(eee18 * R_pow(eee2,2)) - (eee38 - eee32) * eee32/eee52;//  # weereet lsigma, lsigma
    out(j, 7)= (((eee33 - eee62)/eee73 + 
                  (eee36 * eee1 * eee29/eee74 - 1/eee19)/eee2) * eee44 - (eee80/eee19 + 
                                                                     eee89/eee19) * eee47/eee2)/eee18 + eee44 * eee32 * eee47/(eee28 * eee2); //  # weereet lsigma, lxi
    out(j, 8)= -((eee71 * eee26/eee28 - (eee57/eee19 + 
                                       eee22/eee19)/eee18) * eee9 * eee8/eee2); //  # weereet lsigma, lkappa1
    out(j, 9)= -((eee72 * 
                    eee26/eee28 - (eee59/eee19 + eee23/eee19)/eee18) * eee10/eee2); //# weereet lsigma, lkappa2
    out(j, 10)= -((eee32 * eee50/eee28 + (eee25/(eee11 * eee19) - 
                                        eee76/eee19)/eee18) * eee8/eee2); // # weereet lsigma logitp
    out(j, 11)= -(eee68 * 
                   eee32 * eee12/eee38); // # weereet lsigma, logitpi
    out(j, 12)= ((eee70 * eee2/eee73 + 
                   eee43 - (eee40 + eee6 - eee42) * eee1 * eee29/eee74) * eee44 + 
                  (eee80/eee30 + eee89/eee30) * eee45 * eee47)/eee18 - eee44 * eee54 * eee45 * eee47/eee28;  //# weereet lxiee lxi
    out(j, 13)= -(((eee57/eee30 + 
                     eee22/eee30)/eee18 - eee54 * eee13 * eee26/eee28) * eee45 * 
                   eee9 * eee8); // # weereet lxi, lkappa1
    out(j, 14)= -(((eee59/eee30 + eee23/eee30)/eee18 - 
                    eee54 * eee14 * eee26/eee28) * eee45 * eee10); //# weereet lxi, lkappa2
    out(j, 15)= -(((eee76/eee30 - 
                     eee25/(eee11 * eee30))/eee18 - eee54 * eee50/eee28) * eee45 * eee8); // # weereet lxi, logitp
    out(j, 16)= eee68 * eee54 * eee45 * eee12/eee18; // # weereet lxi, logitpi
    out(j, 17)= -(((eee13 + eee90 * 
                     eee26)/eee18 - R_pow(eee7,(2 * eee9)) * eee9 * eee8 * eee26/eee28) * 
                   eee9 * eee8 * eee26); // # weereet lkapap1, lkapap1
    out(j, 18)= eee86; //   # weereet lkapap1, lkapap2
    out(j, 19)= -((eee39/eee18 - 
                    eee50 * eee13 * eee8/eee28) * eee9 * eee8 * eee26); //# weereet lkapap1, logitp
    out(j, 20)= (eee51 * eee13/eee16 + eee13/eee34) * eee9 *eee8 * eee12 * eee26/eee18; //  # weereet lkapap1, logitpi
    out(j, 21) =  -(((eee14 + eee92)/eee18 - R_pow(eee7,(2 * eee10)) * 
                      eee10 * eee26/eee28) * eee10 * eee26); //# weereet lkapap2, lkapap2
    out(j, 22)= (eee50 * eee14/eee28 + eee14/(eee16 * R_pow(eee11,2))) * eee10 *eee8 * eee26; //# weereet lkapap2, logitp
    out(j, 23)= (eee51 * eee14/eee16 + eee14/eee34) *eee10 * eee12 * eee26/eee18; // # weereet lkapap2, logitpi
    out(j, 24)= -(((((eee48 - eee39) * eee8 - eee14)/eee11 + 
                     eee39)/eee18 - eee50 * eee56 * eee8/eee28) * eee8); // # weereet logitp, logitp
    out(j, 25)= eee50 * eee11 * eee8 * eee12/eee28; // # weereet logitp, logitpi
    out(j, 26)= -(eee51 * 
                   (1 - (1/eee16 + eee64) * eee12) * eee12/eee16); // # weereet logitpi, logitpi
	}
}
     
   return out;
}
    

// //' Discrete Extended generalized Pareto distribution of type 3 (zideGPD3) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each zideGPD parameter
// //' @param X1 a design matrix for the zideGPD log scale parameter
// //' @param X2 a design matrix for the zideGPD logshape parameter
// //' @param X3 a design matrix for the zideGPD log delta
// //' @param X3 a design matrix for the zideGPD logit pi
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return zidegpd3d0 a scalar, the negative log-liklihood
// //' @return zidegpd3d12 a matrix, first then second derivatives w.r.t. zideGPD3 parameters
// //' @return zidegpd3d34 a matrix, third then fourth derivatives w.r.t. zideGPD3 parameters (Not given)
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double zidegpd3d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
  
  
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec ldeltavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec logitpivec = X4 * Rcpp::as<arma::vec>(pars[3]);
  int nobs = yvec.size();
  
    if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
   ldeltavec = ldeltavec.elem(dupid);
   logitpivec = logitpivec.elem(dupid);
  }
  
  
  double y, lsigma, lxi, ldelta, logitpi;
  double e1,e2,e3, e4, e5,e6,e7,e8,e9,e10,e11,e12,e13;
  double ee1,ee2,ee3, ee4, ee5, ee6,ee7, ee8, ee9;
  double hi; 
  double lo;
  double nllh=0.0;
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    ldelta = ldeltavec[j];
    logitpi = logitpivec[j];
    if(y>0){ 
	e1=1.0/exp(lxi);
    e2 =  exp(lxi) / exp(lsigma);
    e3= 1+ (y+1)*e2;
    e4= R_pow(1/e3, e1);
    e5= 1-e4; //H(y+1)
    e6= exp(ldelta)+1; 
    e7= R_pow(e4, e6)/ exp(ldelta); //[bar H(y+1)]^(ldelta+1)//ldelta
    e8=e4/ exp(ldelta);  //[bar H(y+1)]//ldelta
	hi= e5+e7-e8;
	e9= 1+ y*e2;
	e10= R_pow(1/e9, e1);
    e11= 1-e10; //H(y)
    e12= R_pow(e10, e6)/ exp(ldelta); //[bar H(y)]^(ldelta+1)//ldelta
	e13=e10/ exp(ldelta);  //[bar H(y+1)]//ldelta
	lo= e11+e12-e13;
	nllh += -log((1/(1+exp(logitpi)))*(hi - lo)); 
	 
}else{
	 ee1=1.0/exp(lxi);
    ee2 =  exp(lxi) / exp(lsigma);
    ee3= 1+ ee2;
    ee4= R_pow(1/ee3, ee1);
    ee5= 1-ee4; //H(y+1)
    ee6= exp(ldelta)+1;
    ee7= R_pow(ee4, ee6)/ exp(ldelta); //[bar H(y+1)]^(ldelta+1)//ldelta
    ee8=ee4/ exp(ldelta);  //[bar H(y+1)]//ldelta
	ee9= ee5+ee7-ee8;
	nllh += -log(exp(logitpi)/(1+exp(logitpi))+ (1/(1+exp(logitpi)))*ee9);
}
   //if (!ISNA(nllh)){
    //nllh = 1e20;
    // break;
//}
 //Rprintf("hi %f ",hi); 
   // Rprintf("lo %f \n",lo);
  }
  return(nllh);
}


// //' @rdname zidegpd3d0
// [[Rcpp::export]]
arma::mat zidegpd3d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3,  arma::mat X4,arma::vec yvec , const arma::uvec dupid, int dcate)
{
  
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec ldeltavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec logitpivec = X4 * Rcpp::as<arma::vec>(pars[3]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs,  14);
  
   if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
    ldeltavec = ldeltavec.elem(dupid);
    logitpivec = logitpivec.elem(dupid);
 }
  
  double y, lsigma, lxi, ldelta, logitpi;
  double ee1, ee2, ee3,  ee5,ee7,ee8,  ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21,ee22, ee23, ee24, ee26,  ee27,   ee28, ee29,ee30;
  double ee31, ee32, ee33, ee34, ee35, ee36, ee136, ee37, ee38, ee39, ee40;
  double ee41,  ee45, ee46, ee47, ee48;
  double ee49, ee50,ee51,  ee56, ee57;
  double ee58, ee59, ee60, ee61;
  double  ee65, ee66;
  double ee70, ee71,ee72,  ee75, ee76;
  double ee77, ee78,ee80, ee81,ee82, ee83;
  double ee84, ee85, ee86, ee87, ee89, ee90, ee91;
  double ee92,ee93, ee94, ee101;
  double ee102,  ee106, ee111, ee112;
  double ee116, ee117, ee118, ee121,ee122,   ee126;
  double ee127,ee130, ee132, ee135,  ee140, ee141;
  double ee142, ee143, ee144, ee146, ee148, ee150;
  double ee151, ee152, ee153, ee158, ee160;
  double  ee165, ee170, ee173, ee175;
  double ee180, ee182, ee183, ee184, ee185;
  double ee186, ee190, ee192, ee194, ee195;
  
  
  double eee1, eee2, eee3,eee4, eee5, eee6, eee7, eee8;
  double eee9, eee10;
  double eee11, eee12, eee13, eee14, eee15, eee16;
  double eee17,  eee20,eee21,  eee22,eee23, eee24;
  double eee25,eee26,eee27,eee28,eee29, eee30, eee31,eee32;
  double eee34, eee35,eee36,eee37,eee38, eee40,  eee41;
  double eee42, eee43, eee44, eee45, eee46, eee47,  eee49, eee50,  eee52, eee53, eee55; 
  double eee56, eee58,  eee59, eee60, eee61,  eee63, eee64, eee65,  eee66, eee68, eee69;
  double eee71, eee72,  eee73, eee76, eee78,  eee81, eee82,  eee85, eee87;
  double eee91, eee94,  eee95;   
  
  
  
  for (int j=0; j < nobs; j++) {
    
    y  = yvec[j];
      lsigma  = lsigmavec[j];
      lxi  = lxivec[j];
      ldelta  = ldeltavec[j];
      logitpi  = logitpivec[j];
   
    if(y>0){   
     ee1 = exp(lxi);
    ee2 = exp(lsigma);
    ee3 = 1 + y;
    ee5 = ee3 * ee1/ee2;
    ee7 = y * ee1/ee2;
    ee8 = exp(ldelta);
    ee9 = 1/ee1;
    ee10 = ee5 + 1;
    ee11 = 1 + ee7;
    ee12 = 1 + ee8;
    ee13 = ee12/ee1;
    ee14 = R_pow(ee10,ee9);
    ee15 = R_pow(ee11,ee9);
    ee16 = 1/ee14;
    ee17 = 1/ee15;
    ee18 = 1 + ee9;
    ee19 = exp(logitpi);
    ee20 = R_pow(ee10,ee13);
    ee21 = R_pow(ee11,ee13);
    ee22 = log1p(ee5);
    ee23 = log1p(ee7);
    ee24 = ee13 + 1;
    ee26 = 1/ee20 + ee17;
    ee27 = ee16 + 1/ee21;
    ee28 = 1 + ee19;
    ee29 = (ee26 - ee27)/ee8;
    ee30 = ee29 + ee17;
    ee31 = R_pow(ee10,ee18);
    ee32 = R_pow(ee11,ee18);
    ee33 = ee3/ee31;
    ee34 = y/ee32;
    ee35 = ee9 - 1;
    ee36 = 2/ee1;
    ee37 = R_pow(ee10,ee24);
    ee38 = R_pow(ee11,ee24);
    ee39 = ee30 - ee16;
    ee40 = ee22/ee20;
    ee41 = ee23/ee21;
    ee45 = ee30 - (ee39 * ee19/ee28 + ee16);
    ee46 = ee13 - 1;
    ee47 = ee3/ee37;
    ee48 = 2 * ee13;
    ee49 = y/ee38;
    ee50 = (R_pow(ee10,ee35) * ee3/ee2 - ee14 * ee22/ee1)/R_pow(ee10,ee36);
    ee51 = (ee33 - ee34)/ee2;
    ee56 = (ee23/ee15 - ee22/ee14)/ee1;
    ee57 = (y * R_pow(ee11,ee35)/ee2 - ee15 * ee23/ee1)/R_pow(ee11,ee36);
    ee58 = (ee41 - ee40)/ee1;
    ee59 = 2 * ee18;
    ee60 = ee49 - ee47;
    ee61 = 2 * ee24;
    ee65 = R_pow(ee10,ee46) * ee3/ee2 - ee20 * ee22/ee1;
    ee66 = R_pow(ee10,ee48);
    ee70 = R_pow(ee11,ee48);
    ee71 = ee60/ee2;
    ee72 = 1 - ee19/ee28;
    ee75 = y * R_pow(ee11,ee46)/ee2 - ee21 * ee23/ee1;
    ee76 = (ee51 + ((ee40 - ee41)/ee1 + ee71) * ee12 + ee56)/ee8;
    ee77 = (ee27 - ee26)/ee8;
    ee78 = ee65/ee66;
    ee80 = ee76 + ee51 + ee56;
    ee81 = ee45 * ee2;
    ee82 = ee47 - ee49;
    ee83 = ee77 + ee58;
    ee84 = ee75/ee70;
    ee85 = 1/ee31;
    ee86 = 1/ee32;
    ee87 = ee45 * ee28;
    ee89 = (ee50 + (ee84 - ee78) * ee12 - ee57)/ee8 + ee50;
    ee90 = (ee82 * ee12 + ee34 - ee33)/ee8;
    ee91 = ee90 + ee34;
    ee92 = ee9 - ee59;
    ee93 = (ee12/ee37 - ee85)/ee8;
    ee94 = (ee12/ee38 - ee86)/ee8;
    ee101 = ee14 * ee18 * ee3 * ee1/ee2 - ee31 * ee22/ee1;
    ee102 = R_pow(ee10,ee59);
    //ee103 = ee24 - ee61;
    ee106 = R_pow(ee11,ee59);
    ee111 = y * ee18 * ee15 * ee1/ee2 - ee32 * ee23/ee1;
    ee112 = ee89 - ee57;
    ee116 = (ee93 - ee85) * ee3 - y * (ee94 - ee86);
    ee117 = ee91 - ee33;
    ee118 = ee80 - ee80 * ee19/ee28;
    ee121 = ee83 - ee83 * ee19/ee28;
    ee122 = ee58 - ee29;
    ee126 = ((ee50 + ee16) * ee22 + ee1 * (ee34 - ee33)/ee2 - 
                (ee57 + ee17) * ee23)/ee1;
    ee127 = ee89 - (ee112 * ee19/ee28 + ee57);
    ee130 = ee117 * ee19/ee28 + ee33;
    ee132 = ee65 * ee22/ee66;
    ee135 = R_pow(ee81,2);
    ee136 = R_pow(ee10,ee92);
    ee140 = ee13 - ee61;
    ee141 = R_pow(ee11,ee92);
    ee142 = R_pow(ee3,2);
    ee143 = ee58 - (ee122 * ee19/ee28 + ee29);
    ee144 = (y * ee111/ee106 - ee101 * ee3/ee102)/ee2;
    ee146 = R_pow(ee22,2)/ee20;
    ee148 = ee23 * ee75/ee70;
    ee150 = R_pow(ee23,2)/ee21;
    ee151 = R_pow(y,2);
    ee152 = ee91 - ee130;
    ee153 = R_pow(ee87,2);
    ee158 = ee24 * ee20 * ee3 * ee1/ee2 - ee37 * ee12 * ee22/ee1;
    //ee159 = R_pow(ee10,ee103);
    ee160 = R_pow(ee10,ee61);
    //ee164 = R_pow(ee11,ee103);
    ee165 = R_pow(ee11,ee61);
    ee170 = y * ee24 * ee21 * ee1/ee2 - ee12 * ee38 * ee23/ee1;
    ee173 = (((ee158 * ee3/ee160 - y * ee170/ee165)/ee2 + (ee82 * 
                                                              ee1/ee2 + ee12 * (ee148 - ee132) + ee41 - ee40)/ee1) * 
                ee12 + ee126 + ee144)/ee8 + ee126 + ee144;
    ee175 = ((ee132 - ee148) * ee12 + ee1 * ee60/ee2 + ee40 - 
                ee41)/ee1 + ((ee78 - ee84) * ee12 + ee57 - ee50)/ee8;
   // ee177 = (((ee159 * ee3 * ee22 - y * ee164 * ee23)/ee2 + 
   //              (ee150 - ee146)/ee1) * ee12 + ee40 - ee41)/ee1 + ee71;
    ee180 = (ee136 * ee142 - ee151 * ee141) * ee18 * ee1/ee2;
  //  ee181 = ee101/ee102;
    ee182 = (ee12 * ee22/ee37 - ee1/ee37) * ee3;
    ee183 = R_pow(ee10,ee140);
    ee184 = ee136 * ee18;
    ee185 = ee39 * ee72;
    ee186 = ee18 * ee141;
    ee190 = R_pow(ee11,ee140);
    ee192 = (ee8 * (ee146 - ee150)/ee1 + ee40 - ee41)/ee1;
   // ee193 = ee111/ee106;
    ee194 = (y * (ee23/ee32 - ee1/ee32) - ee3 * (ee22/ee31 - 
                                                    ee1/ee31))/ee1;
    ee195 = y * (ee12 * ee23/ee38 - ee1/ee38);
    //# first derivatives
    out(j, 0) = -(ee116 * ee72/ee81); //# werete lsigma
    out(j, 1) =  -(ee118/ee45); // # werete lxi
    out(j, 2) = -(ee121/ee45);// # werete ldelta
    out(j, 3) =  ee185 * ee19/ee87; // # logitpi
    //# second derivatives
    out(j, 4) =   -(((((ee24 * ee183 * ee12 - 
                          ee184)/ee8 - ee184) * ee142 - ee151 * ((ee24 * ee12 * 
                                                                    ee190 - ee186)/ee8 - ee186)) * ee1/(ee45 * R_pow(ee2,2)) - 
                       ee116 * (ee90 + ee81 + ee34 - ee130)/ee135) * ee72); //# weret (lsigma, lsigma)
    out(j, 5) =  -((((((ee182 - ee195)/ee1 + ee24 * ee1 * (ee151 * 
                                                             ee190 - ee183 * ee142)/ee2 + ee47 - ee49) * ee12 + 
                        ee180 + ee194 + ee34 - ee33)/ee8 + ee180 + ee194 + 
                       ee34 - ee33) * ee72 - ee152 * ee118/ee45)/ee81); // # weret (lsigma, xi)
    out(j, 6) = -((((ee12 * ee60 + ee33 - ee34)/ee8 + (ee195 - 
                                                         ee182)/ee1) * ee72 - ee152 * ee121/ee45)/ee81); // # weret (lsigma, ldelta)
    out(j, 7) =  (ee117/ee87 - ee152 * ee39 * ee28/ee153) * ee72 * ee19/ee2; // # weret (lsigma, logitpi)
    out(j, 8) = -((ee173 - (ee173 * 
                              ee19/ee28 + ee127 * ee118/ee45))/ee45); // # weret (lxi, lxi)
    out(j, 9) = -((ee175 - 
                      (ee175 * ee19/ee28 + ee127 * ee121/ee45))/ee45); // # weret (lxi, lkappa)
    
    out(j, 10) =  (ee112/ee87 - ee127 * ee39 * ee28/ee153) *ee72 * ee19; // # weret (lxi, logitpi)
    
    out(j, 11) = -((ee192 - (ee121 * ee143/ee45 + (ee192 - 
                                                     ee77) * ee19/ee28 + ee77))/ee45); //  # weret (lkappa, lkappa)
    
    out(j, 12) =  (ee122/ee87 - 
                     ee39 * ee143 * ee28/ee153) * ee72 * ee19; //  # weret (lkappa, logitpi)
    
    out(j, 13) = ee39 * R_pow(ee72,2) * ee19/ee87; // # weret (logitpi, logitpi)
}
else {
     	   
    eee1 = exp(lxi);
   eee2 = exp(lsigma);
   eee3 =eee1/eee2;
  eee4 = 1 +eee3;
    eee5 = exp(ldelta);
    eee6 = 1/eee1;
    eee7 = 1 + eee5;
    eee8 = eee7/eee1;
    eee9 = R_pow(eee4,eee6);
    eee10 = 1/eee9;
    eee11 = R_pow(eee4,eee8);
    eee12 = log1p(eee3);
    eee13 = 1 + eee6;
    eee14 = exp(logitpi);
    eee15 = 1/eee11;
    eee16 = R_pow(eee4,eee13);
    eee17 = (eee15 - eee10)/eee5;
    eee20 = eee17 + 1 + eee14 - eee10;
    eee21 = eee8 + 1;
    eee22 = R_pow(eee4,eee21);
    eee23 = eee11 * eee1;
    eee24 = 1/eee16;
    eee25 = eee16 * eee2;
    eee26 = eee9 * eee1;
    eee27 = eee12/eee23;
    eee28 = eee22 * eee2;
    eee29 = R_pow(eee4,(eee6 - 1));
    eee30 = 1/eee25;
    eee31 = eee7/eee22;
    eee32 = eee29/eee2;
    eee34 = eee9 * eee12/eee1;
    eee35 = eee12/eee26;
    eee36 = (eee31 - eee24)/eee5;
    eee37 = eee20 * eee2;
    eee38 = (eee32 - eee34)/R_pow(eee4,(2/eee1));
    eee40 = R_pow(eee4,(eee8 - 1));
    eee41 = 1 + eee14;
    eee42 = 1/eee28;
    eee43 = eee36 - eee24;
    eee44 = eee11 * eee12;
    eee45 = (eee40/eee2 - eee44/eee1) * eee7;
    eee46 = (eee7 * (eee27 - eee42) + eee30 - eee35)/eee5;
    eee47 = R_pow(eee23,2);
    eee49 = eee13 * eee9 * eee1;
    eee50 = (eee10 - eee15)/eee5;
    eee52 = eee46 + eee30 - eee35;
    eee53 = eee45/R_pow(eee4,(2 * eee8));
    eee55 = eee50 - eee27;
    eee56 = 1 - eee20/eee41;
    eee58 = (eee38 - eee53)/eee5 + eee38;
    eee59 = R_pow(eee37,2);
    eee60 = R_pow(eee25,2);
    eee61 = R_pow(eee26,2);
    eee63 = eee49/eee2 - eee16 * eee12/eee1;
    eee64 = eee17 + eee27;
    eee65 = 1/eee41;
    eee66 = 2 * eee13;
    eee68 = eee21 * eee11 * eee1;
    eee69 = R_pow(eee28,2);
    eee71 = eee56/eee20 + eee65;
    eee72 = 1/eee22;
    eee73 = 2 * eee21;
    eee76 = (eee45 + eee11) * eee1 * eee12/eee47;
    eee78 = eee63 * eee2/eee60;
   // eee79 = eee63/R_pow(eee4,eee66);
    eee81 = eee68/eee2 - eee7 * eee22 * eee12/eee1;
    eee82 = (eee25 - eee49)/eee60;
    eee85 = (eee32 + eee9 - eee34) * eee1 * eee12/eee61;
  //  eee86 = R_pow(eee20,2);
    eee87 = eee13 * R_pow(eee4,(eee6 - eee66));
    eee91 = eee7 * eee40 * eee1 * eee12/eee47;
    eee94 = eee29 * eee1 * eee12/eee61;
    eee95 = 1/eee23;
    //# first derivatives
    out(j, 0) =  -(eee43/eee37);//    # werete lsigma
    out(j, 1) =  -(eee52/eee20); //# werete lxi
    out(j, 2) =  -(eee55/eee20);// # werete ldelta
    out(j, 3) =   -(eee56 * eee14/eee20);// # werete logitpi
    
    //# second derivatives
    out(j, 4) = -(((eee21 * 
                      eee7 * R_pow(eee4,(eee8 - eee73)) - eee87)/eee5 - eee87) * eee1/(eee20 * 
                                                                          R_pow(eee2,2)) - eee43 * (eee36 + eee37 - eee24)/eee59);//  # weret (lsigma, lsigma)
    out(j, 5) =  -(((((eee91 - 
                         eee72)/eee2 + (eee28 - eee68)/eee69) * eee7 - (eee82 + (eee94 - 
                                                                            eee24)/eee2))/eee5 - ((eee52 * eee43/eee20 + eee94 - eee24)/eee2 + 
                                                                                                 eee82))/eee20); //# weret (lsigma, xi)
    out(j, 6) = -(((eee24 - eee31)/eee5 + eee72 - 
                    (eee43 * eee55/eee20 + eee91))/eee37);// # weret (lsigma, lkappa)
    out(j, 7) = eee43 * eee71 * eee14/eee37;// # weret (lsigma, logitpi)
    out(j, 8) =  -((((eee81 * 
                        eee2/eee69 + eee42 - eee76) * eee7 + eee85 - (eee78 + eee30))/eee5 + 
                      eee85 - (eee58 * eee52/eee20 + eee78 + eee30))/eee20); //# weret (lxi, lxi)
    out(j, 9) =  -((eee76 + 
                       (eee53 - eee38)/eee5 - (eee58 * eee55/eee20 + eee42))/eee20); //# weret (lxi, lkappa)
    
    out(j, 10) =  eee58 * eee71 * eee14/eee20; //# weret (lxi, logitpi)
    out(j, 11) = -(((eee11 * eee5 * eee12/eee47 + eee95) * eee12 + 
                     eee64 * eee55/eee20 - eee50)/eee20); //   # weret (lkappa, lkappa)
    out(j, 12) =   -(eee71 * 
                       eee64 * eee14/eee20); //# weret ( lkappa, logitpi)
    out(j, 13) =   -(eee56 * (1 - (1/eee20 + eee65) * eee14) * eee14/eee20); //# weret (logitpi, logitpi)
    }
}
     
   return out;
}
    


// //' Discrete Extended generalized Pareto distribution of type 4 (zideGPD4) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each zideGPD parameter
// //' @param X1 a design matrix for the zideGPD log scale parameter
// //' @param X2 a design matrix for the zideGPD log shape parameter
// //' @param X3 a design matrix for the zideGPD log kappa
// //' @param X4 a design matrix for the zideGPD log delta
// //' @param X5 a design matrix for the zideGPD logit pi
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return zidegpd1d0 a scalar, the negative log-liklihood
// //' @return zidegpd4d12 a matrix, first then second derivatives w.r.t. zideGPD4 parameters
// //' @return zidegpd4d34 a matrix, third then fourth derivatives w.r.t. zideGPD1 parameters (Not given)
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double zidegpd4d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4 , const arma::mat& X5, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec ldeltavec = X4 * Rcpp::as<arma::vec>(pars[3]);
  arma::vec logitpivec = X5 * Rcpp::as<arma::vec>(pars[4]);
  int nobs = yvec.size();
   
    if (dcate == 1) {
   lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
     lkappavec = lkappavec.elem(dupid);
    ldeltavec = ldeltavec.elem(dupid);
    logitpivec = logitpivec.elem(dupid);
  }
  
  double y, lsigma, lxi,lkappa, ldelta, logitpi;
  double e1,e2,e3, e4, e5,e6,e7,e8,e9,e10,e11,e12,e13;
  double ee1,ee2,ee3, ee4, ee5, ee6,ee7, ee8, ee9;
  double hi; 
  double lo;
  double nllh=0.0;
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    lsigma = lsigmavec[j];
    lxi = lxivec[j];
    lkappa = lkappavec[j];
    ldelta = ldeltavec[j];
    logitpi = logitpivec[j];
    if(y>0){ 
	e1=1.0/exp(lxi);
    e2 =  exp(lxi) / exp(lsigma);
    e3= 1+ (y+1)*e2;
    e4= R_pow(1/e3, e1);
    e5= 1-e4; //H(y+1)
    e6= exp(ldelta)+1; 
    e7= R_pow(e4, e6)/ exp(ldelta); //[bar H(y+1)]^(ldelta+1)//ldelta
    e8=e4/ exp(ldelta);  //[bar H(y+1)]//ldelta
	hi= R_pow((e5+e7-e8),exp(lkappa)/2);
	e9= 1+ y*e2;
	e10= R_pow(1/e9, e1);
    e11= 1-e10; //H(y)
    e12= R_pow(e10, e6)/ exp(ldelta); //[bar H(y)]^(ldelta+1)//ldelta
	e13=e10/ exp(ldelta);  //[bar H(y+1)]//ldelta
	lo= R_pow((e11+e12-e13),exp(lkappa)/2);
    
	nllh += -log((1/(1+exp(logitpi)))*(hi - lo)); 
	 
}else{
	ee1=1.0/exp(lxi);
    ee2 =  exp(lxi) / exp(lsigma);
    ee3= 1+ ee2;
    ee4= R_pow(1/ee3, ee1);
    ee5= 1-ee4; //H(y+1)
    ee6= exp(ldelta)+1; 
    ee7= R_pow(ee4, ee6)/ exp(ldelta); //[bar H(y+1)]^(ldelta+1)//ldelta
    ee8=ee4/ exp(ldelta);  //[bar H(y+1)]//ldelta
	ee9= R_pow((ee5+ee7-ee8),exp(lkappa)/2);  
	nllh += -log(exp(logitpi)/(1+exp(logitpi))+ (1/(1+exp(logitpi)))*ee9);
}
   //if (!ISNA(nllh)){
    //nllh = 1e20;
    // break;
//}
 //Rprintf("hi %f ",hi); 
   // Rprintf("lo %f \n",lo);
  }
  return(nllh);
}


// //' @rdname zidegpd4d0
// [[Rcpp::export]]
arma::mat zidegpd4d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3,  arma::mat X4, arma::mat X5, arma::vec yvec , const arma::uvec dupid, int dcate)
{
  
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lxivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec ldeltavec = X4 * Rcpp::as<arma::vec>(pars[3]);
  arma::vec logitpivec = X5 * Rcpp::as<arma::vec>(pars[4]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs,  20);
  
   if (dcate == 1) {
    lsigmavec = lsigmavec.elem(dupid);
    lxivec = lxivec.elem(dupid);
    lkappavec = lkappavec.elem(dupid);
    ldeltavec = ldeltavec.elem(dupid);
    logitpivec = logitpivec.elem(dupid);
  }
  
  double y, lsigma, lxi,lkappa, ldelta, logitpi;
  double ee1, ee2, ee3, ee4, ee5,ee7,  ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21,ee22, ee23, ee24, ee25,ee26,     ee28, ee30;
  double ee31, ee32, ee33, ee34, ee35, ee36, ee136, ee37, ee38, ee39, ee40;
  double ee41,  ee42,ee43,ee44,ee45;
  double ee49, ee50,ee51,ee52,ee53,ee54,ee55,  ee56, ee57;
  double ee58, ee59, ee60, ee61, ee63;
  double  ee65, ee67,ee69;
  double ee70, ee71,ee72 ,ee73, ee75, ee76;
  double ee77, ee78,ee79,ee80, ee81,ee82, ee83;
  double ee84, ee85, ee86, ee87,ee88, ee89, ee90, ee91;
  double ee92, ee94,ee95,ee96,ee97,ee98, ee99,ee100,ee101;
  double ee104, ee105, ee106,ee107, ee108,ee109,ee110, ee112,ee114,ee115;
  double ee116, ee117, ee118, ee119,ee120,ee121,ee122,ee123,ee124,ee125, ee128,  ee126;
  double ee130, ee132, ee134,ee138,ee139,  ee140, ee141;
  double ee142, ee143, ee144,ee145, ee146,ee147, ee148, ee149,ee150;
  double ee151, ee152, ee153, ee154,ee155,ee158,ee159, ee160, ee161,ee164;
  double  ee168,   ee172,ee173, ee174,  ee178;
  double  ee182, ee183, ee184, ee185;
  double ee186, ee187,ee191, ee194, ee196,ee197,ee202, ee203,ee207;
  double ee211, ee212,ee213,ee216, ee219, ee220, ee221,ee222,ee224, ee227, ee231,ee233;
  double ee235,ee238, ee242,ee244,ee245, ee249, ee253,ee257, ee258,ee262, ee263;
  double eee1, eee2, eee3,eee4, eee5, eee6, eee7, eee8;
  double eee9, eee10;
  double eee11, eee12, eee13, eee14,  eee16;
  double eee17, eee18, eee19,eee20,eee21,  eee22,eee23, eee24;
  double eee25,eee26,eee27,eee28,eee29, eee30, eee31,eee32;
  double eee33,eee34, eee35,eee36,eee38, eee39,eee40,  eee41;
  double eee42, eee43, eee44, eee45,  eee47,  eee48,eee49, eee50,eee51,  eee52, eee53,eee54, eee55; 
  double eee57,eee58,  eee59, eee60, eee61,  eee63, eee64, eee65,  eee67, eee68, eee70;
  double eee71, eee72,  eee73,eee74, eee76, eee77, eee79, eee81, eee82,eee83, eee84,  eee87;
  double eee90,eee91,eee92, eee93, eee94,eee96, eee98, eee99,eee100,eee103 ,eee105,eee108,eee111;
  double eee114,eee116,eee119,eee121,eee122,eee123,eee127,eee128;   
  
  
  
  for (int j=0; j < nobs; j++) {
    
    y  = yvec[j];
      lsigma  = lsigmavec[j];
      lxi  = lxivec[j];
      lkappa = lkappavec[j];
      ldelta  = ldeltavec[j];
      logitpi  = logitpivec[j];
   
    if(y>0){   
    ee1 = exp(lxi);
    ee2 = exp(lsigma);
    ee3 = exp(ldelta);
    ee4 = 1/ee1;
    ee5 = 1 + y;
    ee7 = ee5 * ee1/ee2;
    ee9 = y * ee1/ee2;
    ee10 = ee7 + 1;
    ee11 = 1 + ee9;
    ee12 = 1 + ee3;
    ee13 = ee12/ee1;
    ee14 = R_pow(ee10,ee4);
    ee15 = R_pow(ee11,ee4);
    ee16 = 1/ee14;
    ee17 = 1/ee15;
    ee18 = exp(lkappa);
    ee19 = ee18/2;
    ee20 = R_pow(ee10,ee13);
    ee21 = R_pow(ee11,ee13);
    ee22 = 1/ee20;
    ee23 = 1/ee21;
    ee24 = (ee22 - ee16)/ee3;
    ee25 = (ee23 - ee17)/ee3;
    ee26 = exp(logitpi);
    ee28 = ee24 + 1 - ee16;
    ee30 = ee25 + 1 - ee17;
    ee31 = 1 + ee4;
    ee32 = ee19 - 1;
    ee33 = R_pow(ee28,ee19);
    ee34 = R_pow(ee30,ee19);
    ee35 = log1p(ee7);
    ee36 = log1p(ee9);
    ee37 = 1 + ee26;
    ee38 = ee13 + 1;
    ee39 = R_pow(ee10,ee31);
    ee40 = R_pow(ee11,ee31);
    ee41 = R_pow(ee28,ee32);
    ee42 = R_pow(ee30,ee32);
    ee43 = ee4 - 1;
    ee44 = ee33 - ee34;
    ee45 = 2/ee1;
    ee49 = ee33 - (ee44 * ee26/ee37 + ee34);
    ee50 = 1/ee39;
    ee51 = 1/ee40;
    ee52 = R_pow(ee10,ee38);
    ee53 = R_pow(ee11,ee38);
    ee54 = ee20 * ee1;
    ee55 = ee21 * ee1;
    ee56 = ee13 - 1;
    ee57 = ee35/ee54;
    ee58 = ee36/ee55;
    ee59 = R_pow(ee10,ee43);
    ee60 = R_pow(ee11,ee43);
    ee61 = 1 - ee26/ee37;
    ee63 = ee59 * ee5/ee2;
    ee65 = ee14 * ee35/ee1;
    ee67 = ee15 * ee36/ee1;
    ee69 = y * ee60/ee2;
    ee70 = 2 * ee13;
    ee71 = (ee63 - ee65)/R_pow(ee10,ee45);
    ee72 = ee39 * ee2;
    ee73 = ee40 * ee2;
    ee75 = (ee69 - ee67)/R_pow(ee11,ee45);
    ee76 = ee14 * ee1;
    ee77 = ee12/ee52;
    ee78 = ee12/ee53;
    ee79 = ee15 * ee1;
    ee80 = (ee77 - ee50)/ee3;
    ee81 = (ee78 - ee51)/ee3;
    ee82 = ee80 - ee50;
    ee83 = ee81 - ee51;
    ee84 = log(ee28);
    ee85 = log(ee30);
    ee86 = ee5/ee72;
    ee87 = ee35/ee76;
    ee88 = ee36/ee79;
    ee89 = y/ee73;
    ee90 = R_pow(ee10,ee56);
    ee91 = R_pow(ee11,ee56);
    ee92 = y * ee83;
    ee94 = ee82 * ee41 * ee5;
    ee95 = ee52 * ee2;
    ee96 = ee20 * ee35;
    ee97 = ee53 * ee2;
    ee98 = ee21 * ee36;
    ee99 = ee92 * ee42;
    ee100 = ee94/2;
    ee101 = (ee90 * ee5/ee2 - ee96/ee1) * ee12;
    ee104 = ee12 * (y * ee91/ee2 - ee98/ee1);
    ee105 = 0.5 * (ee33 * ee84);
    ee106 = 0.5 * (ee34 * ee85);
    ee107 = ee19 - 2;
    ee108 = ee99/2;
    ee109 = ee101/R_pow(ee10,ee70);
    ee110 = ee104/R_pow(ee11,ee70);
    ee112 = (ee71 - ee109)/ee3 + ee71;
    ee114 = (ee75 - ee110)/ee3 + ee75;
    ee115 = ee24 + ee57;
    ee116 = ee25 + ee58;
    ee117 = ee49 * ee2;
    ee118 = ee5/ee95;
    ee119 = y/ee97;
    ee120 = ee49 * ee37;
    ee121 = ee100 - ee108;
    ee122 = (ee12 * (ee57 - ee118) + ee86 - ee87)/ee3;
    ee123 = (ee12 * (ee58 - ee119) + ee89 - ee88)/ee3;
    ee124 = (ee16 - ee22)/ee3;
    ee125 = (ee17 - ee23)/ee3;
    ee126 = ee105 - ee106;
    ee128 = ee112 * ee41/2;
    ee130 = ee114 * ee42/2;
    ee132 = ee122 + ee86 - ee87;
    ee134 = ee123 + ee89 - ee88;
    ee136 = ee41 * ee115/2;
    ee138 = ee42 * ee116/2;
    ee139 = ee124 - ee57;
    ee140 = ee125 - ee58;
    ee141 = R_pow(ee28,ee107);
    ee142 = R_pow(ee30,ee107);
    ee143 = ee132 * ee41;
    ee144 = ee134 * ee42;
    ee145 = ee41 * ee139;
    ee146 = ee42 * ee140;
    ee147 = 2 * ee31;
    ee148 = ee143/2;
    ee149 = ee144/2;
    ee150 = ee145/2;
    ee151 = ee146/2;
    ee152 = ee128 - ee130;
    ee153 = ee148 - ee149;
    ee154 = R_pow(ee54,2);
    ee155 = R_pow(ee55,2);
    ee158 = ee14 * ee31 * ee5 * ee1;
    ee159 = ee150 - ee151;
    ee160 = ee138 - ee136;
    ee161 = 2 * ee38;
    ee164 = y * ee31 * ee15 * ee1;
    ee168 = ee128 - (ee152 * ee26/ee37 + ee130);
    ee172 = R_pow(ee117,2);
    ee173 = ee100 - (ee121 * ee26/ee37 + ee108);
    ee174 = ee44 * ee61;
    ee178 = ee160 * ee26/ee37 + ee136 - ee138;
    ee182 = ee105 - (ee126 * ee26/ee37 + ee106);
    ee183 = R_pow(ee120,2);
    ee184 = R_pow(ee72,2);
    ee185 = R_pow(ee76,2);
    ee186 = R_pow(ee73,2);
    ee187 = R_pow(ee79,2);
    ee191 = ee158/ee2 - ee39 * ee35/ee1;
    ee194 = ee4 - ee147;
    ee196 = ee164/ee2 - ee40 * ee36/ee1;
    ee197 = R_pow(ee95,2);
   // ee199 = ee174/ee49 - 1;
    ee202 = ee38 * ee20 * ee5 * ee1;
    ee203 = R_pow(ee97,2);
    ee207 = ee41 + ee41 * ee18 * ee84/2;
    ee211 = ee42 + ee42 * ee18 * ee85/2;
    ee212 = 1/ee52;
    ee213 = 1/ee53;
    ee216 = y * ee38 * ee21 * ee1;
    ee219 = (ee101 + ee20) * ee1 * ee35/ee154;
    ee220 = (ee191 * ee2/ee184 + 1/ee72) * ee5;
    ee221 = ee82 * ee141;
    ee222 = ee83 * ee142;
    ee224 = (ee72 - ee158)/ee184 + (ee59 * ee1 * ee35/ee185 - 
                                       ee50)/ee2;
    ee227 = (ee63 + ee14 - ee65) * ee1 * ee35/ee185;
  //  ee228 = ee191/R_pow(ee10,ee147);
    ee231 = (ee104 + ee21) * ee1 * ee36/ee155;
    ee233 = ee202/ee2 - ee52 * ee12 * ee35/ee1;
    ee235 = (ee73 - ee164)/ee186 + (ee60 * ee1 * ee36/ee187 - 
                                       ee51)/ee2;
    ee238 = (ee15 + ee69 - ee67) * ee1 * ee36/ee187;
    ee242 = ee90 * ee12 * ee1 * ee35/ee154;
    ee244 = R_pow(ee10,ee194) * ee31;
    ee245 = ee31 * R_pow(ee11,ee194);
    ee249 = ee12 * ee91 * ee1 * ee36/ee155;
    ee253 = ee13 - ee161;
    //ee254 = ee38 - ee161;
   // ee256 = ee196/R_pow(ee11,ee147);
    ee257 = 1/ee54;
    ee258 = 1/ee55;
    ee262 = ee216/ee2 - ee12 * ee53 * ee36/ee1;
    ee263 = y * (1/ee73 + ee2 * ee196/ee186); 
    //# first derivatives
    out(j, 0) = -(ee121 * ee61 * ee18/ee117); // # werete lsigma
    out(j, 1) = -(ee153 * 
                    ee61 * ee18/ee49); // # werete lxi
    out(j, 2) =-(ee126 * ee61 * ee18/ee49); // # werete lkappa
    out(j, 3) = -(ee159 * ee61 * ee18/ee49); // # werete ldelta
    out(j, 4) = ee174 *ee26/ee120; //  # werete logitpi

  //  # second derivatives
    out(j, 5) = -(((0.5 * 
                      ((((ee38 * R_pow(ee10,ee253) * ee12 - ee244)/ee3 - ee244) * 
                          ee41 * ee1 + R_pow(ee82,2) * ee141 * ee32) * R_pow(ee5,2)) - 0e5 * 
                      (R_pow(y,2) * (((ee38 * ee12 * R_pow(ee11,ee253) - ee245)/ee3 - ee245) * 
                                ee42 * ee1 + R_pow(ee83,2) * ee142 * ee32)))/(ee49 * R_pow(ee2,2)) - 
                     (ee173 * ee18 + ee117) * ee121/ee172) * ee61 * ee18); //# weret (lsigma, lsigma)
    
    out(j, 6) = -((0.5 * ((((((ee242 - ee212)/ee2 + (ee95 - ee202)/ee197) * 
                               ee12 - ee224)/ee3 - ee224) * ee41 + ee132 * ee82 * 
                             ee141 * ee32/ee2) * ee5) - (ee153 * ee173 * ee18/ee117 + 
                                                           0.5 * (y * (((((ee249 - ee213)/ee2 + (ee97 - ee216)/ee203) * 
                                                                           ee12 - ee235)/ee3 - ee235) * ee42 + ee134 * ee83 * 
                                                                         ee142 * ee32/ee2)))) * ee61 * ee18/ee49); // # weret (lsigma, xi)
    out(j, 7) = -((0.5 * 
                     (ee82 * ee207 * ee5) - (ee173 * ee126 * ee18/ee49 + 
                                               0.5 * (ee92 * ee211))) * ee61 * ee18/ee117); // # weret (lsigma, lkappa)
    out(j, 8) = -((0.5 * 
                     ((ee221 * ee139 * ee32 + ee41 * ((ee50 - ee77)/ee3 + 
                                                        ee212 - ee242)) * ee5) - (ee173 * ee159 * ee18/ee49 + 
                                                                                    0.5 * (y * (ee222 * ee140 * ee32 + ee42 * ((ee51 - 
                                                                                                                                  ee78)/ee3 + ee213 - ee249))))) * ee61 * ee18/ee117); //# weret (lsigma, ldelta)
    out(j, 9) = (ee121/ee120 - ee173 * ee44 * ee37/ee183) *ee61 * ee18 * ee26/ee2; //  # weret (lsigma, logitpi)
    out(j, 10) =  -((0.5 * (((((ee233 * ee2/ee197 + 
                                  1/ee95) * ee5 - ee219) * ee12 + ee227 - ee220)/ee3 + 
                               ee227 - ee220) * ee41 + ee112 * ee132 * ee141 * ee32) - 
                       (ee168 * ee153 * ee18/ee49 + 0.5 * (((ee238 + ee12 * 
                                                               (y * (1/ee97 + ee2 * ee262/ee203) - ee231) - ee263)/ee3 + 
                                                              ee238 - ee263) * ee42 + ee134 * ee114 * ee142 * ee32))) * 
                      ee61 * ee18/ee49); //# weret (lxi, lxi)
    
    out(j, 11) = -((0.5 * (ee112 * ee207) - 
                      (ee168 * ee126 * ee18/ee49 + 0.5 * (ee211 * ee114))) * 
                     ee61 * ee18/ee49); //# weret (lxi, lkappa)
    
    out(j, 12)=  -((0.5 * ((ee219 + (ee109 - 
                                       ee71)/ee3 - ee118) * ee41 + ee112 * ee141 * ee139 * ee32) - 
                      (ee168 * ee159 * ee18/ee49 + 0.5 * ((ee231 + (ee110 - 
                                                                      ee75)/ee3 - ee119) * ee42 + ee114 * ee142 * ee140 * 
                                                            ee32))) * ee61 * ee18/ee49); // # weret (lxi, ldelta)
    
    
    out(j, 13) = (ee152/ee120 - 
                    ee168 * ee44 * ee37/ee183) * ee61 * ee18 * ee26; // # weret (lxi, logitpi)
    
    out(j, 14) =  -(((0.25 * (ee33 * R_pow(ee84,2)) - 
                        (ee182 * ee126/ee49 + 0.25 * (ee34 * R_pow(ee85,2)))) * ee18 + 
                       ee105 - ee106) * ee61 * ee18/ee49); //# weret (lkappa, lkappa)
    
    out(j, 15) = -((ee150 + 
                      (0.25 * (ee145 * ee84) - (ee159 * ee182/ee49 + 0.25 * 
                                                  (ee146 * ee85))) * ee18 - ee151) * ee61 * ee18/ee49); // # weret (lkappa, ldelta)
    out(j, 16) =   (ee126/ee120 - ee44 * ee182 * ee37/ee183) *ee61 * ee18 * ee26; //   # weret (lkappa, logitpi)
    out(j, 17) = -((ee178 * ee159 * 
                      ee18/ee49 + 0.5 * (((ee20 * ee3 * ee35/ee154 + ee257) * 
                                            ee35 - ee124) * ee41 - ee141 * ee115 * ee139 * ee32) - 
                      0.5 * (((ee21 * ee3 * ee36/ee155 + ee258) * ee36 - ee125) * 
                               ee42 - ee142 * ee116 * ee140 * ee32)) * ee61 * ee18/ee49); //  # weret (ldelta, ldelta)
    out(j, 18) = (ee178 * ee44 * ee37/ee183 + ee160/ee120) *ee61 * ee18 * ee26; // # weret (ldelta, logitpi)
    out(j, 19) = ee44 *R_pow(ee61,2) * ee26/ee120; // # w.r.t (logitpi, logitpi)
}
else {  
    eee1 = exp(lxi);
    eee2 = exp(lsigma);
    eee3 = eee1/eee2;
    eee4 = 1 + eee3;
    eee5 = exp(ldelta);
    eee6 = 1/eee1;
    eee7 = R_pow(eee4,eee6);
    eee8 = 1 + eee5;
    eee9 = 1/eee7;
    eee10 = eee8/eee1;
    eee11 = R_pow(eee4,eee10);
    eee12 = exp(lkappa);
    eee13 = 1/eee11;
    eee14 = (eee13 - eee9)/eee5;
    eee16 = eee14 + 1 - eee9;
    eee17 = eee12/2;
    eee18 = log1p(eee3);
    eee19 = exp(logitpi);
    eee20 = 1 + eee6;
    eee21 = R_pow(eee4,eee20);
    eee22 = R_pow(eee16,eee17);
    eee23 = eee22 + eee19;
    eee24 = eee17 - 1;
    eee25 = eee10 + 1;
    eee26 = R_pow(eee4,eee25);
    eee27 = 1/eee21;
    eee28 = R_pow(eee16, eee24);
    eee29 = eee11 * eee1;
    eee30 = eee21 * eee2;
    eee31 = eee18/eee29;
    eee32 = 2 * eee23;
    eee33 = eee7 * eee1;
    eee34 = R_pow(eee4, (eee6 - 1));
    eee35 = 1/eee30;
    eee36 = eee34/eee2;
    eee38 = eee7 * eee18/eee1;
    eee39 = eee18/eee33;
    eee40 = eee8/eee26;
    eee41 = eee26 * eee2;
    eee42 = (eee40 - eee27)/eee5;
    eee43 = (eee36 - eee38)/R_pow(eee4,(2/eee1));
    eee44 = eee42 - eee27;
    eee45 = log(eee16);
    eee47 = R_pow(eee4,(eee10 - 1));
    eee48 = 1/eee41;
    eee49 = eee23 * eee2;
    eee50 = eee11 * eee18;
    eee51 = (eee47/eee2 - eee50/eee1) * eee8;
    eee52 = 1 + eee19;
    eee53 = (eee8 * (eee31 - eee48) + eee35 - eee39)/eee5;
    eee54 = (eee9 - eee13)/eee5;
    eee55 = R_pow(eee32,2);
    eee57 = eee53 + eee35 - eee39;
    eee58 = eee51/R_pow(eee4,(2 * eee10));
    eee59 = R_pow(eee16,(eee17 - 2));
    eee60 = eee54 - eee31;
    eee61 = 2 * eee49;
    eee63 = (eee43 - eee58)/eee5 + eee43;
    eee64 = R_pow(eee16, (2 * eee24));
    eee65 = eee14 + eee31;
    eee67 = R_pow(eee29,2);
    eee68 = R_pow(eee16,(eee12 - 1));
    eee70 = eee20 * eee7 * eee1;
    eee71 = 1 - eee23/eee52;
    eee72 = R_pow(eee61,2);
    eee73 = R_pow(eee30,2);
    eee74 = R_pow(eee33,2);
    eee76 = eee70/eee2 - eee21 * eee18/eee1;
    eee77 = 2 * eee20;
    eee79 = eee44 * eee28 * eee12;
    eee81 = eee25 * eee11 * eee1;
    eee82 = R_pow(eee41,2);
    eee83 = eee68 * eee12;
    eee84 = eee68/eee32;
    eee87 = eee28 * eee71/eee32 + eee28/(2 * eee52);
    eee90 = eee28 + 0.5 * (eee28 * eee12 * eee45);
    eee91 = eee28/2;
    eee92 = 1/eee26;
    eee93 = 2 * eee25;
    eee94 = eee63 * eee57;
 // eee95 = eee63 * eee44;
    eee96 = eee57 * eee44; 
    eee98 = eee57 * eee28 * eee12; 
    eee99 = eee44 * eee64;
    eee100 = eee44 * eee59;
    eee103 = (eee51 + eee11) * eee1 * eee18/eee67;
    eee105 = eee90/eee32 - eee83 * eee45/eee55;
    eee108 = (eee91 - eee84) * eee12 * eee45 + eee28;
    eee111 = eee76 * eee2/eee73 + eee35;
  //  eee112 = eee76/R_pow(eee4, eee77);
    eee114 = eee81/eee2 - eee8 * eee26 * eee18/eee1;
    eee116 = (eee30 - eee70)/eee73 + (eee34 * eee1 * eee18/eee74 - eee27)/eee2;
    eee119 = (eee36 + eee7 - eee38) * eee1 * eee18/eee74;
    eee121 = eee28 * eee60 * eee12;
    eee122 = eee22 * eee12;
    eee123 = eee20 * R_pow(eee4,(eee6 - eee77));
    eee127 = eee8 * eee47 * eee1 * eee18/eee67;
    eee128 = 1/eee29;
    
    //# first derivatives
    out(j, 0) = -(eee79/eee61); // # weereetee lsigma
    out(j, 1) = -(eee98/eee32);  //# weereetee lxi
    out(j, 2) = -(0.5 * 
                    (eee122 * eee45/eee23)); // # weereetee lkappa
    out(j, 3) = -(eee121/eee32); // # weereetee ldelta
    out(j, 4) = -(eee71 * 
                    eee19/eee23); // # weereetee logitpi
    
    //# second derivatives
    out(j, 5) = -(((((eee25 * 
                        eee8 * R_pow(eee4,(eee10 - eee93)) - eee123)/eee5 - eee123) * eee28 * 
                      eee1 + R_pow(eee44,2) * eee59 * eee24)/(2 * (eee23 * R_pow(eee2,2))) - 2 * 
                     ((eee79/2 + eee49) * eee44 * eee28/eee72)) * eee12); //# weereet (lsigma, lsigma)
    
    out(j, 6) = -(((((((eee127 - 
                          eee92)/eee2 + (eee41 - eee81)/eee82) * eee8 - eee116)/eee5 - 
                       eee116) * eee28 + eee96 * eee59 * eee24/eee2)/eee32 - eee96 * 
                     eee64 * eee12/(eee55 * eee2)) * eee12); // # weereet (lsigma, xi)
    out(j, 7) = -(0.5 * (eee108 * 
                           eee44 * eee12/eee49)); // # weereet (lsigma, lkappa)
    out(j, 8) = -(((eee100 * eee60 * eee24 + 
                      eee28 * ((eee27 - eee40)/eee5 + eee92 - eee127))/eee32 - eee99 * 
                     eee60 * eee12/eee55) * eee12/eee2); // # weereet (lsigma, ldelta)
    out(j, 9) = eee44 * eee87 * eee12 * eee19/eee49; //  # weereet (lsigma, logitpi)
    out(j, 10) =  -((((((eee114 * 
                           eee2/eee82 + eee48 - eee103) * eee8 + eee119 - eee111)/eee5 + 
                         eee119 - eee111) * eee28 + eee94 * eee59 * eee24)/eee32 - eee94 * 
                       eee64 * eee12/eee55) * eee12); //# weereet (lxi, lxi)
    
     out(j, 11) = -(0.5 * (eee63 * eee108 * 
                            eee12/eee23)); // # weereet (lxi, lkappa)
    
    
    
     out(j, 12)=  -((((eee103 + (eee58 - eee43)/eee5 - 
                        eee48) * eee28 + eee63 * eee59 * eee60 * eee24)/eee32 - eee63 * 
                      eee64 * eee60 * eee12/eee55) * eee12); // # weereet (lxi, ldelta)
    
    
     out(j, 13) = eee63 * eee87 * eee12 * eee19/eee23; // # weereet (lxi, logitpi)
    
     out(j, 14) =  -(0.5 * ((eee22 + (0.5 * eee22 - 
                                      0.5 * (R_pow(eee16,eee12)/eee23)) * eee12 * eee45) * eee12 * eee45/eee23)); //# weereet (lkappa, lkappa)
    
     out(j, 15) = -(eee105 * eee60 * eee12); // # weereet (lkappa, ldelta)
     out(j, 16) =  (0.5 * (eee22 * 
                            eee71/eee23) + 0.5 * (eee22/eee52)) * eee12 * eee19 * eee45/eee23; //   # weereet (lkappa, logitpi)
     out(j, 17) = -(((((eee11 * eee5 * eee18/eee67 + eee128) * 
                        eee18 - eee54) * eee28 - eee59 * eee65 * eee60 * eee24)/eee32 + 
                      eee64 * eee65 * eee60 * eee12/eee55) * eee12); //  # weereet (ldelta, ldelta)
     out(j, 18) = -(eee87 * 
                     eee65 * eee12 * eee19/eee23); // # weereet (ldelta, logitpi)
     out(j, 19) = -(eee71 * (1 - 
                             (1/eee23 + 1/eee52) * eee19) * eee19/eee23); // # weereet (logitpi, logitpi)
	}
}
     
   return out;
}
    

    

    

    

    


 

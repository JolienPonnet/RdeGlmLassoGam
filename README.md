# RdeGlmLassoGam
R-package accompanying the paper "Robust Inference and Modeling of Mean and Dispersion for Generalized Linear Models"
This package provides a robust estimator for the parameters in the double exponential family, the so-called Robust Double Exponential (RDE) Estimator. 
One can use (a combination of) Generalised Linear Models, Generalised Penalized Linear Models or Generalised Additive Models.

Required libraries:
  "CompQuadForm
  "Matrix"
  "robustbase"
  "MASS"
  "stats"
  "utils"
  "cellWise"
  "lassoshooting"
  "glmnet"
  "tidyr"
  
The package consists several files that can briefly be explained as follows: 
  - DBinom_Functions.R contains the density, distribution function, quantile function and random generation for the double binomial distribution with parameters n, mu and theta.
	'dDBinom', 'pDBinom', 'qDBinom', 'rDBinom'
  - DPois_Functions.R contains the density, distribution function, quantile function and random generation for the double Poisson distribution with parameters mu and theta.
	'dDPois', 'pDPois', 'qDPois', 'rDPois'
  - RDE.R contains the main function of this package, called 'RDE' (together with some internal help functions).
	It estimates the RDE coefficients.
  - DispersionFiles.R contains the exported function 'DispersionTest' (together with some internal help functions). 
	It can be used to compare two nested models.
  - glmRDENoZ.R contains the exported function 'glmRDENoZ' (together with some internal help functions).
	It estimates the RDE coefficients in case there is no dispersion part.
  - glmRDENoZ.R contains the exported function 'glmRDENoZ' (together with some internal help functions).
	It estimates the RDE coefficients in case there is no dispersion part.
  - glmRDENoZ.R contains the exported function 'glmRDENoZ' (together with some internal help functions).
	It estimates the RDE coefficients in case there is no dispersion part.
  - CalcAMSE.R contains the exported function 'calcAMSE'. This function calculates the AMSE of the RDE estimator.
  - CalcTuningParam.R contains the exported function 'CalcTuningParam' (together with an internal help function).
	It determines the smallest tuning parameter for the RDE estimator such that the desired AMSE is obtained.
  - EBICBet.R (EBICGam.R) contains the exported function 'EBICBet'('EBICGam').
	It estimates the EBIC value for the mean (dispersion) parameters. 
  - RBICBet.R (RBICGam.R) contains the exported function 'RBICBet'('RBICGam').
  	It estimates the RBIC value for the mean (dispersion) parameters.
  - fit.gam.R, glmrob.R and glmrobMqle.R contain some internal help functions used to obtain good starting values for the RDE estimator. 
  - robustgam.basic.R  contains some internal functions that are used to determine the basis functions for the GAM RDE estimator.

A detailed explanation of the exported functions can be found in their help file.

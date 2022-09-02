##' This function determines the smallest tuning parameter for the RDE estimator such that the desired AMSE is obtained.
##' @param cValues the vector with values for the tuning parameter that should be checked.
##' @param AMSECrit the critical AMSE value.
##' @param y the response vector.
##' @param X the model matrix for the mean.
##' @param Z the model matrix for the dispersion.
##' @param m in case family equals "binomial", this parameter should represent a vector with the number of trials. Default value is NULL.
##' @param family a character string indicating the family. This can be "poisson" or "binomial".
##' @param weights.on.xz1 a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the mean. It is also possible to provide a character string. In case this is "none", all observations get weight 1. In case this is "covMcd", the weights are determined via the function robustbase::covMcd..
##' @param weights.on.xz2 a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the dispersion. It is also possible to provide a character string. In case this is "none", all observations get weight 1. In case this is "covMcd", the weights are determined via the function robustbase::covMcd. Default value is NULL, meaning that the same weights as for the dispersion model are used.
##' @param weightFunction1 a character string indicating which weight function is used to diminish the effect of large residuals in the model for the mean. This can be "Huber" or "Tukey". Default value is "Huber".
##' @param weightFunction2 a character string indicating which weight function is used to diminish the effect of large residuals in the model for the dispersion. This can be "Huber" or "Tukey". Default value is NULL, meaning that the same function as for the mean model is used.
CalcTuningParam <- function(cValues, AMSECrit,
                               y, X, Z, m, family,
                               weights.on.xz1 = "none", weights.on.xz2=NULL,
                               weightFunction1, weightFunction2 = NULL){
  functionforapply <- function(cValue,y,Xmod,Z, m, family,
                               weights.on.xz1, weights.on.xz2,
                               weightFunction1, weightFunction2){
    model = RDE(y=y, X=Xmod, Z=Z, mBin = m,
                           family=family,
                           weights.on.xz1 = weights.on.xz1,
                           weights.on.xz2 = weights.on.xz2,
                           weightFunction1 = weightFunction1,
                           weightFunction2 = weightFunction2,
                           optionList = list(huberC = cValue,
                                             tukeyC = cValue,
                                             tol = 1e-4,
                                             maxIt = 100
                           ))
    AMSEVal = CalcAMSE(y=y, X=Xmod, Z=Z, family=family,
                       m = m,
                       weights.xz1 = model$weights.xz1, weights.xz2 = model$weights.xz2,
                       alpha = c(model$betaEstimate, model$gammaEstimate),
                       cutOff1 = model$c1, weightF1 = model$v_fun1,
                       cutOff2 = model$c2, weightF2 = model$v_fun2)

    return(c(cValue,AMSEVal[1], AMSEVal[2],AMSEVal[3]))
  }
  output = sapply(X=cValues,FUN=functionforapply,y=y,Xmod=X,Z=Z, m=m, family=family,weights.on.xz1 = weights.on.xz1, weights.on.xz2 = weights.on.xz2,
                  weightFunction1 = weightFunction1, weightFunction2 = weightFunction2)
  cBeta <- cValues[which(output[2,]>AMSECrit)[1]]
  cGamma <- cValues[which(output[3,]>AMSECrit)[1]]
  cBetaGamma <- cValues[which(output[4,]>AMSECrit)[1]]
  return(list(tunPar = c(cBeta,cGamma,cBetaGamma),AMSEBeta = output[2,],AMSEGamma = output[3,], AMSEAlpha = output[4,]))
}



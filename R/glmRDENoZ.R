##'This function estimates the RDE coefficients in case there is no dispersion part.
##' @param y the response vector.
##' @param X the model matrix for the mean.
##' @param m in case family equals "binomial", this parameter should represent a vector with the number of trials. Default value is NULL.
##' @param family a character string indicating the family. This can be "poisson" or "binomial".
##' @param beta.ini optional initial value for the beta coefficients.
##' @param weights.on.xz1 a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the mean. It is also possible to provide a character string. In case this is "none", all observations get weight 1. In case this is "covMcd", the weights are determined via the function robustbase::covMcd. The tuning parameter alpha can be provided in the optionList.
##' @param contX indices of the columns in X representing the the continuous variables.  When weights.on.xz1 is equal to "covMcd", this parameter should be specified. Default value is NULL.
##' @param weightFunction1 a character string indicating which weight function is used to diminish the effect of large residuals in the model for the mean. This can be "Huber" or "Tukey". Default value is "Huber".
##' @param optionList list that should contain the tuning parameter huberC in case weightFunction equals "Huber", or the tuning parameter tukeyC in case weightFunction equals "Tukey". Furthermore, the tuning parameter alpha for the function robustbase::covMcd can be provided, with default value 0.7 . Finally, maxIt which is the maximum number of iterations and tol, which is the tolerance can be provided here. Default value is list(huberC = 2, tukeyC = 3, tol = 1e-4, maxIt = 100, alpha = 0.75).
glmRDENoZ <- function(y, X, m = NULL,
                               family,
                               beta.ini,
                               weights.on.xz1 = "none",
                               contX = NULL,
                               weightFunction1 = "Huber",
                               optionList = list(huberC = 2,
                                                 tukeyC = 3,
                                                 tol = 1e-4,
                                                 maxIt = 100,
                                                 alpha = 0.75
                               ))
{

  #################
  # General setup #
  #################

  n <- length(y)
  X <- as.matrix(X)
  if(family == "binomial"){
    if(is.null(m)){
      stop('Please provide the variable m.')
    }
  }

  p = ncol(X)

  ###########
  # Weights #
  ###########
  if(is.numeric(weights.on.xz1)){
    if(length(weights.on.xz1) == n & !(any(weights.on.xz1 < 0)) & !(any(weights.on.xz1 > 1))){
      weights.xz1 <- weights.on.xz1
    }else{
      stop("All 'weights.on.xz1' must be non-negative and smaller or equal to 1")
    }
  }else{
    if(weights.on.xz1 == "none"){
      weights.xz1 = rep(1,n)
    }else{
      if(weights.on.xz1 == "covMcd"){
        if(!("alpha" %in% names(optionList))){
          optionList[["alpha"]] <- 0.75
        }
        if(!is.null(contX)){
          weights.xz1 <-  covMcd(X[,contX], alpha = optionList$alpha)$mcd.wt
        }else{
          stop("When 'weights.on.xz1' is equal to 'covMcd', 'contX' should be specified.")
        }
      }else{
        stop("This value for 'weights.on.xz1' is undefined.")
      }
    }
  }

  if (weightFunction1 == "Huber") {
    cutOff1 <- optionList$huberC
    weightF1 <- function(Arg1) (HuberResid(Arg1, cutOff1) / Arg1)^2
  } else if (weightFunction1 == "Tukey") {
    cutOff1 <- optionList$tukeyC
    weightF1 <- function(Arg1) TukeyResid(Arg1, cutOff1) / Arg1
  }

  #####################
  # Initial estimates #
  #####################
  deltaOld <- rep(0.0, p )
  if(weightFunction1 == "Huber"){
    if (missing(beta.ini)){
      if(family=="binomial"){
        yStart <- cbind(m*y, m-m*y)
      }else{
        yStart <- y
      }

      if(p==1){
        glmrob.ini <- glmrob(yStart ~ 1, method = "Mqle", family = family,
                             control = glmrobMqle.control(acc = 0.01*optionList$tol,
                                                          maxit = 50*optionList$maxIt,
                                                          tcc = optionList$huberC),
                             weights.on.x = weights.xz1,
                             data = data.frame(y, X))
      } else {
        glmrob.ini <- glmrob(yStart~X[,-1],
                             method = "Mqle",
                             family = family,
                             weights.on.x =weights.xz1 ,
                             control = glmrobMqle.control(acc = 0.01*optionList$tol,
                                                          maxit = 50*optionList$maxIt,
                                                          tcc = optionList$huberC))

      }
      betaOld <- as.numeric(glmrob.ini$coefficients)
      deltaOld[1:p] <- betaOld
    } else {
      deltaOld[1:p] <- beta.ini
    }
  }else{
    if (missing(beta.ini)){
      if(family=="binomial"){
        yStart <- cbind(m*y, m-m*y)
      }else{
        yStart <- y
      }

      if(p==1){
        glmrob.ini <- glmrob(yStart ~ 1, method = "Mqle", family = family,
                             control = glmrobMqle.control(acc = 0.0001*optionList$tol,
                                                          maxit = 500*optionList$maxIt,
                                                          tcc = optionList$huberC),
                             weights.on.x = weights.xz1,
                             data = data.frame(y, X))
      } else {
        glmrob.ini <- glmrob(yStart~X[,-1],
                             method = "Mqle",
                             family = family,
                             weights.on.x =weights.xz1 ,
                             control = glmrobMqle.control(acc = 0.0001*optionList$tol,
                                                          maxit = 500*optionList$maxIt,
                                                          tcc = optionList$huberC))

      }
      betaOld <- as.numeric(glmrob.ini$coefficients)
      deltaOld[1:p] <- betaOld
    } else {
      deltaOld[1:p] <- beta.ini
    }
  }
  deltaTemp <- deltaOld
  deltaStart <- deltaOld

  #####################
  # Estimation        #
  #####################

  i <- 1
  while(i <= optionList$maxIt){

    deltaOld <- deltaTemp
    conv_B <- FALSE
    # Update beta
    maxIt.beta <- 4 * optionList$maxIt
    deltaTemp2 <- deltaTemp
    # browser()
    for(j in 1:maxIt.beta){
       betaUpdate <- try(
        CalculateBetaUpdateNOZ2(y = y, X = X,
                                family = family,
                                m = m,
                                alpha = deltaTemp2,
                                cutOff = cutOff1, weightF = weightF1,
                                weights.xz = weights.xz1,
                                optionList = optionList),
        silent = TRUE)
      if(!is.matrix(betaUpdate)){
        deltaTemp2 <- deltaTemp
        break()
      }
      if (any(!is.finite(betaUpdate))) {
        warning("Non-finite coefficients at iteration ", i)
        deltaTemp2 <- deltaTemp
        break()
      }
      relE_B1 <- sqrt(sum(betaUpdate^2)/max(1e-20, sum(deltaTemp[1:p]^2)))
      relE_B2 <- max(abs(betaUpdate))
      conv_B <- (relE_B1 <= optionList$tol)|(relE_B2 <= optionList$tol)
      deltaTemp2[1:p] <- deltaTemp2[1:p] + betaUpdate
      if (conv_B) break()
    }
    deltaTemp <- deltaTemp2
    relE <- sqrt(sum((deltaOld - deltaTemp)^2/deltaTemp^2))
    conv <- relE <= optionList$tol
    if (conv) i <- optionList$maxIt + 1
    i <- i + 1

  }

  betaEstimate <- deltaTemp[1:p]
  EstTarget <- deltaTemp

  ASInfo <- CalculateAsymptoticInfoNoZ2(y = y, X = X, m = m, delta = deltaTemp,
                                    family = family,
                                    cutOff1 = cutOff1,
                                    weightF1 = weightF1,
                                    weights.xz1 = weights.xz1,
                                    optionList = optionList)

  Result <- list(betaEstimate = betaEstimate,
                 thetaStart = deltaStart,
                 cutOff1 = cutOff1,
                 weightF1 = weightF1,
                 weights.xz1 = weights.xz1,
		             M=ASInfo$M,
                 Q=ASInfo$Q,
		             ASVar= ASInfo$ASVar
  )
}

# This function calculates asymptotic information: M matrix, Q matrix and asymptotic variance matrix.
##' @param y the response vector.
##' @param X the model matrix for the mean.
##' @param m in case family equals "binomial", this parameter should represent a vector with the number of trials. Default value is NULL.
##' @param delta the beta parameter estimates.
##' @param family a character string indicating the family. This can be "poisson" or "binomial".
##' @param cutOff1 tuning parameter of the mean model.
##' @param weightF1 the weight function that is used to diminish the effect of large residuals in the model for the mean. 
##' @param weights.xz1 a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the mean.
##' @param optionList list that should contain the tuning parameter huberC in case weightFunction equals "Huber", or the tuning parameter tukeyC in case weightFunction equals "Tukey". Furthermore, the tuning parameter alpha for the function robustbase::covMcd can be provided, with default value 0.7 . Finally, maxIt which is the maximum number of iterations and tol, which is the tolerance can be provided here.
CalculateAsymptoticInfoNoZ2 <- function(y, X, m, delta,
                                    family,
                                    cutOff1, weightF1,
                                    weights.xz1,
                                    optionList){
  n <- length(y)
  p = ncol(X)

  if(family == "binomial"){
    if(is.null(m)){
      stop('Please provide the variable m.')
    }
  }

  beta <- delta[1:p]
  gamma = 0
  # Calculate the AMSE
  if (family == "poisson") {
    mu <- as.numeric(exp(X %*% beta))
    theta <- rep(1,n)

    # Limit the number of considered observations considerably
    minValue <- pmin(pmax(0, floor(mu - cutOff1 * (1 / theta) * mu)))
    maxValue <- pmax(ceiling(mu + cutOff1 * (1/theta) * mu))
    tol <- 1e-12
    minCut <- qpois(p = tol, lambda = mu)
    maxCut <- qpois(p = 1 - tol, lambda = mu)
    minValue <- pmax(minValue, minCut)
    maxValue <- pmin(maxValue, maxCut)

    BnumMat = 0
    BdenMat1 = 0
    BdenMat2 = 0
    BdenMat3 = 0
    GdenMat3 = 0
    BGnumMat12 = 0
    M12 =0
    M21 =0
    Q12 =0
    Q21 =0
    a1 = 0
    a2 = 0

    for (i in 1:n) {
      j <- minValue[i]:maxValue[i]
      Probs <- dDPois(j, mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / (sqrt(mu[i] / theta[i]))
      weight1 <- weightF1(R)
      weight1[which(R == 0)] <- 1

      UMu <- (j-mu[i])/(mu[i]/theta[i])
      UTheta <- 1 / (2 * theta[i]) - mu[i] + j * log(exp(1) * mu[i] / j)
      if(j[1]==0){
        UTheta[1] <- 1 / (2 * theta[i]) - mu[i]
      }
      xxt = tcrossprod(X[i,])
      BnumMat = BnumMat + sum(UMu^2 * Probs)*mu[i]^2*xxt
      BdenMat1 = BdenMat1 + sum(weight1 * UMu^2 * Probs)*weights.xz1[i]*mu[i]^2*xxt
      BdenMat2 = BdenMat2 + sum(UMu * Probs)*mu[i]*t(X[i,])
      BdenMat3 = BdenMat3 + sum(weight1^2 * UMu^2 * Probs)*weights.xz1[i]^2*mu[i]^2*xxt
      a1 = a1 + sum(weight1 * UMu * Probs) * weights.xz1[i] * mu[i] * X[i,]
    }
    Q = matrix(0,p,p); M = matrix(0,p,p)
    Q[1:p,1:p]=BdenMat3-1/n*tcrossprod(a1);
    M[1:p,1:p]=BdenMat1-1/n*a1%*%BdenMat2;
    Q = Q/n
    M = M/n
    Minv = solve(M)
    ASVar = Minv %*% Q %*% t(Minv)
  }else if(family == "binomial"){
    linpred <- as.numeric(X %*% beta)
    mu <- as.numeric(1 / (1 + exp(-linpred)))
    theta <- rep(1,n)
    der <- exp(-linpred)/(1+exp(-linpred))^2
    # Limit the number of considered observations considerably
    tol <- 1e-12
    minValue <- qbinom(p = tol, size = m, prob = mu)
    maxValue <- qbinom(p = 1 - tol, size = m, prob = mu)

    BnumMat = 0
    BdenMat1 = 0
    BdenMat2 = 0
    BdenMat3 = 0
    M12 =0
    M21 =0
    Q12 =0
    Q21 =0
    a1 = 0
    a2 = 0

    for (i in 1:n) {
      j <- minValue[i]:maxValue[i] / m[i]
      Probs <- dDBinom(j, n = m[i] , mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / sqrt(mu[i] * (1 - mu[i]) / (m[i] * theta[i]))
      weight1 <- weightF1(R)
      weight1[which(R == 0)] <- 1

      UMu <- (m[i] * j * theta[i]) / mu[i] + (theta[i] * m[i] * (1 - j)) / (mu[i] - 1)
      UTheta <- 1 / (2 * theta[i]) + (m[i] * j) * log(mu[i] / j) + m[i] * (1 - j) * log((1 - mu[i]) / (1 - j))

      if(j[1] <= 1e-12){
        UMu[1] <- m[i]*theta[i]/(mu[i]-1)
        UTheta[1] <- 1 / (2 * theta[i]) + m[i] * log(1 - mu[i])
      }
      if(j[length(j)] >= 1 - 1e-12){
        UMu[length(j)] <- m[i]*theta[i]/mu[i]
        UTheta[length(j)] <- 1 / (2 * theta[i]) + m[i] * log(mu[i])
      }
      xxt = tcrossprod(X[i,])

      BnumMat = BnumMat + sum(UMu^2 * Probs)*der[i]^2*xxt
      BdenMat1 = BdenMat1 + sum(weight1 * UMu^2 * Probs)*weights.xz1[i]*der[i]^2*xxt
      BdenMat2 = BdenMat2 + sum(UMu * Probs)*der[i]*t(X[i,])
      BdenMat3 = BdenMat3 + sum(weight1^2 * UMu^2 * Probs)*weights.xz1[i]^2*der[i]^2*xxt
      a1 = a1 + sum(weight1 * UMu * Probs) * weights.xz1[i] * der[i] * X[i,]
    }
    Q = matrix(0,p,p); M = matrix(0,p,p)
    Q[1:p,1:p]=BdenMat3-1/n*tcrossprod(a1);
    M[1:p,1:p]=BdenMat1-1/n*a1%*%BdenMat2;
    Q = Q/n
    M = M/n

    Minv = solve(M)
    ASVar = Minv %*% Q %*% t(Minv)
  }else{
    stop("Family not yet implemented.")
  }
  return(list(M = M,
              Q = Q,
              ASVar = ASVar))
}

# This function is used to update the beta parameter.
##' @param y the response vector.
##' @param X the model matrix for the mean.
##' @param family a character string indicating the family. This can be "poisson" or "binomial".
##' @param m in case family equals "binomial", this parameter should represent a vector with the number of trials.
##' @param alpha a vector with estimated beta values
##' @param cutOff tuning constant used in the mean model
##' @param weightF the weight function is used to diminish the effect of large residuals in the model for the mean. 
##' @param weights.xz It is a numeric vector, specifying how points (potential outliers) in xz-space are downweighted while modelling the mean. 
##' @param optionList list that should contain the tuning parameter huberC in case weightFunction equals "Huber", or the tuning parameter tukeyC in case weightFunction equals "Tukey". Furthermore, the tuning parameter alpha for the function robustbase::covMcd can be provided, with default value 0.75 . Finally, maxIt which is the maximum number of iterations and tol, which is the tolerance can be provided here. 
CalculateBetaUpdateNOZ2 <- function(y, X,
                                    family,
                                    m,
                                    alpha,
                                    cutOff, weightF,
                                    weights.xz,
                                    optionList){

  # Data matrix X rows correspond to observations
  n <- length(y)
  p = ncol(X)

  beta <- alpha[1:p]
  gamma <- 0

  if (family == "poisson") {
    mu <- as.numeric(exp(X %*% beta))
    theta <- rep(1,n)
    dLogMu <- (y - mu) / (mu / theta)
    dMuBeta <- X * mu # contains at row i the vector mu * X[i,]
    dLogTheta <- 1 / (2 * theta) - mu + y * log(exp(1) * mu / y)
    tInd <- which(y == 0)
    if (length(tInd) > 0) dLogTheta[tInd] <- 1 / (2 * theta[tInd]) - mu[tInd]
    PearsonResid <- (y - mu) / sqrt(mu / theta)
  } else if (family == "binomial") {
    linpred <- as.numeric(X %*% beta)
    mu <- as.numeric(exp(linpred) / (1 + exp(linpred)))
    theta <- rep(1,n)
    dLogMu <- (m * y * theta) / mu + (theta * m * (1 - y)) / (mu - 1)
    dMuBeta <- X * (mu / (1 + exp(linpred))) # contains at row i the vector vec[i] * X[i,]
    dLogTheta <- 1 / (2 * theta) + (m * y) * log(mu / y) + m * (1 - y) * log((1 - mu) / (1 - y))

    #Adapt for boundary values
    tInd <- which(y <= 1e-12)
    if(length(tInd)>0) {
      dLogTheta[tInd] <- 1/(2*theta[tInd]) + m[tInd] * log(1 - mu[tInd])
    }
    tInd <- which(y >= 1 - 1e-12)
    if(length(tInd)>0) {
      dLogTheta[tInd] <- 1/(2*theta[tInd]) + m[tInd] * log(mu[tInd])
    }

    # PearsonResid <- (y - mu) / sqrt(mu * (1 - mu) * theta / n)
    PearsonResid <- (y - mu) / sqrt(mu * (1 - mu) / (m * theta))
  } else {
    stop("Requested family not yet implemented.")
  }

  # Calculate the expectation
  ExpTermBeta <- rep(0.0, p)
  B11 <- rep(0.0, n)
  if (family == "poisson") {
    # Limit the number of considered observations considerably
    minValue <- pmax(0, floor(mu - cutOff * (1 / theta) * mu))
    maxValue <- ceiling(mu + cutOff * (1/theta) * mu)
    tol <- 1e-12
    minCut <- qpois(p = tol, lambda = mu)
    maxCut <- qpois(p = 1 - tol, lambda = mu)
    minValue <- pmax(minValue, minCut)
    maxValue <- pmin(maxValue, maxCut)

    for (i in 1:n) {
      j <- minValue[i]:maxValue[i]
      Probs <- dDPois(j, mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / (sqrt(mu[i] / theta[i]))
      LLTerm <- R / (sqrt(mu[i]/theta[i]))
      weight <- weightF(R)
      weight[which(R == 0)] <- 1

      ExpTermBeta <- ExpTermBeta + sum(weight * LLTerm * Probs) * (weights.xz[i] * as.matrix(dMuBeta[i,]))
      B11[i] <- weights.xz[i] * sum( weight * LLTerm ^ 2 * Probs) * mu[i] ^ 2
    }

    ExpTermBeta = (1 / n) * ExpTermBeta
    B11 = diag(B11)
  } else if (family == "binomial") {
    # Limit the number of considered observations considerably
    tol <- 1e-12
    minValue <- qbinom(p = tol, size = m, prob = mu)
    maxValue <- qbinom(p = 1 - tol, size = m, prob = mu)
    for (i in 1:n) {
      j <- minValue[i]:maxValue[i] / m[i]
      Probs <- dDBinom(j, n = m[i] , mu = mu[i], theta = theta[i], correction = FALSE)
      R <- (j - mu[i]) / sqrt(mu[i] * (1 - mu[i]) / (m[i] * theta[i]))
      LLTerm <- (m[i] * j * theta[i]) / mu[i] + (theta[i] * m[i] * (1 - j)) / (mu[i] - 1)
      if(j[1] <= tol) LLTerm[1] <- m[i] * theta[i] / (mu[i] - 1)
      if(j[length(j)] >= 1-tol) LLTerm[length(j)] <- m[i] * theta[i] / mu[i]
      weight <- weightF(R)
      weight[which(R == 0)] <- 1

      ExpTermBeta <- ExpTermBeta + sum(weight * LLTerm * Probs) * (weights.xz[i] * as.matrix(dMuBeta[i,]))
      B11[i] <- weights.xz[i] * sum( weight * LLTerm ^ 2 * Probs) * mu[i] ^ 2
    }
    ExpTermBeta = (1 / n) * ExpTermBeta
    B11 = diag(B11)
  } else {
    stop("Requested family not yet implemented.")
  }

  weight <- weightF(PearsonResid)
  weight[which(PearsonResid == 0)] <- 1

  # Calculate the LogLikelihoodValue
  LLBeta <- rep(0.0, p)
  # for (i in 1:n) {
  #   LLBeta <- LLBeta + weight[i] * dLogMu[i] * weights.xz1[i] * as.matrix(dMuBeta[i,])
  # }
  LLBeta <- as.matrix(colSums(dMuBeta * (weight * dLogMu * weights.xz))) # matrix expression for loop above
  LLBeta <- LLBeta - n * ExpTermBeta

  # Calculate the PsiMatrix
  PsiDeriv <- crossprod(X, B11) %*% X
  PsiDeriv <- as.matrix(PsiDeriv)
  return(solve(PsiDeriv,diag(p)) %*% LLBeta)

}

IF.SoR.mean <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsSoR_M.IF=list(), k=4,
                        IFplot=FALSE, IFprint=TRUE,
                        rf=0, compile=TRUE, prewhiten=FALSE, ar.prewhiten.order=1,
                        cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                        ...){
  
  # Function evaluation
  if(isTRUE(evalShape)){
    if(is.null(retVals))
      if(!is.null(returns))
        retVals <- seq(mean(returns)-k*sd(returns), mean(returns)+k*sd(returns), by=0.001) else
          retVals <- seq(0.005-k*0.07, 0.005+k*0.07, by=0.001)
        IFvals <- cbind(retVals, IF.fn(retVals, risk="SoR", returns=returns, parsSoR.IF=parsSoR_M.IF, rf=rf, threshold="mean"))
        colnames(IFvals) <- c("r", "IFvals")
        if(isTRUE(IFplot)){
          plot(IFvals[,1], IFvals[,2], type="l", 
               xlab="r", ylab="IF", col="blue", lwd=1, 
               main="Influence Function - SoRmu",
               panel.first=grid(), cex.lab=1.25)
          abline(h=0, v=0)
        }
        if(IFprint)
          return(IFvals) else{
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt)) 
            stop() 
          }
  }
  
  # Storing the dates
  if(xts::is.xts(returns))
    returns.dates <- zoo::index(returns)
  
  # Adding the robust filtering functionality
  if(cleanOutliers){
    temp.returns <- robust.cleaning(returns, cleanMethod, alpha.robust, eff)
    if(xts::is.xts(returns))
      returns <- xts::xts(temp.returns, returns.dates) else
        returns <- temp.returns
  }
  
  if(compile){
    IF.SoR_M.vector <- as.vector(IF_SoR_mean(returns, rf))
  } else{
    # Computing the mean of the returns
    mu.hat <- mean(returns)
    # Computing the SD of the returns
    sigma.hat <- sqrt(mean((returns-mu.hat)^2))
    # Computing SD- of the returns
    sigma.minus.hat <- sqrt(mean((returns-mu.hat)^2*(returns<=mu.hat)))
    # Computing the SoR estimate
    SoR.hat <- (mu.hat-rf)/sigma.minus.hat
    # Computing the mean parameter for the SoR
    mu1.minus.hat <- mean((returns-mu.hat)*(returns<=mu.hat))
    
    # Computing the IF vector for SoR_M
    IF.SoR_M.vector <- -SoR.hat/2/sigma.minus.hat^2*(returns-mu.hat-rf)^2*(returns<=mu.hat)+
      (1/sigma.minus.hat+SoR.hat*mu1.minus.hat/sigma.minus.hat^2)*(returns-mu.hat-rf)+
      SoR.hat/2
  }
  
  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.SoR_M.vector <- as.numeric(arima(x=IF.SoR_M.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)
  
  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.SoR_M.vector <- xts::xts(IF.SoR_M.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.SoR_M.vector, type="l", main="SoR (Mean Threshold) Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the IF vector for SoR_M
  if(xts::is.xts(returns))
    return(xts::xts(IF.SoR_M.vector, returns.dates)) else
      return(IF.SoR_M.vector)
}

IF.SoR.const <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsSoR_C.IF=list(), k=4,
                         IFplot=FALSE, IFprint=TRUE,
                         const=0, prewhiten=FALSE, ar.prewhiten.order=1,
                         cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                         ...){
  
  # Function evaluation
  if(isTRUE(evalShape)){
    if(is.null(retVals))
      if(!is.null(returns))
        retVals <- seq(mean(returns)-k*sd(returns), mean(returns)+k*sd(returns), by=0.001) else
          retVals <- seq(0.005-k*0.07, 0.005+k*0.07, by=0.001)
        IFvals <- cbind(retVals, IF.fn(x=retVals, risk="SoR", returns=returns, parsSoR.IF=parsSoR_C.IF, const=const, threshold="const"))
        colnames(IFvals) <- c("r", "IFvals")
        if(isTRUE(IFplot)){
          plot(IFvals[,1], IFvals[,2], type="l", 
               xlab="r", ylab="IF", col="blue", lwd=1, 
               main="SoRc",
               panel.first=grid(), cex.lab=1.25)
          abline(h=0, v=0)
        }
        if(IFprint)
          return(IFvals) else{
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt)) 
            stop() 
          }
  }
  
  # Storing the dates
  if(xts::is.xts(returns))
    returns.dates <- zoo::index(returns)
  
  # Adding the robust filtering functionality
  if(cleanOutliers){
    temp.returns <- robust.cleaning(returns, cleanMethod, alpha.robust, eff)
    if(xts::is.xts(returns))
      returns <- xts::xts(temp.returns, returns.dates) else
        returns <- temp.returns
  }
  
  # Length of the vector of returns
  N <- length(returns)
  # Paremeter of SoR_C
  mar.parameter <- const
  # Mean of the returns
  mu.hat <- mean(returns)
  # Computing SD for the returns
  sigma.hat <- sqrt(mean((returns-mu.hat)^2))
  # Computing SD- of the returns
  sigma.minus.hat <- sqrt(sum((returns-mar.parameter)^2*(returns<=mar.parameter))/N)
  # Computing the SoR estimate
  SoR.hat <- (mu.hat-mar.parameter)/sigma.minus.hat
  
  # Computing the IF vector for SoR_C
  IF.SoR_C.vector <- -SoR.hat/2/sigma.minus.hat^2*(returns-mar.parameter)^2*(returns<=mar.parameter)+
    1/sigma.minus.hat*(returns-mu.hat)+
    SoR.hat/2
  
  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.SoR_C.vector <- as.numeric(arima(x=IF.SoR_C.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)
  
  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.SoR_C.vector <- xts::xts(IF.SoR_C.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.SoR_C.vector, type="l", main="SoR (Constant Threshold) Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the IF vector for SoR_C
  if(xts::is.xts(returns))
    return(xts::xts(IF.SoR_C.vector, returns.dates)) else
      return(IF.SoR_C.vector)
}

#' @title Influence Function - Sortino Ratio
#' 
#' @description \code{IF.SoR} returns the data and plots the shape of either the IF or the IF TS for the Sortino Ratio
#'
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsSoR.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param threshold Parameter of threshold is either "mean" or "const". Default is "mean".
#' @param const The threshold if threshold is "const". 
#' @param rf Risk-free interest rate.
#' @param compile Boolean variable to indicate if the IF TS should be computed using compiled code (C++) (TRUE) or not (FALSE).
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm.
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Addtional parameters.
#'
#' @return Influence function of SoR_C.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @export
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.SoR(returns=NULL, evalShape=TRUE, 
#'                 retVals=NULL, parsSoR.IF=list(mu=0.01, lpm2=0.00898, sor.c=0.3337, 
#'                                               ssd=0.0354, smean=-0.0199, sor.mu=0.2929),
#'                  IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.SoR(risk="mean",
#'                  returns=edhec[,"CA"], evalShape=TRUE, 
#'                  retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                  IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.SoR(returns=edhec[,"CA"], evalShape=FALSE, 
#'                 retVals=NULL, nuisance.par=NULL,
#'                 IFplot=TRUE, IFprint=TRUE,
#'                 compile=TRUE, prewhiten=FALSE,
#'                 cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF.SoR <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsSoR.IF=list(mu=0.01, lpm2=0.00898, sor.c=0.3337, 
                                                                                ssd=0.0354, smean=-0.0199, sor.mu=0.2929), k=4,
                   IFplot=FALSE, IFprint=TRUE,
                   compile=TRUE, threshold=c("mean", "const")[1], const=0, rf=0, prewhiten=FALSE, ar.prewhiten.order=1,
                   cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                   ...){
  
  # Checking data if IF risk evaluation
  if(isTRUE(evalShape))
    if(is.null(returns))
      if(is.null(parsSoR.IF))
        stop("For shape evaluation, nuisance parameters must be specified if no returns are provided.")
  
  # Checking data if IF TS evaluation
  if(!isTRUE(evalShape))
    if(is.null(returns))
      stop("Returns must be provided for the IF TS evaluation.")
  
  # Checking the data for the returns
  if(!is.null(returns)){
    if(!any(c(inherits(returns, "matrix"), inherits(returns, "numeric"), inherits(returns, "xts"), inherits(returns, "zoo")))){
      stop("returns should belong to one of the following classes: matrix, numeric, xts, zoo")
    } else if(any(anyNA(returns), any(is.nan(returns)), any(is.infinite(returns)))){
      stop("returns should not have missing, infinite or nan values")
    } else{
      if(inherits(returns, "matrix")){
        if(ncol(returns)>1){
          stop("returns should be a vector")
        }
        # Force to vector if input was a matrix
        returns <- as.numeric(returns)
      }
    }
  }
  
  # Checking the data for k (range parameter)
  if(!inherits(k, "numeric")){
    stop("k should be numeric")
  } else if(any(!k == floor(k), k <= 0)){
    stop("k should be a positive integer")
  }
  
  # Checking data for prewhitening order
  if (!inherits(ar.prewhiten.order, "numeric")) {
    stop("ar.prewhiten.order should be numeric")
  } else if (any(!ar.prewhiten.order == floor(ar.prewhiten.order), ar.prewhiten.order <= 0)) {
    stop("ar.prewhiten.order should be a positive integer")
  }
  
  # Checking data for const value
  if(!inherits(const, "numeric")){
    stop("const should be numeric")
  }
  
  # Checking data for rf value
  if(!inherits(rf, "numeric")){
    stop("rf should be numeric")
  } else if(any(rf < 0, rf > 1)) {
    stop("rf should be a numeric value between 0 and 1.")
  }
  
  # Checking robust cleaning method specified
  if(!(cleanMethod %in% c("locScaleRob", "Boudt")))
    stop("The specified outlier cleaning method is not available.")
  
  # Checking data for the efficiency for robust cleaning with "locScaleRob"
  if(!inherits(eff, "numeric")){
    stop("eff should be numeric")
  } else if(any(eff < 0, eff > 1)) {
    stop("eff should be a numeric value between 0 and 1.")
  }
  
  # Checking data for the efficiency for robust cleaning with "Boudt"
  if(!inherits(alpha.robust, "numeric")){
    stop("alpha.robust should be numeric")
  } else if(any(alpha.robust < 0, alpha.robust > 1)) {
    stop("alpha.robust should be a numeric value between 0 and 1.")
  }
  
  # Function evaluation for mean threshold
  if(threshold=="mean")
    IF.SoR.mean(returns=returns, evalShape=evalShape, retVals=retVals, parsSoR_M.IF=parsSoR.IF, 
                IFplot=IFplot, IFprint=IFprint,
                rf=rf, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, eff=eff, alpha.robust=alpha.robust,
                ...) else if(threshold=="const")
                  IF.SoR.const(returns=returns, evalShape=evalShape, retVals=retVals, parsSoR_C.IF=parsSoR.IF, 
                               IFplot=IFplot, IFprint=IFprint,
                               const=const, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                               cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, eff=eff, alpha.robust=alpha.robust,
                               ...)
}
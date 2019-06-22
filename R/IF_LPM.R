# Function to compute the LPM
LPM <- function(returns, const = 0, order = 1, ...){
  
  # Compute the length of the returns vector
  N <- length(returns)
  
  # Computing the LPM
  return(1/N*sum((const-returns[returns<=const])^order)^(1/order))
}

#' @title Influence Function - Lower Partial Moment (LPM)
#' 
#' @description \code{IF.LPM} returns the data and plots the shape of either the IF or the IF TS for the LPM
#'
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsLPM.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param const Constant threshold.
#' @param order Order of LPM. Can only take values 1 or 2.
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function of LPM.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @export
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.LPM(returns=NULL, evalShape=TRUE, 
#'                 retVals=NULL, parsLPM.IF=list(lpm1=0.0255, lpm2=0.00218),
#'                 IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.LPM(risk="mean",
#'                 returns=edhec[,"CA"], evalShape=TRUE, 
#'                 retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                 IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.LPM(returns=edhec[,"CA"], evalShape=FALSE, 
#'                 retVals=NULL, nuisance.par=NULL,
#'                 IFplot=TRUE, IFprint=TRUE,
#'                 prewhiten=FALSE,
#'                 cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF.LPM <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsLPM.IF=list(lpm1=0.0255, lpm2=0.00218), k=4,
                   IFplot=FALSE, IFprint=TRUE,
                   const=0, order=1, prewhiten=FALSE, ar.prewhiten.order=1,
                   cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                   ...){
  
  # Checking data if IF risk evaluation
  if(isTRUE(evalShape))
    if(is.null(returns))
      if(is.null(parsLPM.IF))
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
  
  # Checking the data for k (range parameter)
  if(!inherits(order, "numeric")){
    stop("order should be numeric")
  } else if(any(!order == floor(order), order <= 0)){
    stop("order should be a positive integer")
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
  
  # Function evaluation
  if(isTRUE(evalShape)){
    if(is.null(retVals))
      if(!is.null(returns))
        retVals <- seq(mean(returns)-k*sd(returns), mean(returns)+k*sd(returns), by=0.001) else
          retVals <- seq(0.005-k*0.07, 0.005+k*0.07, by=0.001)
        IFvals <- cbind(retVals, IF.fn(retVals, risk="LPM", returns, parsLPM.IF, const, order))
        colnames(IFvals) <- c("r", "IFvals")
        if(isTRUE(IFplot)){
          plot(IFvals[,1], IFvals[,2], type="l", 
               xlab="r", ylab="IF", col="blue", lwd=1, 
               main="LPM",
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
  
  if(order==1){
    
    # Computing the IF vector for LPM (order=1 case)
    IF.LPM.vector <- (const - returns)*(returns <= const) - LPM(returns, const=const, order=order)
    
    # Adding the pre-whitening functionality  
    if(prewhiten)
      IF.LPM.vector <- as.numeric(arima(x=IF.LPM.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)
    
    # Adjustment for data (xts)
    if(xts::is.xts(returns))
      IF.LPM.vector <- xts::xts(IF.LPM.vector, returns.dates)
    
    # Plot of the IF TS
    if(isTRUE(IFplot)){
      print(plot(IF.LPM.vector, type="l", main="LPM Estimator Influence Function Transformed Returns"))
    }
    
    # Returning the IF vector for LPM (order=1 case)
    if(xts::is.xts(returns))
      return(xts::xts(IF.LPM.vector, returns.dates)) else
        return(IF.LPM.vector)
    
  } else if (order==2){
    
    # Computing the IF vector for LPM (order=2 case)
    modified.returns <- (returns - const)^2 * (returns <= const)
    LPM.stored <- LPM(returns, const=const, order=order)
    modified.returns <- modified.returns - LPM.stored^2
    modified.returns <- modified.returns / 2 / LPM.stored
    
    # Adding the pre-whitening functionality  
    if(prewhiten)
      modified.returns <- as.numeric(arima(x=modified.returns, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)
    
    # Adjustment for data (xts)
    if(xts::is.xts(returns))
      modified.returns <- xts::xts(modified.returns, returns.dates)
    
    # Plot of the IF TS
    if(isTRUE(IFplot)){
      print(plot(modified.returns, type="l", main="LPM Estimator Influence Function Transformed Returns", ylab="IF"))
    }
    
    # Stop if no printing of the TS
    if(!IFprint){
      opt <- options(show.error.messages=FALSE)
      on.exit(options(opt)) 
      stop() 
    }
    
    # Returning the IF vector for LPM (k=2 case)
    if(xts::is.xts(returns))
      return(xts::xts(modified.returns, returns.dates)) else
        return(modified.returns)
    
  } else{
    
    # Stop the computation (invalid k parameter)
    stop("Influence Function of LPM is only available for order = 1 or 2", call. = FALSE)
  }
}
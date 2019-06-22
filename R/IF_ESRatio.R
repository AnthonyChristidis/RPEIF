#' @title Influence Function - Expected Shortfall (ES) Ratio
#' 
#' @description \code{IF.ESratio} returns the data and plots the shape of either the IF or the IF TS for the Expected Shortfall Ratio.
#'
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsESratio.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param alpha Tail Probability.
#' @param rf Risk-free interest rate.
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function of ESratio
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @export
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.ESratio(returns=NULL, evalShape=TRUE, 
#'                     retVals=NULL, parsESratio.IF=list(mu=0.01, q.alpha=-0.0541, 
#'                     ES.alpha=0.0777, ES.ratio=0.129),
#'                     IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.ESratio(risk="mean",
#'                     returns=edhec[,"CA"], evalShape=TRUE, 
#'                     retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                     IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.ESratio(returns=edhec[,"CA"], evalShape=FALSE, 
#'                     retVals=NULL, nuisance.par=NULL,
#'                     IFplot=TRUE, IFprint=TRUE,
#'                     prewhiten=FALSE,
#'                     cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF.ESratio <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsESratio.IF=list(mu=0.01, q.alpha=-0.0541, ES.alpha=0.0777, ES.ratio=0.129), k=4,
                       IFplot=FALSE, IFprint=TRUE,
                       alpha=0.1, rf=0, prewhiten=FALSE, ar.prewhiten.order=1,
                       cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                       ...){
  
  # Checking data if IF risk evaluation
  if(isTRUE(evalShape))
    if(is.null(returns))
      if(is.null(parsESratio.IF))
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
  
  # Checking data for alpha value
  if(!inherits(alpha, "numeric")){
    stop("alpha should be numeric")
  } else if(any(alpha < 0, alpha > 1)) {
    stop("alpha should be a numeric value between 0 and 1.")
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
  
  # Function evaluation
  if(isTRUE(evalShape)){
    if(is.null(retVals))
      if(!is.null(returns))
        retVals <- seq(mean(returns)-k*sd(returns), mean(returns)+k*sd(returns), by=0.001) else
          retVals <- seq(0.005-k*0.07, 0.005+k*0.07, by=0.001)
        IFvals <- cbind(retVals, IF.fn(retVals, risk="ESratio", returns, parsESratio.IF, alpha, rf))
        colnames(IFvals) <- c("r", "IFvals")
        if(isTRUE(IFplot)){
          plot(IFvals[,1], IFvals[,2], type="l", 
               xlab="r", ylab="IF", col="blue", lwd=1, 
               main="ES Ratio ",
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
  
  # Computing the mean of the returns
  mu.hat <- mean(returns)
  # Computing the SD of the returns
  sigma.hat <- mean((returns-mu.hat)^2)
  # Computing the VaR of the returns
  VaR.hat <- -quantile(returns, alpha)
  # Storing the negative value of the VaR based on the desired alpha
  quantile.alpha <- -VaR.hat
  # Computing the ES of the returns
  ES.hat <- -mean(returns[returns<=-VaR.hat])
  # Computing the ESratio of the returns
  ESratio.hat <- (mu.hat-rf)/ES.hat
  
  # Computing the final ESratio estimate
  IF.ESratio.vector <- (returns-mu.hat-rf)/ES.hat-ESratio.hat/ES.hat*(1/alpha*((-returns-VaR.hat)*(returns<=quantile.alpha))+VaR.hat-ES.hat)
  
  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.ESratio.vector <- as.numeric(arima(x=IF.ESratio.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)
  
  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.ESratio.vector <- xts::xts(IF.ESratio.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.ESratio.vector, type="l", main="ES Ratio Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the final ESratio estimate
  if(xts::is.xts(returns))
    return(xts::xts(IF.ESratio.vector, returns.dates)) else
      return(IF.ESratio.vector)
}
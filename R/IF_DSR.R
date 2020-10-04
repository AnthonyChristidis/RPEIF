#' @title Influence Function - Downside Sharpe Ratio
#' 
#' @description \code{IF.DSR} returns the data and plots the shape of either the IF or the IF TS for the Downside Sharpe Ratio.
#'
#' @param returns Returns data of the asset or portfolio. This can be a numeric or an xts object.
#' @param evalShape Evaluation of the shape of the IF risk or performance measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param nuisPars Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param rf Risk-free interest rate.
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether outliers are cleaned with a robust location and scale estimator.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "locScaleRob" (default) and "Boudt" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm.
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Addtional parameters.
#'
#' @return Influence function of DSR.
#' 
#' @details 
#' For further details on the usage of the \code{nuisPars} argument, please refer to Section 3.1 for the \code{RPEIF} vignette.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @export
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.DSR(returns=NULL, evalShape=TRUE, 
#'                 retVals=NULL, nuisPars=NULL,
#'                 IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.DSR(returns=edhec[,"CA"], evalShape=TRUE, 
#'                 retVals=seq(-0.1, 0.1, by=0.001), nuisPars=NULL,
#'                 IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with prewhitening) with a plot of IF TS
#' outIF <- IF.DSR(returns=edhec[,"CA"], evalShape=FALSE, 
#'                 retVals=NULL, nuisPars=NULL,
#'                 IFplot=TRUE, IFprint=TRUE,
#'                 prewhiten=FALSE)
#'
IF.DSR <- function(returns=NULL, evalShape=FALSE, retVals=NULL, nuisPars=NULL, k=4,
                   IFplot=FALSE, IFprint=TRUE,
                   rf=0, prewhiten=FALSE, ar.prewhiten.order=1,
                   cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                   ...){
  
  # Evaluation of nuisance parameters
  if(is.null(nuisPars))
    nuisPars <- nuisParsFn() else{
      if(!is.null(nuisPars$mu)){
        nuis.mu <- nuisPars$mu} else{
          nuis.mu <- 0.01}
      if(!is.null(nuisPars$sd)){
        nuis.sd <- nuisPars$sd} else{
          nuis.sd <- 0.05}
      if(!is.null(nuisPars$c)){
        nuis.c <- nuisPars$c} else{
          nuis.c <- 0}
      if(!is.null(nuisPars$alpha)){
        nuis.alpha <- nuisPars$alpha} else{
          nuis.alpha <- 0.1}
      if(!is.null(nuisPars$beta)){
        nuis.beta <- nuisPars$beta} else{
          nuis.beta <- 0.1}
      nuisPars <- nuisParsFn(nuis.mu, nuis.sd, nuis.c, nuis.alpha, nuis.beta)
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
  
  # Function evaluation
  if(isTRUE(evalShape)){
    if(is.null(retVals))
      if(!is.null(returns))
        retVals <- seq(mean(returns)-k*sd(returns), mean(returns)+k*sd(returns), by=0.001) else
          retVals <- seq(0.005-k*0.07, 0.005+k*0.07, by=0.001)
        IFvals <- cbind(retVals, IF.fn(retVals, estimator="DSR", returns=returns, nuisPars , rf=rf, threshold="mean"))
        colnames(IFvals) <- c("r", "IFvals")
        if(isTRUE(IFplot)){
          plot(IFvals[,1], IFvals[,2], type="l", 
               xlab="r", ylab="IF", col="blue", lwd=1, 
               main="DSR",
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
  
  # Computing the IF vector for DSR
  IF.DSR.vector <- -SoR.hat/2/sigma.minus.hat^2*(returns-mu.hat-rf)^2*(returns<=mu.hat)+
    (1/sigma.minus.hat+SoR.hat*mu1.minus.hat/sigma.minus.hat^2)*(returns-mu.hat-rf)+
    SoR.hat/2
  
  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.DSR.vector <- as.numeric(arima(x=IF.DSR.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)
  
  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.DSR.vector <- xts::xts(IF.DSR.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.DSR.vector, type="l", main="DSR Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the IF vector for DSR
  if(xts::is.xts(returns))
    return(xts::xts(IF.DSR.vector, returns.dates)) else
      return(IF.DSR.vector)
}
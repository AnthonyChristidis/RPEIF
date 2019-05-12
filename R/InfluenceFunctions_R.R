#' @useDynLib IFs
#' @importFrom Rcpp sourceCpp
NULL

#' @title Influence Function for Available Risk Measures
#'
#' @description \code{IF} returns the data and plots the shape of either the IF or the IF TS for a risk measure specified.
#' 
#' @param risk Risk measure.
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param nuisance.par Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param compile Boolean variable to indicate if the IF TS should be computed using compiled code (C++) (TRUE) or not (FALSE).
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the IF TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters passed on to influence function of risk measure.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @examples
#' # Plot of IF using the wrapper function
#' outIF <- IF(risk="mean",
#'             returns=NULL, evalShape=TRUE, retVals=NULL, nuisance.par=list(mu=0.005),
#'             IFplot=TRUE, IFprint=TRUE)
#' 
#' #' # Loading data (hedge funds returns)
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#'                     
#' # Plot of IF using wrapper function and with a specified TS 
#' outIF <- IF(risk="mean",
#'             returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'             IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF(risk="mean",
#'             returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'             IFplot=TRUE, IFprint=TRUE,
#'             compile=TRUE, prewhiten=FALSE,
#'             cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF <- function(risk,
               returns=NULL, evalShape=FALSE, retVals=NULL, nuisance.par=NULL, k=4,
               IFplot=FALSE, IFprint=TRUE,
               compile=TRUE, prewhiten=FALSE, ar.prewhiten.order=1,
               cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
               ...){
  
  # Available risk measures
  risk.available <- c("mean", "SD", "VaR", "ES", "SR", "SoR", "ESratio", "VaRratio", "RachR", "LPM", "Omega", "SSD")
  
  # Checking if the specified risk measure is available
  if(!(risk %in% risk.available))
    stop("The specified risk measure is not available.")
  
  # Computation for the specified risk measure
  switch(risk,
         mean = IF.mean(returns=returns, evalShape=evalShape, retVals=retVals, parsMean.IF=nuisance.par, k=k,
                        IFplot=IFplot, IFprint=IFprint,
                        compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                        cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                        ...),
         SD = IF.SD(returns=returns, evalShape=evalShape, retVals=retVals, parsSD.IF=nuisance.par, k=k,
                    IFplot=IFplot, IFprint=IFprint,
                    compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                    cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                    ...),
         VaR = IF.VaR(returns=returns, evalShape=evalShape, retVals=retVals, parsVaR.IF=nuisance.par, k=k,
                      IFplot=IFplot, IFprint=IFprint,
                      compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                      ...),
         ES = IF.ES(returns=returns, evalShape=evalShape, retVals=retVals, parsES.IF=nuisance.par, k=k,
                    IFplot=IFplot, IFprint=IFprint,
                    compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                    cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                    ...),
         SR = IF.SR(returns=returns, evalShape=evalShape, retVals=retVals, parsSR.IF=nuisance.par, k=k,
                    IFplot=IFplot, IFprint=IFprint,
                    compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                    cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                    ...),
         SoR = IF.SoR(returns=returns, evalShape=evalShape, retVals=retVals, parsSoR.IF=nuisance.par, k=k,
                      IFplot=IFplot, IFprint=IFprint,
                      compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                      ...),
         ESratio = IF.ESratio(returns=returns, evalShape=evalShape, retVals=retVals, parsESratio.IF=nuisance.par, k=k,
                              IFplot=IFplot, IFprint=IFprint,
                              compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                              cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                              ...),
         VaRratio = IF.VaRratio(returns=returns, evalShape=evalShape, retVals=retVals, parsVaRratio.IF=nuisance.par, k=k,
                                IFplot=IFplot, IFprint=IFprint,
                                compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                                cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                              ...),
         RachR = IF.RachR(returns=returns, evalShape=evalShape, retVals=retVals, parsRachR.IF=nuisance.par, k=k,
                          IFplot=IFplot, IFprint=IFprint,
                          compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                          cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                            ...),
         LPM = IF.LPM(returns=returns, evalShape=evalShape, retVals=retVals, parsLPM.IF=nuisance.par, k=k,
                      IFplot=IFplot, IFprint=IFprint,
                      compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                      ...),
         Omega = IF.Omega(returns=returns, evalShape=evalShape, retVals=retVals, parsOmega.IF=nuisance.par, k=k,
                          IFplot=IFplot, IFprint=IFprint,
                          compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                          cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                          ...),
         SSD = IF.SSD(returns=returns, evalShape=evalShape, retVals=retVals, parsSSD.IF=nuisance.par, k=k,
                      IFplot=IFplot, IFprint=IFprint,
                      compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                      ...))
}

#' @title Influence Function - Mean
#'
#' @description \code{IF.mean} returns the data and plots the shape of either the IF or the IF TS for the mean.
#'
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsMean.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param compile Boolean variable to indicate if the IF TS should be computed using compiled code (C++) (TRUE) or not (FALSE).
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function for the specified risk measure.
#'
#' @export
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.mean(returns=NULL, evalShape=TRUE, retVals=NULL, parsMean.IF=list(mu=0.005),
#'                  IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.mean(risk="mean",
#'                  returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                  IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.mean(returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'                  IFplot=TRUE, IFprint=TRUE,
#'                  compile=TRUE, prewhiten=FALSE,
#'                  cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF.mean <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsMean.IF=list(mu=0.005), k=4,
                    IFplot=FALSE, IFprint=TRUE, 
                    compile=TRUE, prewhiten=FALSE, ar.prewhiten.order=1,
                    cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                    ...){
  
  # Checking data if IF risk evaluation
  if(isTRUE(evalShape))
    if(is.null(returns))
      if(is.null(parsMean.IF))
        stop("For shape evaluation, nuisance parameters must be specified if no returns are provided.")
  
  # Checking data if IF TS evaluation
  if(!isTRUE(evalShape))
    if(is.null(returns))
      stop("Returns must be provided for the IF TS evaluation.")
  
  # Checking the data for the returns
  if(!is.null(returns)){
    if(!any(c(inherits(returns, "matrix"), inherits(returns, "numeric"), inherits(returns, "xts"), inherits(returns, "zoo")))){
      stop("returns should belong to one of the following classes: matrix, numeric, xts, zoo, xts, zoo")
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
          retVals <- seq(parsMean.IF$mu-k*0.07, parsMean.IF$mu+k*0.07, by=0.001)
    IFvals <- cbind(retVals, IF.fn(retVals, risk="mean", returns, parsMean.IF))
    colnames(IFvals) <- c("r", "IFvals")
    if(isTRUE(IFplot)){
      plot(IFvals[,1], IFvals[,2], type="l", 
           xlab="r", ylab="IF", col="blue", lwd=1, 
           main="Influence Function - Mean",
           panel.first=grid())
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
    IF.mean.vector <- as.vector(IF_mean(returns))
  } else{   
    # Computing the mean of the returns
    mu.hat <- mean(returns)
    
    # Computing the IF vector for the mean of the returns
    IF.mean.vector <- returns - mu.hat
  }
  
  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.mean.vector <- as.numeric(arima(x=IF.mean.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)
  
  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.mean.vector <- xts::xts(IF.mean.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.mean.vector, type="l", main="Mean Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }

  # Returning the IF vector for the mean of the returns
  return(IF.mean.vector)
}

#' @title Influence Function - Standard Deviation 
#' 
#' @description \code{IF.SD} returns the data and plots the shape of either the IF or the IF TS for the standard deviation
#'
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsSD.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param compile Boolean variable to indicate if the IF TS should be computed using compiled code (C++) (TRUE) or not (FALSE).
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function of the standard deviation.
#'
#' @export
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.SD(returns=NULL, evalShape=TRUE, retVals=NULL, parsSD.IF=list(mu=0.005, sd=0.07),
#'                IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.SD(risk="mean",
#'                returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.SD(returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'                IFplot=TRUE, IFprint=TRUE,
#'                compile=TRUE, prewhiten=FALSE,
#'                cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
#'
IF.SD <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsSD.IF=list(mu=0.005, sd=0.07), k=4,
                  IFplot=FALSE, IFprint=TRUE,
                  compile=TRUE, prewhiten=FALSE, ar.prewhiten.order=1,
                  cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                  ...){
  
  # Checking data if IF risk evaluation
  if(isTRUE(evalShape))
    if(is.null(returns))
      if(is.null(parsSD.IF))
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
          retVals <- seq(parsSD.IF$mu-k*parsSD.IF$sd, parsSD.IF$mu+k*parsSD.IF$sd, by=0.001)
    IFvals <- cbind(retVals, IF.fn(retVals, risk="SD", returns, parsSD.IF))
    colnames(IFvals) <- c("r", "IFvals")
    if(isTRUE(IFplot)){
      plot(IFvals[,1], IFvals[,2], type="l", 
           xlab="r", ylab="IF", col="blue", lwd=1, 
           main="Influence Function - SD",
           panel.first=grid())
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
    IF.SD.vector <- as.vector(IF_SD(returns))
  } else{
    # Computing hte mean of the returns
    mu.hat <- mean(returns)
    # Computing the standard deviation of the returns
    sd.hat <- sd(returns)
    
    # Computing the IF vector for the standard deviation
    IF.SD.vector <- ((returns-mu.hat)^2-sd.hat^2)/2/sd.hat
  }
  
  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.SD.vector <- as.numeric(arima(x=IF.SD.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)

  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.SD.vector <- xts::xts(IF.SD.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.SD.vector, type="l", main="SD Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the IF vector for the standard deviation
  if(xts::is.xts(returns))
    return(xts::xts(IF.SD.vector, returns.dates)) else
      return(IF.SD.vector)
}

#' @title Influence Function - Value at Risk (VaR)
#' 
#' @description \code{IF.VaR} returns the data and plots the shape of either the IF or the IF TS for the Value at Risk
#'
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsVaR.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param alpha The tail probability of interest.
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function of the VaR.
#'
#' @export
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.VaR(returns=NULL, evalShape=TRUE, retVals=NULL, parsVaR.IF=list(q.alpha=-0.0847, fq.alpha=2.507),
#'                 IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.VaR(risk="mean",
#'                 returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                 IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.VaR(returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'                 IFplot=TRUE, IFprint=TRUE,
#'                 compile=TRUE, prewhiten=FALSE,
#'                 cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF.VaR <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsVaR.IF=list(q.alpha=-0.0847, fq.alpha=2.507), k=4,
                   IFplot=FALSE, IFprint=TRUE,
                   alpha=0.05, prewhiten=FALSE, ar.prewhiten.order=1,
                   cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                   ...){

  # Checking data if IF risk evaluation
  if(isTRUE(evalShape))
    if(is.null(returns))
      if(is.null(parsVaR.IF))
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
  
  # Checking data for alpha
  if(!inherits(alpha, "numeric")){
    stop("alpha should be numeric")
  } else if(any(alpha < 0, alpha > 1)) {
    stop("alpha should be a numeric value between 0 and 1.")
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
    IFvals <- cbind(retVals, IF.fn(retVals, risk="VaR", returns, parsVaR.IF, alpha))
    colnames(IFvals) <- c("r", "IFvals")
    if(isTRUE(IFplot)){
      plot(IFvals[,1], IFvals[,2], type="l", 
           xlab="r", ylab="IF", col="blue", lwd=1, 
           main="Influence Function - VaR",
           panel.first=grid())
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
  
  # Fitting a density function to the returns
  density.fit <- approxfun(density(returns))
  # Finding the quantile of the density fit based on the desired tail probability
  quantile.alpha <- quantile(returns, alpha)

  # Computing the IF vector for the VaR
  IF.VaR.vector <- ((returns<=quantile.alpha)-alpha)/density.fit(quantile.alpha)

  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.VaR.vector <- as.numeric(arima(x=IF.VaR.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)

  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.VaR.vector <- xts::xts(IF.VaR.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.VaR.vector, type="l", main="VaR Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the IF vector for the VaR
  if(xts::is.xts(returns))
    return(xts::xts(IF.VaR.vector, returns.dates)) else
      return(IF.VaR.vector)
}

#' @title Influence Function - Expected Shortfall (ES)
#' 
#' @description \code{IF.ES} returns the data and plots the shape of either the IF or the IF TS for the ES
#'
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsES.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param alpha.ES Tail Probability.
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function of the ES.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @export
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.ES(returns=NULL, evalShape=TRUE, retVals=NULL, parsES.IF=list(q.alpha=-0.0847, es.alpha=0.273),
#'                IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.ES(risk="mean",
#'                returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.ES(returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'                IFplot=TRUE, IFprint=TRUE,
#'                compile=TRUE, prewhiten=FALSE,
#'                cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF.ES <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsES.IF=list(q.alpha=-0.0847, es.alpha=0.273), k=4,
                  IFplot=FALSE, IFprint=TRUE,
                  alpha.ES=0.05, prewhiten=FALSE, ar.prewhiten.order=1,
                  cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                  ...){

  # Checking data if IF risk evaluation
  if(isTRUE(evalShape))
    if(is.null(returns))
      if(is.null(parsES.IF))
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
  
  # Checking data for alpha.ES value
  if(!inherits(alpha.ES, "numeric")){
    stop("alpha.ES should be numeric")
  } else if(any(alpha.ES < 0, alpha.ES > 1)) {
    stop("alpha.ES should be a numeric value between 0 and 1.")
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
    IFvals <- cbind(retVals, IF.fn(retVals, risk="ES", returns, parsES.IF, alpha.ES))
    colnames(IFvals) <- c("r", "IFvals")
    if(isTRUE(IFplot)){
      plot(IFvals[,1], IFvals[,2], type="l", 
           xlab="r", ylab="IF", col="blue", lwd=1, 
           main="Influence Function - Expected Shortfall",
           panel.first=grid())
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
  
  # Fitting a density function to the returns
  pdf.fit <- approxfun(density(returns))
  # Finding the quantile of the density fit based on the desired tail probability
  quantile.alpha <- quantile(returns, alpha.ES)
  # Computing the ES parameter based on the desired tail probability
  ES.tail <- -mean(returns[returns<=quantile.alpha])

  # Computing the IF vector for the ES
  IF.ES.vector <- 1/alpha.ES*((-returns+quantile.alpha)*(returns<=quantile.alpha))-quantile.alpha-ES.tail
  
  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.ES.vector <- as.numeric(arima(x=IF.ES.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)

  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.ES.vector <- xts::xts(IF.ES.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.ES.vector, type="l", main="ES Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the IF vector for the ES
  if(xts::is.xts(returns))
    return(xts::xts(IF.ES.vector, returns.dates)) else
      return(IF.ES.vector)
}

#' @title Influence Function - Sharpe Ratio (SR)
#' 
#' @description \code{IF.SR} returns the data and plots the shape of either the IF or the IF TS for the SR
#'
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsSR.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param rf Risk-free interest rate.
#' @param compile Boolean variable to indicate if the IF TS should be computed using compiled code (C++) (TRUE) or not (FALSE).
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function of the SR.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @export
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.SR(returns=NULL, evalShape=TRUE, retVals=NULL, parsSR.IF=list(mu.e=0.01, sd=0.05, sr=0.69),
#'                IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.SR(risk="mean",
#'                returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.SR(returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'                IFplot=TRUE, IFprint=TRUE,
#'                compile=TRUE, prewhiten=FALSE,
#'                cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF.SR <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsSR.IF=list(mu.e=0.01, sd=0.05, sr=0.69), k=4,
                  IFplot=FALSE, IFprint=TRUE,
                  rf=0, compile=TRUE, prewhiten=FALSE, ar.prewhiten.order=1,
                  cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                  ...){
  
  # Checking data if IF risk evaluation
  if(isTRUE(evalShape))
    if(is.null(returns))
      if(is.null(parsSR.IF))
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
          retVals <- seq(parsSR.IF$mu.e-k*parsSR.IF$sd, parsSR.IF$mu.e+k*parsSR.IF$sd, by=0.001)
    IFvals <- cbind(retVals, IF.fn(retVals, risk="SR", returns, parsSR.IF, rf))
    colnames(IFvals) <- c("r", "IFvals")
    if(isTRUE(IFplot)){
      plot(IFvals[,1], IFvals[,2], type="l", 
           xlab="r", ylab="IF", col="blue", lwd=1, 
           main="Influence Function - Sharpe Ratio",
           panel.first=grid())
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
    IF.SR.vector <- as.vector(IF_SR(returns, rf))
  } else{
    # Computing the mean of the returns
    mu.hat <- mean(returns)
    # Computing the SD of the returns
    sd.hat <- sd(returns)
    
    # Computing the IF vector for the SR
    IF.SR.vector <- 1/sd.hat*(returns-mu.hat)-1/2*mu.hat/sd.hat^3*((returns-mu.hat)^2-sd.hat^2)
  }

  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.SR.vector <- as.numeric(arima(x=IF.SR.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)

  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.SR.vector <- xts::xts(IF.SR.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.SR.vector, type="l", main="SR Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the IF vector for the SR
  if(xts::is.xts(returns))
    return(xts::xts(IF.SR.vector, returns.dates)) else
      return(IF.SR.vector)
}

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
           main="Influence Function - Sortino Ratio",
           panel.first=grid())
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
           main="Influence Function - Sortino Ratio",
           panel.first=grid())
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
#' outIF <- IF.SoR(returns=NULL, evalShape=TRUE, retVals=NULL, parsSoR.IF=list(mu=0.01, lpm2=0.00898, sor.c=0.3337, 
#'                                                                             ssd=0.0354, smean=-0.0199, sor.mu=0.2929),
#'                  IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.SoR(risk="mean",
#'                  returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                  IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.SoR(returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'                  IFplot=TRUE, IFprint=TRUE,
#'                  compile=TRUE, prewhiten=FALSE,
#'                  cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF.SoR <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsSoR.IF=list(mu=0.01, lpm2=0.00898, sor.c=0.3337, 
                                                                                ssd=0.0354, smean=-0.0199, sor.mu=0.2929), k=4,
                   IFplot=FALSE, IFprint=TRUE,
                   threshold=c("mean", "const")[1], const=0, rf=0, prewhiten=FALSE, ar.prewhiten.order=1,
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
#' outIF <- IF.ESratio(returns=NULL, evalShape=TRUE, retVals=NULL, parsESratio.IF=list(mu=0.01, q.alpha=-0.0541, ES.alpha=0.0777, ES.ratio=0.129),
#'                     IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.ESratio(risk="mean",
#'                     returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                     IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.ESratio(returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'                     IFplot=TRUE, IFprint=TRUE,
#'                     compile=TRUE, prewhiten=FALSE,
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
           main="Influence Function - ES ratio ",
           panel.first=grid())
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

#' @title Influence Function - Value at Risk (VaR) Ratio 
#' 
#' @description \code{IF.VaRratio} returns the data and plots the shape of either the IF or the IF TS for the VaR Ratio.
#'
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsVaRratio.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param alpha The tail probability of interest.
#' @param rf Risk-free interest rate.
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function of the VaRratio.
#'
#' @export
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.VaRratio(returns=NULL, evalShape=TRUE, retVals=NULL, parsVaRratio.IF=list(mu=0.1, q.alpha=-0.0541, fq.alpha=3.99, VaR.ratio=0.185),
#'                      IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.VaRratio(risk="mean",
#'                      returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                      IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.VaRratio(returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'                      IFplot=TRUE, IFprint=TRUE,
#'                      compile=TRUE, prewhiten=FALSE,
#'                      cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF.VaRratio <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsVaRratio.IF=list(mu=0.1, q.alpha=-0.0541, fq.alpha=3.99, VaR.ratio=0.185), k=4,
                        IFplot=FALSE, IFprint=TRUE,
                        alpha=0.05, rf=0, prewhiten=FALSE, ar.prewhiten.order=1,
                        cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                        ...){
  
  # Checking data if IF risk evaluation
  if(isTRUE(evalShape))
    if(is.null(returns))
      if(is.null(parsVaRratio.IF))
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
    IFvals <- cbind(retVals, IF.fn(retVals, risk="VaRratio", returns, parsVaRratio.IF, alpha, rf))
    colnames(IFvals) <- c("r", "IFvals")
    if(isTRUE(IFplot)){
      plot(IFvals[,1], IFvals[,2], type="l", 
           xlab="r", ylab="IF", col="blue", lwd=1, 
           main="Influence Function - VaR Ratio",
           panel.first=grid())
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
  
  # Mean of returns
  mu.hat <- mean(returns)
  # Fitting a density function to the returns
  density.fit <- approxfun(density(returns))
  # Finding the quantile of the density fit based on the desired tail probability
  quantile.alpha <- quantile(returns, alpha)
  # Computing the VaR ratio
  VaRratio.hat <- (mu.hat - rf)/quantile.alpha
  
  # Computing the IF vector for the VaR
  IF.VaRratio.vector <- -(returns-mu.hat)/(-quantile.alpha) - (VaRratio.hat/quantile.alpha)*((returns<=quantile.alpha)-alpha)/density.fit(quantile.alpha)
  
  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.VaRratio.vector <- as.numeric(arima(x=IF.VaRratio.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)

  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.VaRratio.vector <- xts::xts(IF.VaRratio.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.VaRratio.vector, type="l", main="VaR Ratio Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the IF vector for the VaR
  if(xts::is.xts(returns))
    return(xts::xts(IF.VaRratio.vector, returns.dates)) else
      return(IF.VaRratio.vector)
}

#' @title Influence Function - Rachev Ratio
#' 
#' @description \code{IF.Rachev} returns the data and plots the shape of either the IF or the IF TS for the Rachev Ratio.
#'
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsRachev.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param alpha Lower tail probability.
#' @param beta Upper tail probability.
#' @param rf Risk-free interest rate.
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function of Rachev Ratio.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @export
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.RachR(returns=NULL, evalShape=TRUE, retVals=NULL, parsRachev.IF=list(q.alpha=-0.0541, es.alpha=0.0777, q.beta=0.0741, eg.beta=0.0977, rach.r=1.257),
#'                   IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.RachR(risk="mean",
#'                   returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                   IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.RachR(returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'                   IFplot=TRUE, IFprint=TRUE,
#'                   compile=TRUE, prewhiten=FALSE,
#'                   cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF.RachR <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsRachev.IF=list(q.alpha=-0.0541, es.alpha=0.0777, q.beta=0.0741, eg.beta=0.0977, rach.r=1.257), k=4,
                      IFplot=FALSE, IFprint=TRUE,
                      alpha=0.1, beta=0.1, rf=0, prewhiten=FALSE, ar.prewhiten.order=1,
                      cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                      ...){
  
  # Checking data if IF risk evaluation
  if(isTRUE(evalShape))
    if(is.null(returns))
      if(is.null(parsRachev.IF))
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
  
  # Checking data for beta value
  if(!inherits(beta, "numeric")){
    stop("beta should be numeric")
  } else if(any(beta < 0, beta > 1)) {
    stop("beta should be a numeric value between 0 and 1.")
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
    IFvals <- cbind(retVals, IF.fn(retVals, risk="Rachev", returns, parsRachev.IF, alpha, beta, rf))
    colnames(IFvals) <- c("r", "IFvals")
    if(isTRUE(IFplot)){
      plot(IFvals[,1], IFvals[,2], type="l", 
           xlab="r", ylab="IF", col="blue", lwd=1, 
           main="Influence Function - Rachev Ratio",
           panel.first=grid())
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
  
  # Computing the VaR of the returns (lower tail)
  VaR.hat.lower <- -quantile(returns, alpha)
  # Storing the negative value of the VaR based on the desired alpha (lower tail)
  quantile.lower <- -VaR.hat.lower
  # Computing the ES of the returns (lower tail)
  ES.lower <- -mean(returns[returns<=-VaR.hat.lower])
  
  # Computing the VaR of the returns (upper tail)
  n.upper <- floor((1-beta)*length(returns))
  sorted.returns <- sort(as.numeric(returns))
  VaR.hat.upper <- sorted.returns[n.upper]
  # Storing the negative value of the VaR based on the desired alpha (upper tail)
  quantile.upper <- VaR.hat.upper
  # Computing the ES of the returns (upper tail)
  ES.upper <- mean(returns[returns>=quantile.upper])
  
  # Computing the final ESratio estimate
  IF.Rachev.vector <- (1/ES.lower)*((1/beta)*(returns-rf-quantile.upper)*(returns>=quantile.upper)+quantile.upper-ES.upper) - 
    (ES.upper/ES.lower^2)*(1/alpha * (-returns + rf - quantile.lower)*(returns<=quantile.lower)+quantile.lower-ES.lower)

  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.Rachev.vector <- as.numeric(arima(x=IF.Rachev.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)
  
  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.Rachev.vector <- xts::xts(IF.Rachev.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.Rachev.vector, type="l", main="Rachev Ratio Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the final ESratio estimate
  if(xts::is.xts(returns))
    return(xts::xts(IF.Rachev.vector, returns.dates)) else
      return(IF.Rachev.vector)
}

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
#' outIF <- IF.LPM(returns=NULL, evalShape=TRUE, retVals=NULL, parsLPM.IF=list(lpm1=0.0255, lpm2=0.00218),
#'                 IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.LPM(risk="mean",
#'                 returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                 IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.LPM(returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'                 IFplot=TRUE, IFprint=TRUE,
#'                 compile=TRUE, prewhiten=FALSE,
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
           main="Influence Function - LPM",
           panel.first=grid())
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

#' @title Influence Function - Omega Ratio
#' 
#' @description \code{IF.OmegaRatio} returns the data and plots the shape of either the IF or the IF TS for the Omega Ratio.
#'
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsOmegaRatio.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param const Constant threshold.
#' @param compile Boolean variable to indicate if the IF TS should be computed using compiled code (C++) (TRUE) or not (FALSE).
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function of Omega Ratio.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @export
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.Omega(returns=NULL, evalShape=TRUE, retVals=NULL, parsOmegaRatio.IF=list(lpm1=0.0153, upm1=0.0253, omega=1.652),
#'                   IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.Omega(risk="mean",
#'                   returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                   IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.Omega(returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'                   IFplot=TRUE, IFprint=TRUE,
#'                   compile=TRUE, prewhiten=FALSE,
#'                   cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF.Omega <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsOmegaRatio.IF=list(lpm1=0.0153, upm1=0.0253, omega=1.652), k=4,
                     IFplot=FALSE, IFprint=TRUE,
                     const=0, compile=TRUE, prewhiten=FALSE, ar.prewhiten.order=1,
                     cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                     ...){
  
  # Checking data if IF risk evaluation
  if(isTRUE(evalShape))
    if(is.null(returns))
      if(is.null(parsOmegaRatio.IF))
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
    IFvals <- cbind(retVals, IF.fn(retVals, risk="OmegaRatio", returns, parsOmegaRatio.IF, const))
    colnames(IFvals) <- c("r", "IFvals")
    if(isTRUE(IFplot)){
      plot(IFvals[,1], IFvals[,2], type="l", 
           xlab="r", ylab="IF", col="blue", lwd=1, 
           main="Influence Function - Omega Ratio",
           panel.first=grid())
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
    IF.Omega.vector <- as.vector(IF_OmegaRatio(returns, const))
  } else{
    # Returning length of returns vector
    N <- length(returns)
    
    # Computing Omega+
    Omega_p <- sum(returns[returns>=const]-const)/N
    # Computing Omega-
    Omega_m <- sum(const-returns[returns<=const])/N
    # Computing the IF vector for Omega Ratio
    IF.Omega.vector <- ((returns - const) * (returns >= const) - Omega_p)/ Omega_m
    IF.Omega.vector <- IF.Omega.vector - Omega_p / Omega_m^2 * ((const - returns) * (returns <= const) - Omega_m) 
  }

  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.Omega.vector <- as.numeric(arima(x=IF.Omega.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)

  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.Omega.vector <- xts::xts(IF.Omega.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.Omega.vector, type="l", main="Omega Ratio Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the IF vector for Omega Ratio
  if(xts::is.xts(returns))
    return(xts::xts(IF.Omega.vector, returns.dates)) else
      return(IF.Omega.vector)
}

#' @title Influence Function - Semi-Standard Deviation (SSD)
#' 
#' @description \code{IF.SSD} returns the data and plots the shape of either the IF or the IF TS for the SSD
#'
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsSSD.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param rf Risk-free interest rate.
#' @param compile Boolean variable to indicate if the IF TS should be computed using compiled code (C++) (TRUE) or not (FALSE).
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through a robust filter.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "Boudt" and "locScaleRob" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function of SSD.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @export
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' outIF <- IF.SSD(returns=NULL, evalShape=TRUE, retVals=NULL, parsSSD.IF=list(mu=0.005, ssd=0.0495, smean=-0.0279),
#'                 IFplot=TRUE, IFprint=TRUE)
#'
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#' 
#' # Plot of IF a specified TS 
#' outIF <- IF.SSD(risk="mean",
#'                 returns=edhec[,"CA"], evalShape=TRUE, retVals=seq(-0.1, 0.1, by=0.001), nuisance.par=NULL,
#'                 IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF.SSD(returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
#'                 IFplot=TRUE, IFprint=TRUE,
#'                 compile=TRUE, prewhiten=FALSE,
#'                 cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF.SSD <- function(returns=NULL, evalShape=FALSE, retVals=NULL, parsSSD.IF=list(mu=0.005, ssd=0.0495, smean=-0.0279), k=4,
                   IFplot=FALSE, IFprint=TRUE,
                   rf=0, compile=TRUE, prewhiten=FALSE, ar.prewhiten.order=1,
                   cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
                   ...){

  # Checking data if IF risk evaluation
  if(isTRUE(evalShape))
    if(is.null(returns))
      if(is.null(parsSSD.IF))
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
          retVals <- seq(parsSSD.IF$mu-k*0.07, parsSSD.IF$mu+k*0.07, by=0.001)
    IFvals <- cbind(retVals, IF.fn(retVals, risk="SSD", returns, parsSSD.IF, rf))
    colnames(IFvals) <- c("r", "IFvals")
    if(isTRUE(IFplot)){
      plot(IFvals[,1], IFvals[,2], type="l", 
           xlab="r", ylab="IF", col="blue", lwd=1, 
           main="Influence Function - SSD",
           panel.first=grid())
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
    IF.SSD.vector <- as.vector(IF_SSD(returns, rf))
  } else {
    # Computing the mean of the returns
    mu.hat <- mean(returns)
    # Computing SD- of the returns
    sigma.minus.hat <- sqrt(mean((returns-mu.hat)^2*(returns<=mu.hat)))
    
    # Computing the IF vector for SSD
    IF.SSD.vector <- (returns - mu.hat)^2 * (returns <= mu.hat)
    IF.SSD.vector <- IF.SSD.vector - 2 * mean((returns-mu.hat) * (returns <= mu.hat)) * (returns - mu.hat)
    IF.SSD.vector <- IF.SSD.vector - sigma.minus.hat^2
    IF.SSD.vector <- IF.SSD.vector / 2 / sigma.minus.hat
  }

  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.SSD.vector <- as.numeric(arima(x=IF.SSD.vector, order=c(ar.prewhiten.order,0,0), include.mean=TRUE)$residuals)
  
  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.SSD.vector <- xts::xts(IF.SSD.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.SSD.vector, type="l", main="SSD Estimator Influence Function Transformed Returns", ylab="IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the IF vector for SSD
  if(xts::is.xts(returns))
    return(xts::xts(IF.SSD.vector, returns.dates)) else
      return(IF.SSD.vector)
}

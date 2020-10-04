#'
#' @importFrom stats approxfun arima density dnorm pnorm qnorm quantile sd
#' @importFrom graphics plot grid abline

#' @title Influence Function for Available Risk and Performance Measures 
#'
#' @description \code{IF} returns the data and plots the shape of either the IF or the IF TS for a specified estimator.
#' 
#' @param estimator The estimator of interest.
#' @param returns Returns data of the asset or portfolio. This can be a numeric or an xts object.
#' @param evalShape Evaluation of the shape of the IF risk or performance measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param nuisPars Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether outliers are cleaned with a robust location and scale estimator.
#' @param cleanMethod Robust method used to clean outliers from the TS. The choices are "locScaleRob" (default) and "Boudt" for the function. 
#' @param alpha.robust Tuning parameter for the quantile of the "Boudt" robust data cleaning algorithm, using the minimum covariance determinant estimator (MCD).
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters passed on to influence function of risk or performance measure.
#' 
#' @details 
#' For further details on the usage of the \code{nuisPars} argument, please refer to Section 3.1 for the \code{RPEIF} vignette.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @examples
#' # Plot of IF using the wrapper function
#' outIF <- IF(estimator="mean",
#'             returns=NULL, evalShape=TRUE, retVals=NULL, nuisPars=list(mu=0.005),
#'             IFplot=TRUE, IFprint=TRUE)
#' 
#' #' # Loading data (hedge funds returns)
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#'                     
#' # Plot of IF using wrapper function and with a specified TS 
#' outIF <- IF(estimator="mean",
#'             returns=edhec[,"CA"], evalShape=TRUE, 
#'             retVals=seq(-0.1, 0.1, by=0.001), nuisPars=NULL,
#'             IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with prewhitening) with a plot of IF TS
#' outIF <- IF(estimator="mean",
#'             returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisPars =NULL,
#'             IFplot=TRUE, IFprint=TRUE,
#'             compile=TRUE, prewhiten=FALSE)
#'
IF <- function(estimator,
               returns=NULL, evalShape=FALSE, retVals=NULL, nuisPars =NULL, k=4,
               IFplot=FALSE, IFprint=TRUE,
               prewhiten=FALSE, ar.prewhiten.order=1,
               cleanOutliers=FALSE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05,
               ...){
  
  # Available estimators
  estimator.available <- c("mean", "SD", "VaR", "ES", "SR", "SoR", "DSR", "ESratio", "VaRratio", "RachevRatio", "LPM", "Omega", "SemiSD")
  
  # Checking if the specified estimator is available
  if(!(estimator %in% estimator.available))
    stop("The specified estimator is not available.")
  
  # Computation for the specified estimator
  switch(estimator,
         mean = IF.mean(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                        IFplot=IFplot, IFprint=IFprint,
                        prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                        cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                        ...),
         SD = IF.SD(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                    IFplot=IFplot, IFprint=IFprint,
                    prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                    cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                    ...),
         VaR = IF.VaR(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                      IFplot=IFplot, IFprint=IFprint,
                      prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                      ...),
         ES = IF.ES(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                    IFplot=IFplot, IFprint=IFprint,
                    prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                    cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                    ...),
         SR = IF.SR(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                    IFplot=IFplot, IFprint=IFprint,
                    prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                    cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                    ...),
         SoR = IF.SoR(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                      IFplot=IFplot, IFprint=IFprint,
                      prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                      ...),
         DSR = IF.DSR(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                      IFplot=IFplot, IFprint=IFprint,
                      prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                      ...),
         ESratio = IF.ESratio(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                              IFplot=IFplot, IFprint=IFprint,
                              prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                              cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                              ...),
         VaRratio = IF.VaRratio(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                                IFplot=IFplot, IFprint=IFprint,
                                prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                                cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                              ...),
         RachevRatio = IF.RachevRatio(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                                      IFplot=IFplot, IFprint=IFprint,
                                      prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                                      ...),
         LPM = IF.LPM(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                      IFplot=IFplot, IFprint=IFprint,
                      prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                      ...),
         Omega = IF.Omega(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                          IFplot=IFplot, IFprint=IFprint,
                          prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                          cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                          ...),
         SemiSD = IF.SemiSD(returns=returns, evalShape=evalShape, retVals=retVals, nuisPars =nuisPars , k=k,
                            IFplot=IFplot, IFprint=IFprint,
                            prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                            cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                            ...))
}


#' @useDynLib IFs, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL
#'
#' @import stats
#' @import graphics

#' @title Influence Function for Available Risk Measures
#'
#' @description \code{IF} returns the data and plots the shape of either the IF or the IF TS for a risk measure specified.
#' 
#' @param risk Risk measure.
#' @param returns Vector of the returns of the asset or portfolio.
#' @param evalShape Evaluation of the shape of the IF risk measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param parsNuisance.IF Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
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
#'             returns=NULL, evalShape=TRUE, retVals=NULL, parsNuisance.IF=list(mu=0.005),
#'             IFplot=TRUE, IFprint=TRUE)
#' 
#' #' # Loading data (hedge funds returns)
#' data(edhec, package="PerformanceAnalytics")
#' colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                     "GM", "LS", "MA", "RV", "SS", "FoF") 
#'                     
#' # Plot of IF using wrapper function and with a specified TS 
#' outIF <- IF(risk="mean",
#'             returns=edhec[,"CA"], evalShape=TRUE, 
#'             retVals=seq(-0.1, 0.1, by=0.001), parsNuisance.IF=NULL,
#'             IFplot=TRUE, IFprint=TRUE)
#' 
#' # Computing the IF of the returns (with outlier cleaning and prewhitening) with a plot of IF TS
#' outIF <- IF(risk="mean",
#'             returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, parsNuisance.IF=NULL,
#'             IFplot=TRUE, IFprint=TRUE,
#'             compile=TRUE, prewhiten=FALSE,
#'             cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99)
#'
IF <- function(risk,
               returns=NULL, evalShape=FALSE, retVals=NULL, parsNuisance.IF=NULL, k=4,
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
         mean = IF.mean(returns=returns, evalShape=evalShape, retVals=retVals, parsNuisance.IF=parsNuisance.IF, k=k,
                        IFplot=IFplot, IFprint=IFprint,
                        compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                        cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                        ...),
         SD = IF.SD(returns=returns, evalShape=evalShape, retVals=retVals, parsNuisance.IF=parsNuisance.IF, k=k,
                    IFplot=IFplot, IFprint=IFprint,
                    compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                    cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                    ...),
         VaR = IF.VaR(returns=returns, evalShape=evalShape, retVals=retVals, parsNuisance.IF=parsNuisance.IF, k=k,
                      IFplot=IFplot, IFprint=IFprint,
                      compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                      ...),
         ES = IF.ES(returns=returns, evalShape=evalShape, retVals=retVals, parsNuisance.IF=parsNuisance.IF, k=k,
                    IFplot=IFplot, IFprint=IFprint,
                    compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                    cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                    ...),
         SR = IF.SR(returns=returns, evalShape=evalShape, retVals=retVals, parsNuisance.IF=parsNuisance.IF, k=k,
                    IFplot=IFplot, IFprint=IFprint,
                    compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                    cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                    ...),
         SoR = IF.SoR(returns=returns, evalShape=evalShape, retVals=retVals, parsNuisance.IF=parsNuisance.IF, k=k,
                      IFplot=IFplot, IFprint=IFprint,
                      compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                      ...),
         ESratio = IF.ESratio(returns=returns, evalShape=evalShape, retVals=retVals, parsNuisance.IF=parsNuisance.IF, k=k,
                              IFplot=IFplot, IFprint=IFprint,
                              compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                              cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                              ...),
         VaRratio = IF.VaRratio(returns=returns, evalShape=evalShape, retVals=retVals, parsNuisance.IF=parsNuisance.IF, k=k,
                                IFplot=IFplot, IFprint=IFprint,
                                compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                                cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                              ...),
         RachR = IF.RachR(returns=returns, evalShape=evalShape, retVals=retVals, parsNuisance.IF=parsNuisance.IF, k=k,
                          IFplot=IFplot, IFprint=IFprint,
                          compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                          cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                            ...),
         LPM = IF.LPM(returns=returns, evalShape=evalShape, retVals=retVals, parsNuisance.IF=parsNuisance.IF, k=k,
                      IFplot=IFplot, IFprint=IFprint,
                      compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                      ...),
         Omega = IF.Omega(returns=returns, evalShape=evalShape, retVals=retVals, parsNuisance.IF=parsNuisance.IF, k=k,
                          IFplot=IFplot, IFprint=IFprint,
                          compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                          cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                          ...),
         SSD = IF.SSD(returns=returns, evalShape=evalShape, retVals=retVals, parsNuisance.IF=parsNuisance.IF, k=k,
                      IFplot=IFplot, IFprint=IFprint,
                      compile=compile, prewhiten=prewhiten, ar.prewhiten.order=ar.prewhiten.order,
                      cleanOutliers=cleanOutliers, cleanMethod=cleanMethod, alpha.robust=alpha.robust, eff=eff,
                      ...))
}


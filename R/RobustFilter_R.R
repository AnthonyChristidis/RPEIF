
# Implementation for Robust filtering for the outliers of the IF time series

# Implementation of the robust filter for the IF TS
robust.cleaning <- function(IF.vector, robust.method=c("locScaleRob", "Boudt")[1], 
                            alpha.robust=0.05, normal.efficiency=0.99){
  
  if(robust.method=="locScaleRob"){ # Data cleaning using the Martin method
    
    if(!requireNamespace("RobStatTM", quietly = TRUE)) {
      stop("Package \"RobStatTM\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
  
    # Length of the IF TS
    n <- length(IF.vector)
    # Casting as a numerci
    IF.numeric <- as.numeric(IF.vector)
    # Location and scale parameters
    mu <- RobStatTM::locScaleM(IF.numeric)$mu
    s <- RobStatTM::locScaleM(IF.numeric)$disper
    
    # Efficiency parameter
    if(normal.efficiency==0.95)
      efficiency.param <- 3 else if(normal.efficiency==0.99)
        efficiency.param <- 3.568 else if(normal.efficiency==0.999)
          efficiency.param <- 4.21 else
            warning("Invalid value for the normal distribution efficiency.")
    
    # Computing the limits
    uplim <- rep(mu + s*efficiency.param, n) 
    dnlim <- rep(mu - s*efficiency.param, n)
    
    # Computing the filtered time series
    xcl <- ifelse(IF.numeric >= uplim, uplim, IF.numeric)
    xcl <- ifelse(IF.numeric <= dnlim, dnlim, xcl)
    
    return(xcl)
    
  } else if(robust.method=="Boudt"){ # Data cleaning using the Boudt method
    
    return(PerformanceAnalytics::Return.clean(IF.vector, alpha=alpha.robust, method = "boudt"))
    
  } else
    warning("Invalid name for robust method.") # Otherwise invalid method name for robust data cleaning
  
}




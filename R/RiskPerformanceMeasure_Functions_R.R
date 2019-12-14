
IF.fn <- function(x, risk, returns, nuisance.par, ...){
  
  # Available risk measures
  risk.available <- c("mean", "SD", "VaR", "ES", "SR", "SoR", "ESratio", "VaRratio", "Rachev", "LPM", "OmegaRatio", "SemiSD")
  
  # Checking if the specified risk measure is available
  if(!(risk %in% risk.available))
    stop("The specified risk measure is not available.")
  
  # Plot for the specified risk measure
  switch(risk,
         mean = IF.mean.fn(x, returns, nuisance.par, ...),
         SD = IF.SD.fn(x, returns, nuisance.par, ...),
         VaR = IF.VaR.fn(x, returns, nuisance.par, ...),
         ES = IF.ES.fn(x, returns, nuisance.par, ...),
         SR = IF.SR.fn(x, returns, nuisance.par, ...),
         SoR = IF.SoR.fn(x, returns, nuisance.par, ...),
         ESratio = IF.ESratio.fn(x, returns, nuisance.par, ...),
         VaRratio = IF.VaRratio.fn(x, returns, nuisance.par, ...),
         Rachev = IF.Rachev.fn(x, returns, nuisance.par, ...),
         LPM = IF.LPM.fn(x, returns, nuisance.par, ...),
         OmegaRatio = IF.OmegaRatio.fn(x, returns, nuisance.par, ...),
         SemiSD = IF.SemiSD.fn(x, returns, nuisance.par, ...))
}


IF.mean.fn <- function(x, returns, parsMean.IF){
  
  # Computing the nuisance parameters
  if(is.null(returns))
    mu.hat <- parsMean.IF$mu else
      mu.hat <- mean(returns)
  
  # Computing the IF for the mean of the returns
  IF.mean <- x - mu.hat
  
  # Return value for the risk measure
  return(IF.mean)
}


IF.SD.fn <- function(x, returns, parSemiSD.IF){
  
  # Computing the nuisance parameters
  if(is.null(returns)){
    mu.hat <- parSemiSD.IF$mu
    sd.hat <- parSemiSD.IF$sd
  } else{
      # Computing hte mean of the returns
      mu.hat <- mean(returns)
      # Computing the standard deviation of the returns
      sd.hat <- sd(returns)
    }

  # Computing the IF for the standard deviation
  IF.SD <- ((x-mu.hat)^2-sd.hat^2)/2/sd.hat
  
  # Return value for the risk measure
  return(IF.SD)
}


IF.VaR.fn <- function(x, returns, parsVaR.IF, alpha = 0.05){
  
  # IF for null returns
  if(is.null(returns)){
    # Fitting a density function to the returns
    fq.alpha <- parsVaR.IF$fq.alpha
    # Finding the quantile of the density fit based on the desired tail probability
    quantile.alpha <- parsVaR.IF$q.alpha
    
    # Computing the IF for the VaR
    IF.VaR <- ((x<=quantile.alpha)-alpha)/fq.alpha
    
    # Return value for the risk measure
    return(as.numeric(IF.VaR))
  } else{
    # Fitting a density function to the returns
    density.fit <- approxfun(density(returns))
    # Finding the quantile of the density fit based on the desired tail probability
    quantile.alpha <- quantile(returns, alpha)
    
    # Computing the IF for the VaR
    IF.VaR <- ((x<=quantile.alpha)-alpha)/density.fit(quantile.alpha)
    
    # Return value for the risk measure
    return(as.numeric(IF.VaR))
  }
}


IF.ES.fn <- function(x, returns, parsES.IF, alpha.ES=0.05){
  
  # IF for null returns
  if(is.null(returns)){
    # Finding the quantile of the density fit based on the desired tail probability
    quantile.alpha <- parsES.IF$q.alpha
    # Computing the ES parameter based on the desired tail probability
    ES.tail <- parsES.IF$es.alpha
    
    # Computing the IF for the ES
    IF.ES <- 1/alpha.ES*((-x+quantile.alpha)*(x<=quantile.alpha))-quantile.alpha-ES.tail
    
    # Return value for the risk measure
    return(as.numeric(IF.ES))
  } else{
    # Finding the quantile of the density fit based on the desired tail probability
    quantile.alpha <- quantile(returns, alpha.ES)
    # Computing the ES parameter based on the desired tail probability
    ES.tail <- -mean(returns[returns<=quantile.alpha])
    
    # Computing the IF for the ES
    IF.ES <- 1/alpha.ES*((-x+quantile.alpha)*(x<=quantile.alpha))-quantile.alpha-ES.tail
    
    # Return value for the risk measure
    return(as.numeric(IF.ES))
  }
}


IF.SR.fn <- function(x, returns, parsSR.IF, rf=0){
  
  # Computing the nuisance parameters
  if(is.null(returns)){
    mu.hat <- parsSR.IF$mu
    sd.hat <- parsSR.IF$sd
    sr <- parsSR.IF$sr
    # Return value is parameters are given
    IF.SR <- -sr/(2*sd.hat^2)*(x-mu.hat)^2+1/sd.hat*(x-mu.hat)+sr/2
    return(IF.SR)
  } else{
    # Computing the mean of the returns
    mu.hat <- mean(returns)
    # Computing the SD of the returns
    sd.hat <- sd(returns)
  }
  
  # Computing the IF vector for the SR
  IF.SR <- 1/sd.hat*(x-mu.hat)-1/2*mu.hat/sd.hat^3*((x-mu.hat)^2-sd.hat^2)
  
  # Return value for the risk measure
  return(as.numeric(IF.SR))
}

IF.SoR.fn <- function(x, returns, parsSoR.IF, rf=0, const=0, threshold=c("mean", "const")[1]){
  
  # Case where we want mean threshold
  if(threshold=="mean")
    return(IF.SoR_M.fn(x=x, returns=returns, parsSoR_M.IF=parsSoR.IF, rf=rf)) else if(threshold=="const")
      return(IF.SoR_C.fn(x=x, returns=returns, parsSoR_C.IF=parsSoR.IF, const=const))
}

IF.SoR_M.fn <- function(x, returns, parsSoR_M.IF, rf=0){
  
  # IF for null returns
  if(is.null(returns)){
    # Parameters for null returns
    mu.hat <- parsSoR_M.IF$mu
    sor.mu <- parsSoR_M.IF$sor.mu
    SemiSD <- parsSoR_M.IF$semisd
    smean <- parsSoR_M.IF$smean
    # IF computation for null returns
    IF.SoR_M <- -sor.mu * (x - mu.hat)^2*(x<=mu.hat)/(2*SemiSD^2) + (x - mu.hat)*(1/SemiSD + sor.mu*smean/SemiSD^2) + sor.mu/2
    # Return value for null returns
    return(IF.SoR_M)
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
    IF.SoR_M <- -SoR.hat/2/sigma.minus.hat^2*(x-mu.hat-rf)^2*(x<=mu.hat)+
      (1/sigma.minus.hat+SoR.hat*mu1.minus.hat/sigma.minus.hat^2)*(x-mu.hat-rf)+
      SoR.hat/2
    
    # Return value for the risk measure
    return(IF.SoR_M)
  }
}

IF.SoR_C.fn <- function(x, returns, parsSoR_C.IF, const=0){
  
  # IF for null returns
  if(is.null(returns)){
    # Parameters for null returns
    mu.hat <- parsSoR_C.IF$mu
    sor.c <- parsSoR_C.IF$sor.c
    lpm2 <- parsSoR_C.IF$lpm2
    # IF value for null returns
    IF.SoR_C <- -sor.c*(x-const)^2*(x<=const)/(2*lpm2) + (x-mu.hat)/sqrt(lpm2) + sor.c/2
    # Return value for null returns
    return(IF.SoR_C)
  } else{
    
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
    IF.SoR_C <- -SoR.hat/2/sigma.minus.hat^2*(x-mar.parameter)^2*(x<=mar.parameter)+
      1/sigma.minus.hat*(x-mu.hat)+
      SoR.hat/2
    
    # Return value for the risk measure
    return(IF.SoR_C)
  }
}

IF.ESratio.fn <- function(x, returns, parsESratio.IF, alpha=0.1, rf=0){
  
  # IF for null returns
  if(is.null(returns)){
    # Parameters for null returns
    mu.hat <- parsESratio.IF$mu
    q.alpha <- parsESratio.IF$q.alpha
    ES.alpha <- parsESratio.IF$es.alpha
    ES.ratio <- parsESratio.IF$ES.ratio
    # IF computation for null returns
    IF.ESratio <- (x - mu.hat)/ES.alpha - ES.ratio/ES.alpha*((-x + q.alpha)*(x<=q.alpha)/alpha - q.alpha - ES.alpha)
    # Return value for null returns
    return(IF.ESratio)
  } else{
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
    IF.ESratio <- (x-mu.hat-rf)/ES.hat-ESratio.hat/ES.hat*(1/alpha*((-x-VaR.hat)*(x<=quantile.alpha))+VaR.hat-ES.hat)
    
    # Return value for the risk measure
    return(as.numeric(IF.ESratio))
  }
}

IF.VaRratio.fn <- function(x, returns, parsVaRratio.IF, alpha=0.1, rf=0){
  
  # IF for null returns
  if(is.null(returns)){
    
    # Parameters for null returns
    mu.hat <- parsVaRratio.IF$mu
    VaR.alpha <- -parsVaRratio.IF$q.alpha
    VaR.ratio <- parsVaRratio.IF$VaR.ratio
    fq.alpha <- parsVaRratio.IF$fq.alpha
    # IF computation for null returns
    IF.VaRratio <- (x - mu.hat)/VaR.alpha - (VaR.ratio/VaR.alpha) * ((x <= -VaR.alpha)-alpha)/fq.alpha
    # Return value for null returns
    return(IF.VaRratio)
    
  } else{
    # Mean of returns
    mu.hat <- mean(returns)
    # Fitting a density function to the returns
    density.fit <- approxfun(density(returns))
    # Finding the quantile of the density fit based on the desired tail probability
    quantile.alpha <- quantile(returns, alpha)
    # Computing the VaR ratio
    VaRratio.hat <- (mu.hat - rf)/quantile.alpha
    
    # Computing the IF vector for the VaR
    IF.VaRratio <- -(x-mu.hat)/(-quantile.alpha) - (VaRratio.hat/quantile.alpha)*((x<=quantile.alpha)-alpha)/density.fit(quantile.alpha)
    
    # Return value for the risk measure
    return(as.numeric(IF.VaRratio))
  }
}

IF.Rachev.fn <- function(x, returns, parsRachev.IF, alpha=0.1, beta=0.1, rf=0){
  
  # IF for null returns
  if(is.null(returns)){
    # Parameters for null returns
    q.alpha <- parsRachev.IF$q.alpha
    es.alpha <- parsRachev.IF$es.alpha
    q.beta <- parsRachev.IF$q.beta
    eg.beta <- parsRachev.IF$eg.beta
    rach.r <- parsRachev.IF$rach.r
    # IF computation for null returns
    IF.Rachev <- (1/es.alpha)*((x>=q.beta)*(x-q.beta)/beta + q.beta - eg.beta) - 
      (rach.r/es.alpha)*(-(x<=q.alpha)*(x-q.alpha)/alpha - q.alpha - es.alpha)
    # Return value for null returns
    return(IF.Rachev)
    
  } else{
    # Storing the dates
    if(xts::is.xts(returns))
      returns.dates <- zoo::index(returns)
    
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
    IF.Rachev <- (1/ES.lower)*((1/beta)*(x-rf-quantile.upper)*(x>=quantile.upper)+quantile.upper-ES.upper) - 
      (ES.upper/ES.lower^2)*(1/alpha * (-x + rf - quantile.lower)*(x<=quantile.lower)+quantile.lower-ES.lower)
    
    # Return value for the risk measure
    return(as.numeric(IF.Rachev))
  }
}


IF.LPM.fn <- function(x, returns, parsLPM.IF, const=0, order=1){
  
  # Nuisance paramters if returns null
  if(is.null(returns)){
    if(order==1){
      # Computing the IF vector for LPM (order=1 case)
      IF.LPM <- (const - x)*(x <= const) - parsLPM.IF$lpm1
    } else if (order==2){
      
      # Computing the IF vector for LPM (order=2 case)
      IF.LPM <- (x - const)^2 * (x <= const)
      LPM.stored <- parsLPM.IF$lpm2
      IF.LPM <- IF.LPM - LPM.stored^2
    }
    # Return value for the risk measure
    return(IF.LPM)
    
  } else{
    if(order==1){
      # Computing the IF vector for LPM (order=1 case)
      IF.LPM <- (const - x)*(x <= const) - LPM(x, const=const, order=order)
      } else if (order==2){
    
      # Computing the IF vector for LPM (order=2 case)
      IF.LPM <- (x - const)^2 * (x <= const)
      LPM.stored <- LPM(x, const=const, order=order)
      IF.LPM <- IF.LPM - LPM.stored^2
      IF.LPM <- IF.LPM / 2 / LPM.stored
      }
    # Return value for the risk measure
    return(IF.LPM)
    }
}


IF.OmegaRatio.fn <- function(x, returns, parsOmega.IF, const=0){
  
  # IF for null returns
  if(is.null(returns)){
    # Parameters for null returns
    lpm1 <- parsOmega.IF$lpm1
    upm1 <- parsOmega.IF$upm1
    omega <- parsOmega.IF$omega
    # IF computation for null returns
    IF.Omega <- 1/lpm1*((x-const)*(x>=const)-upm1) - omega/lpm1*((const-x)*(x<=const)-lpm1)
    # Return value for null returns
    return(IF.Omega)
  } else{
    # Returning length of returns vector
    N <- length(returns)
    # Computing Partial moments 
    LPM.c <- LPM(returns, const=const, order=1) 
    UPM.c <- UPM(returns, const=const, order=1)
    
    # Computing the IF for the Omega Ratio
    IF.Omega <- 1/LPM.c*((x-const)*(x>=const)-UPM.c) - ((UPM.c/LPM.c)/LPM.c)*((const-x)*(x<=const)-LPM.c)
    
    # Return value for the risk measure
    return(IF.Omega)
  }
}


IF.SemiSD.fn <- function(x, returns, parsSemiSD.IF, rf=0){
  
  if(is.null(returns)){
    mu.hat <- parsSemiSD.IF$mu
    SemiSD <- parsSemiSD.IF$semisd
    smean <- parsSemiSD.IF$smean
    # Return value if nuisance parameters are given
    IF.SemiSD <- ((x - mu.hat)^2 * (x <= mu.hat) - 2*smean * (x-mu.hat) - SemiSD^2)/(2*SemiSD)
    return(IF.SemiSD)
  } else{
    # Computing the mean of the returns
    mu.hat <- mean(returns)
    # Computing SD- of the returns
    sigma.minus.hat <- sqrt(mean((returns-mu.hat)^2*(returns<=mu.hat)))
  }
  
  # Computing the IF vector for SemiSD
  IF.SemiSD <- (x - mu.hat)^2 * (x <= mu.hat)
  IF.SemiSD <- IF.SemiSD - 2 * mean((x-mu.hat) * (x <= mu.hat)) * (x - mu.hat)
  IF.SemiSD <- IF.SemiSD - sigma.minus.hat^2
  IF.SemiSD <- IF.SemiSD / 2 / sigma.minus.hat
  
  # Return value for the risk measure
  return(IF.SemiSD)
}















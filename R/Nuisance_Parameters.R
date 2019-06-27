#' @title Nuisance Parameters Computation
#' 
#' @description \code{nuis.pars} returns the value of the nuisance parameters used in the evaluation of the shape of influence functions for risk and performance measures.
#'
#' @param mu Mean parameter.
#' @param sd Standard deviation parameter.
#' @param c Constant value for threshold.
#' @param alpha Parameters for the lower tail quantile.
#' @param beta Parameter for the upper tail quantile.
#'
#' @return List of nuisance parameters.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @export
#'
#' @examples
#' # Nuisance parameters using default values
#' defaultNuisance <- nuis.parsFn()
#'
#' # Nuisance parameters using specified values
#' specifiedNuisance <- nuis.parsFn(mu=0.02, sd=0.1, c=0.01, alpha=0.05, beta=0.1)
#'

nuis.parsFn <- function(mu=0.01, sd=0.05, c=0, alpha=0.1, beta=0.1){
  
  # Check format of input
  if(!is.numeric(c(mu, sd, c, alpha, beta)))
    stop("All function arguments must be numeric.")
  
  # Check input for alpha parameter
  if(alpha < 0 || alpha >1)
    stop("The argument alpha must be a numeric value between 0 and 1.")
  
  # Check input for beta parameter
  if(beta < 0 || beta >1)
    stop("The argument beta must be a numeric value between 0 and 1.") 
    
  # Nuisance Parameters for Risk Estimators
  
  # SSD
  ssd <- sd/sqrt(2)
  smean <- -dnorm(0)*sd

  # LPM1
  d <- (c-mu)/sd
  lpm1 <- (d*pnorm(d)+dnorm(d))*sd

  # LPM2
  d <- (c-mu)/sd
  lpm2 <- ((d^2+1)*pnorm(d)+d*dnorm(d))*sd^2

  # ES
  q.alpha <- mu+qnorm(alpha)*sd
  es.alpha <- -mu+dnorm(qnorm(alpha))/alpha*sd

  # VaR
  q.alpha <- mu+qnorm(alpha)*sd
  fq.alpha <- dnorm(q.alpha, mu, sd)

  # Nuisance Parameters for Performance Estimators

  # SR 
  sr <- mu/sd   

  # SoR_c
  d <- (c-mu)/sd
  lpm2 <- ((d^2+1)*pnorm(d)+d*dnorm(d))*sd^2
  sor.c <- mu/sqrt(lpm2)

  # SoR_mu
  ssd <- sd/sqrt(2)
  smean <- -dnorm(0)*sd
  sor.mu <- mu/ssd

  # ESratio
  q.alpha <- mu+qnorm(alpha)*sd
  es.alpha <- -mu+dnorm(qnorm(alpha))/alpha*sd
  ES.ratio <- mu/es.alpha

  # VaRratio
  q.alpha <- mu+qnorm(alpha)*sd
  VaR.ratio <- -mu/q.alpha
  fq.alpha <- dnorm(q.alpha,mu,sd) 

  # RachR
  q.alpha <- mu+qnorm(alpha)*sd
  es.alpha <- -mu+dnorm(qnorm(alpha))/alpha*sd
  q.beta <- mu+qnorm(1-beta)*sd
  eg.beta <- mu+dnorm(qnorm(1-beta))/beta*sd
  rach.r <- eg.beta/es.alpha

  # Omega
  d <- (c-mu)/sd
  lpm1 <- (d*pnorm(d)+dnorm(d))*sd
  upm1 <- lpm1+mu-c
  omega <- upm1/lpm1
  
  # Return list for nuisance parameters
  return(list(mu=mu, sd=sd, c=c, alpha=alpha, beta=beta, ssd=ssd, smean=smean, lpm1=lpm1, lpm2=lpm2, 
              q.alpha=q.alpha, es.alpha=es.alpha, fq.alpha=fq.alpha, sr=sr, sor.c=sor.c, sor.mu=sor.mu,
              es.alpha=es.alpha, ES.ratio=ES.ratio, VaR.ratio=VaR.ratio, q.beta=q.beta, eg.beta=eg.beta, 
              rach.r=rach.r, omega=omega, upm1=upm1))
}
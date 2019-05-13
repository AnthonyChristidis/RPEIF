#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"

// [[Rcpp::export]]
arma::vec IF_mean(const arma::vec & returns){
  
  // Computing the mean of the returns
  double mu_hat = mean(returns);
  
  // Computing the IF vector for the mean of the returns
  arma::vec IF_mean_vector = returns - mu_hat;
  
  // Returning the IF vector for the mean of the returns
  return(IF_mean_vector);
}

// [[Rcpp::export]]
arma::vec IF_SD(const arma::vec & returns){
  
  // Computing hte mean of the returns
  double mu_hat = arma::mean(returns);
  // Computing the standard deviation of the returns
  double sd_hat = arma::stddev(returns);
  
  // Computing the IF vector for the standard deviation
  arma::vec IF_SD_vector = (pow((returns-mu_hat),2)-pow(sd_hat,2))/2/sd_hat;
  
  // Returning the IF vector for the standard deviation
  return(IF_SD_vector);
}

// [[Rcpp::export]]
arma::vec IF_SR(const arma::vec & returns, const double & rf){
  
  // Computing the mean of the returns
  double mu_hat = arma::mean(returns);
  // Computing the SD of the returns
  double sd_hat = arma::stddev(returns);
  //Computin the IF vector for the SR
  arma::vec IF_SR_vector = 1/sd_hat*(returns-mu_hat)-1/2*mu_hat/pow(sd_hat,3)*(pow((returns-mu_hat),2)-pow(sd_hat,2));
  
  // Returning the IF vector for the SR
  return(IF_SR_vector);
}

// [[Rcpp::export]]
arma::vec IF_SoR_mean(const arma::vec & returns, const double & rf){
  
  // Computing the mean of the returns
  double mu_hat = arma::mean(returns);
  // Computing the SD of the returns
  double sigma_hat = sqrt(mean(pow((returns-mu_hat),2)));
  
  // Computing SD- of the returns
  double sigma_minus_hat = sqrt(mean(pow((returns-mu_hat),2)%(returns<=mu_hat)));
  // Computing the SoR estimate
  double SoR_hat = (mu_hat-rf)/sigma_minus_hat;
  // Computing the mean parameter for the SoR
  double mu1_minus_hat = mean((returns-mu_hat)%(returns<=mu_hat));
  
  // Computing the IF vector for SoR_M
  arma::vec IF_SoR_M_vector = -SoR_hat/2/pow(sigma_minus_hat,2)*pow((returns-mu_hat-rf),2)%(returns<=mu_hat)+
    (1/sigma_minus_hat+SoR_hat*mu1_minus_hat/pow(sigma_minus_hat,2))*(returns-mu_hat-rf)+
    SoR_hat/2;
  
  // Returning the IF vector for SoR_M
  return(IF_SoR_M_vector);
}

// [[Rcpp::export]]
arma::vec IF_OmegaRatio(const arma::vec & returns, const double & cst){
  
  // Returning length of returns vector
  arma::uword N = returns.n_elem;
  
  // Computing Omega+
  double Omega_p = sum(returns(returns>=cst)-cst)/N;
  // Computing Omega-
  double Omega_m = sum(cst-returns(returns<=cst))/N; 
  // Computing the IF vector for Omega Ratio
  arma::vec IF_Omega_vector = ((returns - cst) % (returns >= cst) - Omega_p)/ Omega_m;
  IF_Omega_vector = IF_Omega_vector - Omega_p / pow(Omega_m,2) * ((cst - returns) % (returns <= cst) - Omega_m);
  
  // Returning the IF vector for Omega Ratio
  return(IF_Omega_vector);
}

// [[Rcpp::export]]
arma::vec IF_SSD(const arma::vec & returns, const double & rf){
  
  // Computing the mean of the returns
  double mu_hat = arma::mean(returns);
  // Computing the SD of the returns
  double sigma_hat = sqrt(mean(pow((returns-mu_hat),2)));
  // Computing SD- of the returns
  double sigma_minus_hat = sqrt(mean(pow((returns-mu_hat),2)%(returns<=mu_hat)));
  // Computing the Sortino Ratio
  double SoR_hat = (mu_hat-rf)/sigma_minus_hat;
  // Computing the IF vector for SSD
  arma::vec IF_SSD_vector = pow((returns - mu_hat),2) % (returns <= mu_hat);
  IF_SSD_vector = IF_SSD_vector - 2 * mean((returns-mu_hat) % (returns <= mu_hat)) * (returns - mu_hat);
  IF_SSD_vector = IF_SSD_vector - pow(sigma_minus_hat,2);
  IF_SSD_vector = IF_SSD_vector / 2 / sigma_minus_hat;
  
  // Returning the IF vector for SSD
  return(IF_SSD_vector);
}





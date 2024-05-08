//library('Rcpp')
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

#include <RcppArmadillo.h>

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>
#include <ctime>
#include <Rcpp.h>


// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
using namespace std;






//*******************************************************************//
//                MAIN FUNC                        //
//*******************************************************************//
//' MAPLE
//' 
//' @export
// [[Rcpp::export]]

SEXP MAPLE_CPP(SEXP zscorexin, SEXP zscoreyin, SEXP Sigmain, SEXP samplen1, SEXP samplen2, SEXP Gibbsnumberin, SEXP burninproportion, SEXP initial_betain, SEXP pi_beta_shape_in, SEXP pi_beta_scale_in,
                 SEXP pi_c_shape_in, SEXP pi_c_scale_in, SEXP pi_1_shape_in, SEXP pi_1_scale_in,SEXP pi_0_shape_in, SEXP pi_0_scale_in, SEXP sigma2x_in, SEXP sigma2y_in, SEXP rho_in){// *
  try{
    const int Gibbs_number = Rcpp::as<int>(Gibbsnumberin);
    const double burnin_p = Rcpp::as<double>(burninproportion);
    const int n1 = Rcpp::as<int>(samplen1);  
    const int n2 = Rcpp::as<int>(samplen2); 
    const double lamda_beta1 = Rcpp::as<double>(pi_beta_shape_in); 
    const double lamda_beta2 = Rcpp::as<double>(pi_beta_scale_in); 
    const double lamda_c1 = Rcpp::as<double>(pi_c_shape_in); 
    const double lamda_c2 = Rcpp::as<double>(pi_c_scale_in); 
    const double lamda_21 = Rcpp::as<double>(pi_1_shape_in); 
    const double lamda_22 = Rcpp::as<double>(pi_1_scale_in); 
    const double lamda_31 = Rcpp::as<double>(pi_0_shape_in); 
    const double lamda_32 = Rcpp::as<double>(pi_0_scale_in); 
    const double sigma2x = Rcpp::as<double>(sigma2x_in);
    const double sigma2y = Rcpp::as<double>(sigma2y_in);
    const double rho = Rcpp::as<double>(rho_in);
    const arma::vec zscorex = as<arma::vec>(zscorexin);
    const arma::vec zscorey = as<arma::vec>(zscoreyin);
    const arma::mat Sigma = as<arma::mat>(Sigmain);
    
    const arma::vec initial_beta = as<arma::vec>(initial_betain);
    
    int p1 = zscorex.n_elem, p2 = zscorey.n_elem;
    if (p1 != p2){
      perror("The dimensions of x1 and x2 are not matched");
    }
    int p = p1;
    
    vec zscorexh=zscorex, zscoreyh=zscorey;
    
    double sigmaz1 = 2.0/p;
    
    double pi = 0.1;
    //pi means pi_beta
    double pi_2 = 0.05;
    //pi_2 mean pi_c 
    double pi_gamma_21 = 0.25;
    //pi_1 
    double pi_gamma_31 = 0.005;
    //pi_0
    double sigma2gamma_1 = 1.0/p;
    
    double w=0.1414;

    double alpha = 0;  
    
    vec latent_beta = initial_beta;  
    
    vec pleiotropy_gamma = zeros<vec>(p);
    
    mat Identity = eye<mat>(p,p);
    
    int burnin = ceil(burnin_p*Gibbs_number);
    vec sample_alpha(Gibbs_number);
    vec sample_pi(Gibbs_number);
    vec sample_pi_2(Gibbs_number);
    vec sample_sigmaz1(Gibbs_number);
    vec sample_sigma2gamma_1(Gibbs_number);
    vec sample_w(Gibbs_number);
    vec sample_pi_gamma_21(Gibbs_number);
    vec sample_pi_gamma_31(Gibbs_number);
   
    double sample_w_variance;
    double sample_alpha_variance;
    double U_beta_1,U_gamma_1;
    vec latent_gamma = zeros<vec>(p);
    vec latent_gamma_pleiotropy = zeros<vec>(p);
    vec I = zeros<vec>(p);
    double part_1,part_2,Probability_gamma,randu_number;
    double part_1_pleiotropy,part_2_pleiotropy,Probability_gamma_pleiotropy,randu_number_pleiotropy;
    double part_1_I,part_2_I,Probability_I_indicator,randu_number_I;
    
    double post_sigmaz1_scale;
    double post_sigmaz1_shape;
    double post_sigma2gamma_1_scale;
    double post_sigma2gamma_1_shape;
    
    double sigmaz1_prior_shape=p/10+1;
    double sigmaz1_prior_scale=0.2;
    double sigma2gamma_1_prior_shape=p/5.0+1;
    double sigma2gamma_1_prior_scale=0.2;
    mat left = Sigma - Identity;
    
    
    for(int m=0; m<Gibbs_number; ++m){
      
      double sample_var_pleiotropy_gamma_1 = 1.0/((n2-1)/((1-rho*rho)*sigma2y) + 1/sigma2gamma_1);
    
      for(int k=0; k<p; ++k){
        
        double sample_var_1 = 1.0/((n1-1)/((1-rho*rho)*sigma2x) + (n2-1)*((alpha+w*I(k))*(alpha+w*I(k))/((1-rho*rho)*sigma2y)) -
                              2*rho*sqrt(n1-1)*sqrt(n2-1)*(alpha+w*I(k))/((1-rho*rho)*sqrt(sigma2x)*sqrt(sigma2y))+ 1.0/sigmaz1);
        
        U_beta_1 = ((sqrt(n1-1)*zscorexh(k)-(n1-1)*sum(left.col(k)%latent_beta))/((1-rho*rho)*sigma2x)+((alpha+w*I(k))*(sqrt(n2-1)*zscoreyh(k)-
          (n2-1)*alpha*sum(left.col(k)%latent_beta)-(n2-1)*w*sum(left.col(k)%(latent_beta%I))-(n2-1)*sum(Sigma.col(k)%pleiotropy_gamma)))/((1-rho*rho)*sigma2y)-
         rho*(alpha+w*I(k))*(sqrt(n2-1)*zscorexh(k)-sqrt(n1-1)*sqrt(n2-1)*sum(left.col(k)%latent_beta))/((1-rho*rho)*sqrt(sigma2x)*sqrt(sigma2y))-
          rho*(sqrt(n1-1)*zscoreyh(k)-sqrt(n1-1)*sqrt(n2-1)*alpha*sum(left.col(k)%latent_beta)-sqrt(n1-1)*sqrt(n2-1)*w*sum(left.col(k)%(latent_beta%I))-
          sqrt(n1-1)*sqrt(n2-1)*sum(Sigma.col(k)%pleiotropy_gamma))/((1-rho*rho)*sqrt(sigma2x)*sqrt(sigma2y)))* sample_var_1;
                      	
        part_1 = exp(0.5* U_beta_1*U_beta_1/sample_var_1 + 0.5*log(sample_var_1) - 0.5*log(sigmaz1)+log(pi)+ latent_gamma_pleiotropy(k)*log(pi_gamma_21)+
          (1-latent_gamma_pleiotropy(k))*log(1-pi_gamma_21)+ I(k)*log(pi_2)+ (1-I(k))*log(1-pi_2));
        
        part_2 = exp(log(1-pi)+ latent_gamma_pleiotropy(k)*log(pi_gamma_31)+ (1-latent_gamma_pleiotropy(k))*log(1-pi_gamma_31));
        
        Probability_gamma = part_1/(part_1+part_2);
        
        randu_number = as_scalar(randu(1));
        
        if(randu_number <= Probability_gamma){
          latent_gamma(k) = 1;
          latent_beta(k) = as_scalar(randn(1)*sqrt(sample_var_1)+ U_beta_1);
        } else {
          latent_gamma(k) = 0;
          latent_beta(k) = 0;
        }
        
        U_gamma_1 = ((sqrt(n2-1)*zscoreyh(k)-(n2-1)*alpha*sum(Sigma.col(k)%latent_beta)-(n2-1)*w*sum(Sigma.col(k)%(latent_beta%I))-
          (n2-1)*sum(left.col(k)%pleiotropy_gamma))/((1-rho*rho)*sigma2y)-rho*(sqrt(n2-1)*zscorexh(k)-
          sqrt(n1-1)*sqrt(n2-1)*sum(Sigma.col(k)%latent_beta))/((1-rho*rho)*sqrt(sigma2x)*sqrt(sigma2y))) * sample_var_pleiotropy_gamma_1;
        
        part_1_pleiotropy = exp(0.5* U_gamma_1*U_gamma_1/sample_var_pleiotropy_gamma_1 + 0.5*log(sample_var_pleiotropy_gamma_1) - 0.5*log(sigma2gamma_1)+
          latent_gamma(k)*log(pi_gamma_21)+ (1-latent_gamma(k))*log(pi_gamma_31));
        
        part_2_pleiotropy = exp(latent_gamma(k)*log(1-pi_gamma_21)+ (1-latent_gamma(k))*log(1-pi_gamma_31));
        
        Probability_gamma_pleiotropy = part_1_pleiotropy/(part_1_pleiotropy+part_2_pleiotropy);
        
        randu_number_pleiotropy = as_scalar(randu(1));
        
        if(randu_number_pleiotropy<= Probability_gamma_pleiotropy){
          latent_gamma_pleiotropy(k) = 1;
          pleiotropy_gamma(k) = as_scalar(randn(1)*sqrt(sample_var_pleiotropy_gamma_1)+ U_gamma_1);
         
        } else {
          latent_gamma_pleiotropy(k) = 0;
          pleiotropy_gamma(k) = 0;
        }		
        
        part_1_I = exp(-0.5*(w * w * latent_beta(k)*latent_beta(k)*(n2-1)-2*w*latent_beta(k)*(sqrt(n2-1)*zscoreyh(k)-(n2-1)*alpha*sum(Sigma.col(k)% latent_beta)-
          (n2-1)*w*sum(left.col(k)% (latent_beta % I))- (n2-1)*sum(Sigma.col(k)% pleiotropy_gamma)))/((1-rho*rho)*sigma2y)-rho*sqrt(n2-1)*w*latent_beta(k)*(zscorexh(k)-
          sqrt(n1-1)*sum(Sigma.col(k)%latent_beta))/((1-rho*rho)*sqrt(sigma2x)*sqrt(sigma2y))+ latent_gamma(k)*log(pi_2));
        
        part_2_I = exp(latent_gamma(k)*log(1-pi_2));
        
        Probability_I_indicator = part_1_I/(part_1_I+part_2_I);
        
        randu_number_I = as_scalar(randu(1));
        
        if(randu_number_I<= Probability_I_indicator){
          I(k) = 1;
        } else {
          I(k) = 0;
        }		
      }
      
      double indicator = as_scalar((n2-1)*(latent_beta%I).t()*Sigma*(latent_beta%I))/((1-rho*rho)*sigma2y);
      
      if(indicator >0){
        sample_w_variance = 1.0/(as_scalar((n2-1)*(latent_beta%I).t()*Sigma*(latent_beta%I))/((1-rho*rho)*sigma2y));
        double sample_w_mean = as_scalar(sqrt(n2-1)*(zscoreyh-sqrt(n2-1)*Sigma*latent_beta*alpha-sqrt(n2-1)*Sigma*pleiotropy_gamma).t()*(latent_beta%I)/((1-rho*rho)*sigma2y)-
                                sqrt(n2-1)*rho*(zscorexh-sqrt(n1-1)*Sigma*latent_beta).t()*(latent_beta%I)/((1-rho*rho)*sqrt(sigma2x*sigma2y)))*sample_w_variance;
        
        sample_w(m) = as_scalar(randn(1)*sqrt(sample_w_variance)+ sample_w_mean);
        
      } else {
        sample_w(m) = 0;
      }
      
      w=sample_w(m); 
     
      double indicator_alpha = as_scalar((n2-1)*latent_beta.t()* Sigma * latent_beta)/((1-rho*rho)*sigma2y);
      
      if(indicator_alpha > 0){
        sample_alpha_variance = 1.0/(as_scalar((n2-1)*latent_beta.t()* Sigma * latent_beta)/((1-rho*rho)*sigma2y));
        double sample_alpha_mean = as_scalar((sqrt(n2-1)*latent_beta.t()* zscoreyh-(n2-1)*latent_beta.t()*Sigma*pleiotropy_gamma-
                                    (n2-1)*latent_beta.t()*Sigma*(latent_beta % I)*w)/((1-rho*rho)*sigma2y)-rho*(sqrt(n2-1)*latent_beta.t()*zscorexh-
                                    sqrt(n1-1)*sqrt(n2-1)*latent_beta.t()*Sigma*latent_beta)/((1-rho*rho)*sqrt(sigma2x*sigma2y)))*sample_alpha_variance;
        
        sample_alpha(m) = as_scalar(randn(1)*sqrt(sample_alpha_variance)+ sample_alpha_mean);

      } else {
        sample_alpha(m) = 0;
        
      }
      
      alpha=sample_alpha(m);
      
      double shape1 = sum(latent_gamma)+lamda_beta1;
      
      double shape2 = sum(1.0-latent_gamma)+lamda_beta2;
      
      sample_pi(m)= r_4beta(shape1,shape2,0,1);
      
      pi = sample_pi(m);
      
      double shape1_pleiotropy_21 = sum(latent_gamma % latent_gamma_pleiotropy)+lamda_21;
      
      double shape2_pleiotropy_21 = sum(latent_gamma % (1.0-latent_gamma_pleiotropy))+lamda_22;
      
      sample_pi_gamma_21(m)= r_4beta(shape1_pleiotropy_21,shape2_pleiotropy_21,0,1);
      
      pi_gamma_21 = sample_pi_gamma_21(m);
      
      double shape1_pleiotropy_31 = sum((1-latent_gamma) % latent_gamma_pleiotropy)+lamda_31;
      
      double shape2_pleiotropy_31 = sum((1-latent_gamma) % (1.0-latent_gamma_pleiotropy))+lamda_32;
      
      sample_pi_gamma_31(m)= r_4beta(shape1_pleiotropy_31,shape2_pleiotropy_31,0,1);
      
      pi_gamma_31 = sample_pi_gamma_31(m);
      
      
      
      double shape1_pi_2 = sum(latent_gamma % I)+lamda_c1;
      
      double shape2_pi_2 = sum(latent_gamma % (1.0-I))+lamda_c2;
      
      sample_pi_2(m)= r_4beta(shape1_pi_2,shape2_pi_2,0,1);
      
      pi_2 = sample_pi_2(m);
      
      post_sigmaz1_scale = 2/(as_scalar(sum(latent_gamma % latent_beta % latent_beta)+2*sigmaz1_prior_scale));
      
      post_sigmaz1_shape = 0.5*sum(latent_gamma)+ sigmaz1_prior_shape;
      
      post_sigmaz1_scale = abs(post_sigmaz1_scale);
      
      sample_sigmaz1(m) = 1.0/as_scalar(randg(1, distr_param(post_sigmaz1_shape,post_sigmaz1_scale)));
      
      sigmaz1 = sample_sigmaz1(m);
      
      post_sigma2gamma_1_scale = 2/(as_scalar(sum(latent_gamma_pleiotropy % pleiotropy_gamma % pleiotropy_gamma)+2*sigma2gamma_1_prior_scale));
      
      post_sigma2gamma_1_shape = 0.5*sum(latent_gamma_pleiotropy)+ sigma2gamma_1_prior_shape;
      
      post_sigma2gamma_1_scale = abs(post_sigma2gamma_1_scale);
      
      sample_sigma2gamma_1(m) = 1.0/as_scalar(randg(1, distr_param(post_sigma2gamma_1_shape,post_sigma2gamma_1_scale)));
      
      sigma2gamma_1 = sample_sigma2gamma_1(m);
      
    } 
    
    
    double alpha_estimate = mean(sample_alpha.subvec((burnin-1),(Gibbs_number-1)));
    
    double alpha_sd = stddev(sample_alpha.subvec((burnin-1),(Gibbs_number-1)));
    
    double w_estimate = mean(sample_w.subvec((burnin-1),(Gibbs_number-1))); 
    
    double sigmaz1_estimate = mean(sample_sigmaz1.subvec((burnin-1),(Gibbs_number-1)));
    
    double sigma2gamma_1_estimate = mean(sample_sigma2gamma_1.subvec((burnin-1),(Gibbs_number-1)));
    
    return List::create(Rcpp::Named("alpha") = alpha_estimate,
                        Rcpp::Named("w") = w_estimate,
                        Rcpp::Named("sigmabeta") = sigmaz1_estimate,
                        Rcpp::Named("sigmaeta") = sigma2gamma_1_estimate,
                        Rcpp::Named("alpha_sd") = alpha_sd);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)..." );
  }
  return R_NilValue;
}// end func

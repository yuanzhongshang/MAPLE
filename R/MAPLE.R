#' @title The function of MAPLE method Mendelian Randomization 
#' @description  MAPLE is an efficient Mendelian randomization method with self-adaptive determination of sample structure and multiple pleiotropic effects
#' @param Zscore_1 the Zscore vector of the SNP effect size vector for the exposure
#' @param Zscore_2 the Zscore vector of the SNP effect size vector for the outcome 
#' @param Sigma1in the LD matrix for the SNPs in the exposure GWAS data
#' @param Sigma2in the LD matrix for the SNPs in the outcome GWAS data
#' @param samplen1 the sample size of exposure GWAS 
#' @param samplen2 the sample size of outcome GWAS 
#' @param Gibbsnumber the number of Gibbs sampling iterations with the default to be 1000
#' @param burninproportion  the proportion to burn in from Gibbs sampling iterations, with default to be 0.2 
#' @param pi_beta_shape the prior shape paramter for pi_beta with the default to be 0.5
#' @param pi_beta_scale the prior scale paramter for pi_beta with the default to be 4.5
#' @param pi_c_shape the prior shape paramter for pi_c with the default to be 0.5
#' @param pi_c_scale the prior shape paramter for pi_c with the default to be 9.5
#' @param pi_1_shape the prior shape paramter for pi_1 with the default to be 0.5
#' @param pi_1_scale the prior scale paramter for pi_1 with the default to be 1.5
#' @param pi_0_shape the prior shape paramter for pi_0 with the default to be 0.05
#' @param pi_0_scale the prior scale paramter for pi_0 with the default to be 9.95
#' @param t1 the t1 estimated from exposure data by LDSC
#' @param t2 the t2 estimated from outcome data by LDSC
#' @param t12 the t12 estimated from both exposure and outcome data by LDSC
#' 
#' @return A list of estimated parameters including the p values for the causal effect test 
#' \item{causal_effect}{The estimate of causal effect}
#' \item{causal_pvalue}{The p value for the causal effect}
#' \item{correlated_pleiotropy_effect}{The confounder effect on the outcome}
#' \item{sigmaeta}{The variance estimate for the uncorrelated pleiotropy effect}
#' \item{sigmabeta}{The variance estimate for the SNP effect sizes on the exposure}
#' 
#' 
#' @export

MAPLE <- function(Zscore_1,Zscore_2,Sigma1in,Sigma2in,samplen1,samplen2,Gibbsnumber=1000,burninproportion=0.2,pi_beta_shape=0.5,
pi_beta_scale=4.5,pi_c_shape=0.5,pi_c_scale=9.5,pi_1_shape=0.5,pi_1_scale=1.5,pi_0_shape=0.05,pi_0_scale=9.95,t1,t2,t12){

initial_betain = rep(0,length(Zscore_1))
Sigmain_est = (samplen1*Sigma1in+samplen2*Sigma2in)/(samplen1+samplen2)
sigma2x_est = t1
sigma2y_est = t2
rho_est = t12/sqrt(t1*t2)

if (abs(rho_est) >= 1) {
if (rho_est >= 1){
rho_est = 0.99
} else {
rho_est = -0.99
}
}     
   
re = MAPLE_CPP(zscorexin=Zscore_1,zscoreyin=Zscore_2,Sigmain=Sigmain_est,samplen1,samplen2,Gibbsnumberin=Gibbsnumber,burninproportion=burninproportion,initial_betain=initial_betain,
pi_beta_shape_in=pi_beta_shape,pi_beta_scale_in=pi_beta_scale,pi_c_shape_in=pi_c_shape,pi_c_scale_in=pi_c_scale,pi_1_shape_in=pi_1_shape,pi_1_scale_in=pi_1_scale,
pi_0_shape_in=pi_0_shape,pi_0_scale_in=pi_0_scale,sigma2x_in=sigma2x_est,sigma2y_in=sigma2y_est,rho_in=rho_est)

while(re$alpha == 0){
if (rho_est > 0) {
rho_est = rho_est-0.01
}
if (rho_est < 0) {
rho_est = rho_est + 0.01
}
re = MAPLE_CPP(zscorexin=Zscore_1,zscoreyin=Zscore_2,Sigmain=Sigmain_est,samplen1,samplen2,Gibbsnumberin=Gibbsnumber,burninproportion=burninproportion,initial_betain=initial_betain,
pi_beta_shape_in=pi_beta_shape,pi_beta_scale_in=pi_beta_scale,pi_c_shape_in=pi_c_shape,pi_c_scale_in=pi_c_scale,pi_1_shape_in=pi_1_shape,pi_1_scale_in=pi_1_scale,
pi_0_shape_in=pi_0_shape,pi_0_scale_in=pi_0_scale,sigma2x_in=sigma2x_est,sigma2y_in=sigma2y_est,rho_in=rho_est)
}

pvalue_alpha = 2*(1-pnorm(abs(re$alpha/re$alpha_sd)))
result = list()
result$causal_effect = re$alpha
result$causal_pvalue = pvalue_alpha
result$cause.se = re$alpha_sd
result$correlated_pleiotropy_effect = re$w
result$sigmabeta = re$sigmabeta
result$sigmaeta = re$sigmaeta
return(result)
}

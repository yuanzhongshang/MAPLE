% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MAPLE.R
\name{MAPLE}
\alias{MAPLE}
\title{The function of MAPLE method Mendelian Randomization}
\usage{
MAPLE(
  Zscore_1,
  Zscore_2,
  Sigma1in,
  Sigma2in,
  samplen1,
  samplen2,
  Gibbsnumber = 1000,
  burninproportion = 0.2,
  pi_beta_shape = 0.5,
  pi_beta_scale = 4.5,
  pi_c_shape = 0.5,
  pi_c_scale = 9.5,
  pi_1_shape = 0.5,
  pi_1_scale = 1.5,
  pi_0_shape = 0.05,
  pi_0_scale = 9.95,
  t1,
  t2,
  t12
)
}
\arguments{
\item{Zscore_1}{the Zscore vector of the SNP effect size vector for the exposure}

\item{Zscore_2}{the Zscore vector of the SNP effect size vector for the outcome}

\item{Sigma1in}{the LD matrix for the SNPs in the exposure GWAS data}

\item{Sigma2in}{the LD matrix for the SNPs in the outcome GWAS data}

\item{samplen1}{the sample size of exposure GWAS}

\item{samplen2}{the sample size of outcome GWAS}

\item{Gibbsnumber}{the number of Gibbs sampling iterations with the default to be 1000}

\item{burninproportion}{the proportion to burn in from Gibbs sampling iterations, with default to be 0.2}

\item{pi_beta_shape}{the prior shape paramter for pi_beta with the default to be 0.5}

\item{pi_beta_scale}{the prior scale paramter for pi_beta with the default to be 4.5}

\item{pi_c_shape}{the prior shape paramter for pi_c with the default to be 0.5}

\item{pi_c_scale}{the prior shape paramter for pi_c with the default to be 9.5}

\item{pi_1_shape}{the prior shape paramter for pi_1 with the default to be 0.5}

\item{pi_1_scale}{the prior scale paramter for pi_1 with the default to be 1.5}

\item{pi_0_shape}{the prior shape paramter for pi_0 with the default to be 0.05}

\item{pi_0_scale}{the prior scale paramter for pi_0 with the default to be 9.95}

\item{t1}{the t1 estimated from exposure data by LDSC}

\item{t2}{the t2 estimated from outcome data by LDSC}

\item{t12}{the t12 estimated from both exposure and outcome data by LDSC}
}
\value{
A list of estimated parameters including the p values for the causal effect test 
\item{causal_effect}{The estimate of causal effect}
\item{causal_pvalue}{The p value for the causal effect}
\item{correlated_pleiotropy_effect}{The confounder effect on the outcome}
\item{sigmaeta}{The variance estimate for the uncorrelated pleiotropy effect}
\item{sigmabeta}{The variance estimate for the SNP effect sizes on the exposure}
}
\description{
MAPLE is an efficient Mendelian randomization method with self-adaptive determination of sample structure and multiple pleiotropic effects
}

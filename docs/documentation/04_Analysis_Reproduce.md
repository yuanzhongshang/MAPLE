---
layout: page
title: Analysis Reproduce
description: ~
---

# Generate the simulation data

We randomly selected 100 SNPs to have non-zero effects on the exposure. We obtained genotypes from UK Biobank on chromosome 1, standardized each SNP to have zero mean and unit variance.For null simulations, we performed 500 simulation replicates in each setting to examine type I error control. As the same p-values from different methods may correspond to different type I errors, we evaluated the power based on a false discovery rate (FDR) of 0.05 for fair comparison. Specifically, we conducted 100 alternative simulations along with 400 null simulations, repeating this analysis five times to calculate the average power across these replicates. 

## Generate the individual-level data 

```
##situation: 50% sample overlap with horizontal pleiotropy
library(data.table)
library(MASS)
library(BEDMatrix) 
###parameter setting
nx=50000
ny=50000
PVE_zx = 0.1 
#null simulation: PVE_alpha = 0; alternative simulation: PVE_alpha = 0.0015
PVE_alpha = 0 
PVE_u = 0.05
PVE_c = 0.0025*PVE_zx
pi_1 = 0.2
pi_c = 0.05
w = sqrt(0.02)
K = 100
#correlation between exposure and outcome
rho = 0.4 
#load genotype data of exposure without NA; data format: plink binary (.bim/.fam/.bed)
gene = BEDMatrix(paste0("./exposure_genotype")) 
#load genotype data of outcome without NA; data format: plink binary (.bim/.fam/.bed)
gene2 = BEDMatrix(paste0("./outcome_genotype")) 
#note: the genotype of exposure and outcome have 50% same individual and have same SNPs
seqsum = seq(from = 1, to = ncol(gene), by = 1)
seq1 = sample(seqsum, K, replace = F) 
seq2 = sample(setdiff(seqsum, seq1), 100-pi_1*K, replace = F)
seq3 = sample(seq1, pi_1*K, replace = F)
sequ = c(seq2,seq3)
seqc_n = sample(c(1:K), pi_c*K, replace = F)

Gx0 = gene[, seq1]
Gx1 = scale(Gx0)
beta = rnorm(K, 0, sqrt(PVE_zx/K))

Gy0 = gene2[,seq1]
Gy1 = scale(Gy0)
alpha =sqrt(PVE_alpha/PVE_zx)

betac = beta[seqc_n]
Gyc = Gy1[,seqc_n]

Gyu = gene2[,sequ]
Gyu = scale(Gyu)
eta_u = rnorm(K, 0, sqrt(PVE_u/K))

sigma2x = 1-PVE_zx
sigma2y = 1-PVE_alpha-PVE_u-PVE_c
sigmaxy = rho*sqrt(sigma2x)*sqrt(sigma2y)
mean = c(0,0) 
sigma = matrix(c(sigma2x, sigmaxy, sigmaxy, sigma2y), nrow=2, ncol=2)
# the residual error of 50% individuals in both exposure and 
# outcome is generated from the multi-normal distribution
epsilon = mvrnorm(nx/2, mean, sigma)
# the residual error of other 50% individuals in exposure or outcome
# is generated from the normal distribution
epsilonx_1 = rnorm(nx/2,0,sqrt(sigma2x))
epsilony_1 = rnorm(ny/2,0,sqrt(sigma2y))
epsilonx_2 = epsilon[,1]
epsilony_2 = epsilon[,2]
epsilonx = c(epsilonx_1,epsilonx_2)
epsilony = c(epsilony_1,epsilony_2)
#generate the phenotype
x = Gx1%*%beta+epsilonx
y = Gy1%*%beta*alpha+Gyc%*%betac*w+Gyu%*%eta_u+epsilony
write.table(x,paste0("./individual_xy/exp.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(y,paste0("./individual_xy/out.txt"),sep="\t",quote=F,row.names=F,col.names=F)
```

## Convert the individual-level data into the summary statistics 

```
library(data.table)
library(BEDMatrix) 
#use Gemma to do association between SNP and exposure
GEMMA="./gemma-0.98.1-linux-static"
pfile=paste0("./individual_xy/exp.txt")
gfile=paste0("./exposure_genotype")
setwd("./gemma_exp")
system(paste0(GEMMA," -bfile ",gfile," -p ",pfile," -lm 1"," -o betax"))
#use Gemma to do association between SNP and outcome
setwd("./gemma_out")
pfile=paste0("./individual_xy/out_",i,".txt")
gfile=paste0("./outcome_genotype")
system(paste0(GEMMA," -bfile ",gfile," -p ",pfile," -lm 1"," -o betay"))
#save zscores of exposure and outcome
exp_output = fread(paste0("./gemma_exp/output/betax.assoc.txt"), header=T)
out_output = fread(paste0("./gemma_out/output/betay.assoc.txt"), header=T)
exp <- exp_output[which(exp_output$p_wald<5e-8),]
exp_re <- exp[which(exp$rs %in% out_output$rs),]
out_re <- out_output[which(out_output$rs %in% exp$rs),]
snplist <- out_re$rs
zscorex <- matrix(exp_re$beta/exp_re$se)
zscorey <- matrix(out_re$beta/out_re$se)
write.table(zscorex,paste0("./zscorex/zscorex.txt"),sep="\t",
            quote=F,row.names=F, col.names=F)
write.table(zscorey,paste0("./zscorey/zscorey.txt"),sep="\t",
            quote=F,row.names=F, col.names=F)
#save LD matrices of exposure and outcome
exposure_data <-  BEDMatrix(paste0("./exposure_genotype"))
bim <- fread(paste0("./exposure_genotype.bim"))
name <- bim[,2]
colnames(exposure_data) <- t(name)
Gx <- exposure_data[,t(snplist)]
Gx_sca  = scale(Gx)
n1=dim(Gx_sca)[1]
Sigmax <- t(Gx_sca)%*%Gx_sca/(n1-1)
write.table(round(Sigmax,5),paste0("./Sigma/Sigmax.txt"),sep="\t",
            quote=F,row.names=F, col.names=F)
outcome_data <-  BEDMatrix(paste0("./outcome_genotype"))
bim <- fread(paste0("./outcome_genotype.bim"))
name <- bim[,2]
colnames(outcome_data) <- t(name)
Gy <- outcome_data[,t(snplist)]
Gy_sca  = scale(Gy)
n2=dim(Gy_sca)[1]
Sigmay <- t(Gy_sca)%*%Gy_sca/(n2-1)
write.table(round(Sigmay,5),paste0("./Sigma/Sigmay.txt"),sep="\t",
            quote=F,row.names=F, col.names=F)
```

# Run MAPLE

```
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(magrittr)
library(data.table)
library(MAPLE)
#load the summary data for exposure and outcome
exp = fread(paste0("./gemma_exp/output/betax.assoc.txt"),head=T)
exp_raw = exp[,c("rs","beta","se","af","allele1","allele0","p_wald","n_obs")] 
colnames(exp_raw) = c("SNP","b","se","frq_A1","A1","A2","P","N")
out = fread(paste0("./gemma_exp/output/betay.assoc.txt"),head=T)
out_raw = out[,c("rs","beta","se","af","allele1","allele0","p_wald","n_obs")] 
colnames(out_raw) = c("SNP","b","se","frq_A1","A1","A2","P","N")
paras = est_SS(dat1 = exp_raw,
               dat2 = out_raw,
               trait1.name = "exp",
               trait2.name = "out",
               ldscore.dir = "./eur_w_ld_chr")
#load the z-score
zscorex = fread(paste0("./zscorex/zscorex.txt"),head=F)
Zscore_1 = as.vector(zscorex[[1]])
zscorey = fread(paste0("./zscorey/zscorey.txt"),head=F)
Zscore_2 = as.vector(zscorey[[1]])
#load the LD matrix
Sigmaxin = fread(paste0("./Sigma/Sigmax.txt"),head=F)
Sigma1in = as.matrix(Sigmaxin)
Sigmayin = fread(paste0("./Sigma/Sigmay.txt"),head=F)
Sigma2in = as.matrix(Sigmayin)
#load the sample size
samplen1 = 50000
samplen2 = 50000
#load the nuisance error parameter
t1 = paras$Omega[1,1]
t2 = paras$Omega[2,2]
t12 = paras$Omega[1,2]
result = MAPLE(Zscore_1,Zscore_2,Sigma1in,Sigma2in,samplen1,samplen2,Gibbsnumber=1000,
               burninproportion=0.2,pi_beta_shape=0.5,pi_beta_scale=4.5,
               pi_c_shape=0.5,pi_c_scale=9.5,pi_1_shape=0.5,pi_1_scale=1.5,
               pi_0_shape=0.05,pi_0_scale=9.95,t1,t2,t12)
```

# Run MRAID

```
library(MRAID)
library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
#load the z-score
zscorex = fread(paste0("./zscorex/zscorex.txt"),head=F)
Zscore_1 = as.vector(zscorex[[1]])
zscorey = fread(paste0("./zscorey/zscorey.txt"),head=F)
Zscore_2 = as.vector(zscorey[[1]])
#load the LD matrix
sigma = fread(paste0("./Sigma/Sigmax.txt"),head=F)
Sigmainx = as.matrix(sigma)
sigma = fread(paste0("./Sigma/Sigmay.txt"),head=F)
Sigmainy = as.matrix(sigma)
Sigma1sin = Sigmainx
Sigma2sin = Sigmainy
#load the sample size
samplen1 = 50000
samplen2 = 50000
result = MRAID(Zscore_1, Zscore_2, Sigma1sin, Sigma2sin, samplen1, samplen2)
```

# Run MR-APSS

```
library(data.table)
library(MRAPSS)
library(readr)

exp <- fread(paste0("./gemma_exp/output/betax.assoc.txt"),head=T)
exp_raw <- exp[,c("rs","beta","se","af","allele1","allele0","p_wald","n_obs")] 
colnames(exp_raw) <- c("SNP","b","se","frq_A1","A1","A2","P","N")
out <- fread(paste0("./gemma_out/output/betay.assoc.txt"),head=T)
out_raw <- out[,c("rs","beta","se","af","allele1","allele0","p_wald","n_obs")] 
colnames(out_raw) <- c("SNP","b","se","frq_A1","A1","A2","P","N")

#2. we use function format_data to obtain the formatted data sets for exposure and outcome
exp = format_data(exp_raw,
                  snp_col = "SNP",
                  b_col = "b",
                  se_col = "se",
                  freq_col = "frq_A1",
                  A1_col = "A1",
                  A2_col = "A2",
                  p_col = "P",
                  n_col = "N")

out = format_data(out_raw,
                  snp_col = "SNP",
                  b_col = "b",
                  se_col = "se",
                  freq_col = "frq_A1",
                  A1_col = "A1",
                  A2_col = "A2",
                  p_col = "P",
                  n_col = "N")
#3.Harmonize the formatted data sets and estimate nuisance parameters
paras = est_paras(dat1 = exp,
                  dat2 = out,
                  trait1.name = "exp",
                  trait2.name = "out",
                  ldscore.dir = "./eur_w_ld_chr")
#4. LD clumping
MRdat =  clump(paras$dat,
               IV.Threshold = 5e-04,
               SNP_col = "SNP",
               pval_col = "pval.exp",
               clump_kb = 1000,
               clump_r2 = 0.01,
               bfile = "./all_1000G_EUR_Phase3",
               plink_bin = "./PLINK/plink")
#5.Fit MR-APSS
MRres = MRAPSS(MRdat,
               exposure = "exp",
               outcome = "out",
               C = paras$C,
               Omega = paras$Omega ,
               Cor.SelectionBias = T)
#MRplot(MRres, exposure = "exp", outcome = "out")
c1 = paras$C[1,1]
c2 =  paras$C[2,2]
c12 = paras$C[1,2]
pi0 = MRres$pi0
beta = MRres$beta
beta.se = MRres$beta.se
p = MRres$pvalue
re = data.frame(beta,beta.se,p,c1,c2,c12,pi0)
```

# Run CAUSE

```
library(data.table)
library(cause)
library(dplyr)
#1.Format Data for CAUSE
X1 <- fread(paste0("./gemma_exp/output/betax.assoc.txt"),head=T)
X1 <- X1[,c("rs","beta","se","allele1","allele0","p_wald")]
colnames(X1) <- c("snp","beta_hat_1","seb1","A1","A2","p1")

X2 <- fread(paste0("./gemma_out/output/betay.assoc.txt"),head=T)
X2 <- X2[,c("rs","beta","se","allele1","allele0","p_wald")]
colnames(X2) <- c("snp","beta_hat_2","seb2","A1","A2","p2")

X <- merge(X1,X2,by=c("snp","A1","A2"))
X <- X[,c("snp","beta_hat_1","seb1","p1","beta_hat_2","seb2","A1","A2","p2")]
cause_data <- new_cause_data(X)
X=cause_data
#2.Step 2: Calculate nuisance parameters
set.seed(i+10000)
varlist <- with(X, sample(snp, size=100000, replace=FALSE))
params <- est_cause_params(X, varlist)
#Step 3: LD Pruning
r2_thresh = 0.01
pval_thresh = 1e-3
X_clump <- X %>%
  rename(rsid = snp,
         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     bfile = "./all_1000G_EUR_Phase3",
                     plink_bin =  "./PLINK/plink", 
                     pop = "EUR")
top_vars <- X_clump$rsid
#Step 4: Fit CAUSE
res <- cause(X=X, variants = top_vars, param_ests = params)
cause_z <- with(res$elpd, z[model1=="sharing" & model2=="causal"])
summary(res, ci_size=0.95)
ci_size=0.95
fit <- res[["causal"]]
ix <- which(fit$params == "gamma")
qs <- with(fit$marge_post[[ix]], step_quantile(c(0.5, (1-ci_size)/2, 1-((1-ci_size)/2)),
                                               begin, end, post))
cause_est <- qs[1]
cause_est_lower <- qs[2]
cause_est_upper <- qs[3]
pvalue<-summary(res)$p
rho_est <- params$rho
result <- c(cause_est,pvalue,cause_z,cause_est_lower,cause_est_upper,rho_est)
result <- as.data.frame(result <- t(result))
colnames(result) <- c("causel_est","pvalue","cause_z","cause_est_lower","cause_est_upper",
                      "rho_est")
```

# Run MR-CUE

```
library(MR.CUE)
library(ggplot2)
library(data.table)
rm(list = ls())
filepan <- vector("list", 22)
NumChr = 1
for(i in 1:NumChr){
  filepan[[i]] <- paste0("UK10KCHR", i, "LDhm3.RDS")
}
exp <- fread(paste0("./gemma_exp/output/betax.assoc.txt"),head=T)
exp_raw <- exp[,c("rs","chr","ps","allele1","allele0","beta","se","p_wald")] 
colnames(exp_raw) <- c("SNP","chr","BP","A1","A2","beta","se","pvalue")
out <- fread(paste0("./gemma_out/output/betay.assoc.txt"),head=T)
out_raw <- out[,c("rs","chr","ps","allele1","allele0","beta","se","p_wald")] 
colnames(out_raw) <- c("SNP","chr","BP","A1","A2","beta","se","pvalue")
write.table(exp_raw,"exp_raw.txt",sep="\t",quote=F,row.names=F, col.names=T)
write.table(out_raw,"out_raw.txt",sep="\t",quote=F,row.names=F, col.names=T)
fileexp = "exp_raw.txt"
fileout = "out_raw.txt"
snpinfo = "UK10Ksnpinforhm3.RDS"
#Step 0. Estimate for the overlape samples
ld_r2_thresh = 0.001
lambad = 0.85
pth = 1.96
RhoEst = EstRho(fileexp, fileout, filepan, snpinfo, ld_r2_thresh, lambad, pth)
rho = mean(RhoEst$Rhores)
#Step 1. Format Data for MR-CUE
pva_cutoff = 5e-4
lambad = 0.85
data <- ReadSummaryStat(fileexp, fileout, filepan, snpinfo, pva_cutoff, lambad)
F4gammah <- data$ResF4gammah
F4Gammah <- data$ResF4Gammah
F4se1 <- data$ResF4se1
F4se2 <- data$ResF4se2
F4Rblock <- data$ResF4Rblock
F4SNPs <- data$ResF4SNPchr
#Step 2. Fit MR-CUE
L = length(F4gammah)
coreNum = 20
opt = list(agm = 0, bgm = 0, atau1 = 0, btau1 = 0,
           atau2 = 0, btau2 = 0,
           a = 2, b = L, CluDes = "PropMajor", maxIter = 5000, thin = 10, burnin = 5000)
RealRes = MRCUE(F4gammah, F4Gammah, F4se1, F4se2, F4Rblock, rho, coreNum, opt)
```

# Run MR-RAPS, IVW-R, Weighted median and MR-PRESSO

```

library(data.table)
library(ieugwasr)
library(dplyr)
library(MendelianRandomization)
library(mr.raps)
library(MRPRESSO)
#load summarydata
exp <- fread(paste0("./gemma_exp/output/betax.assoc.txt"),head=T)
exp_raw <- exp[,c("chr","rs","beta","se","allele1","allele0","p_wald","n_obs")] 
colnames(exp_raw) <- c("Chr","SNP","b","se","A1","A2","P","N")
exp_data <- exp_raw[which(exp_raw$P<5e-8),]
# LD clump
data <- as.data.frame(exp_data[,c("SNP","P")])
names(data) <- c("rsid","pval")

clump_result <- ld_clump_local(dat=data,clump_kb = 1000,clump_r2 = 0.01,clump_p = 1,
                               bfile = "./all_1000G_EUR_Phase3", 
                               plink_bin = "./PLINK/plink")
names(clump_result) <- c("SNP","P")
exp_ind = exp_data[which(exp_data$SNP %in% clump_result$SNP),]
out <- fread(paste0("./gemma_out/output/betay.assoc.txt"),head=T)
out_raw <- out[,c("chr","rs","beta","se","allele1","allele0","p_wald","n_obs")] 
colnames(out_raw) <- c("Chr","SNP","b","se","A1","A2","P","N")
x_y <- merge(exp_ind,out_raw,by= c("Chr","SNP"))
x_y <- x_y %>% filter((A1.x==A1.y&A2.x==A2.y)|(A1.x==A2.y&A2.x==A1.y))
x_y <- x_y %>% mutate(b.y = ifelse(A1.x==A1.y,b.y,-b.y))
mr_frame <- x_y[,c("SNP","b.x","se.x","b.y","se.y")]
mr_frame <- as.data.frame(mr_frame)
mr_object <- mr_input(bx = mr_frame$b.x, bxse = mr_frame$se.x, 
                      by = mr_frame$b.y, byse = mr_frame$se.y,
                      snps = mr_frame$SNP)
### MR-RAPS
res_raps <- mr.raps.overdispersed.robust(
  mr_frame$b.x, mr_frame$b.y, mr_frame$se.x, mr_frame$se.y, 
  loss.function = "huber", k = 1.345, 
  initialization = c("l2"), suppress.warning = FALSE, 
  diagnosis = FALSE, niter = 20, tol = .Machine$double.eps^0.5)
### IVW-R
res_ivw = mr_ivw(mr_object)
estimate_ivw<-res_ivw$Estimate
se_ivw<- res_ivw$StdError
pvalue_ivw<-  res_ivw$Pvalue					  
ivw <- data.frame(estimate_ivw,se_ivw,pvalue_ivw)
### MR-median
res_median <- mr_median(mr_object, weighting = "weighted", 
                        iterations = 10000, seed = i+10000)
estimate_median <- res_median$Estimate
se_median <- res_median$StdError
pvalue_median <- res_median$Pvalue
median <- data.frame(estimate_median,se_median,pvalue_median)
### MR-PRESSO
mr_presso <- mr_presso(BetaOutcome = "b.y", BetaExposure = "b.x",
                       SdOutcome = "se.y", SdExposure = "se.x", data = mr_frame, 
                       NbDistribution = 5000, SignifThreshold = 0.05)
```

# Run LHC-MR

```
library(lhcMR)
## File paths needed for the analysis
LD.filepath = "./LDscores_filtered.csv" # LD scores
rho.filepath = "./LD_GM2_2prm.csv" # local/SNP-specfic LD scores
ld = "./eur_w_ld_chr/"
hm3 = "./w_hm3.snplist"
## Read in summary statistics file
X = data.table::fread(paste0("./gemma_exp/output/betax.assoc.txt"),head = T)
X = X[,c("rs","allele1","allele0","beta","se","p_wald","n_obs")]
colnames(X) = c("SNP","A1","A2","BETA","SE","P","N")

Y = data.table::fread(paste0("./gemma_out/output/betay.assoc.txt"),head = T)
Y = Y[,c("rs","chr","ps","allele1","allele0","beta","se","p_wald","n_obs")]
colnames(Y) = c("SNP","CHR","BP","A1","A2","BETA","SE","P","N")
## To make sure we have all the columns needed (including allele data), 
##we merge with Neale-provided variants file
XY = merge(X,Y,by="SNP")
exp = XY[,c("SNP","CHR","BP","A1.x","A2.x","BETA.x","SE.x","P.x","N.x")]
out = XY[,c("SNP","CHR","BP","A1.y","A2.y","BETA.y","SE.y","P.y","N.y")]
colnames(exp) = c("SNP","CHR","BP","A1","A2","BETA","SE","P","N")
colnames(out) = c("SNP","CHR","BP","A1","A2","BETA","SE","P","N")
## Step 1
exposure <- paste0("x",i)
outcome <- paste0("y",i)
trait.names=c(exposure,outcome)
input.files = list(exp,out)
df = merge_sumstats(input.files,trait.names,LD.filepath,rho.filepath)

## Step 2
SP_list = calculate_SP(df,trait.names,run_ldsc=FALSE,run_MR=FALSE,hm3=hm3,ld=ld,nStep = 2,
                       SP_single=3,SP_pair=50,SNP_filter=10)
SP_list$sp_mat[,1] <- runif(1,0,0.5)
SP_list$sp_mat[,2] <- runif(1,0,0.5)
save(SP_list,file=sprintf("./%s_%s.Rdata",exposure,outcome))
## Step 3
load(paste0(exposure,"_",outcome,".Rdata"))
res = lhc_mr(SP_list, trait.names, paral_method="lapply", nCores=6, nBlock=200)
```

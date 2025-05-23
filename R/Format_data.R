#' @title Format GWAS summary data for calculating Sample structure.
#' @description Format GWAS summary data for calculating Sample structure.
#' This function is adapted from the format_data() function in MRCIEU/TwoSampleMR.
#'
#' @md
#' @param dat  data.frame with at least SNP A1 A2 signed statistics for calculating z scores and sample size.
#' @param snps.merge data.frame with SNPs to extract. The data frame must have headers: SNP A1 and A2. For example, the hapmap3 SNP list.
#' @param snps.remove a set of SNPs that needed to be removed. For example, SNPs in MHC region.
#' @param snp_col name of a column with rs numbers. The default is `SNP`.
#' @param b_col   name of a column with effect sizes. The default is `b`.
#' @param se_col name of a column with standard errors. The default is `se`.
#' @param freq_col name of a column with effect allele frequencies. The default is `freq`.
#' @param A1_col  name of a column with effect alleles. The default is `A1`.
#' @param A2_col  name of a column with the other alleles. The default is `A2`.
#' @param p_col   name of a column with p-values. The default is `p`.
#' @param n_col name of a column with sample sizes. The default is `n`.
#' @param info_col name of a column with imputation INFO. The default is `Info`.
#' @param min_freq SNPs with allele frequency less than min_freq will be removed. The default is `0.05`.
#' @param n  sample size. If the column for sample size is not available, users can use argument "n" to specify the total sample size for each SNP.
#' @param chi2_max SNPs with tested chi^2 statistics large than chi2_max will be removed. The default is `80`
#'
#' @export
#' @return a data frame with headers: SNP: rs number; A1: effect allele; A2: the other allele;
#' Z: Z score;  N: sample size;  chi2: chi-square statistics;  P: p-value.
#' @importFrom stats pnorm
#' @importFrom dplyr mutate_if
#' @importFrom magrittr %<>%

format_data <- function(dat,
                        snps.merge = NULL,
                        snps.remove = NULL,
                        snp_col="SNP",
                        b_col= "b",
                        se_col= "se",
                        freq_col= "freq",
                        A1_col="A1",
                        A2_col="A2",
                        p_col="p",
                        n_col="n",
                        info_col="INFO",
                        chi2_max = NULL,
                        min_freq=0.05){
  
  message("Begin formatting .... ")
  message("The raw data set has ", nrow(dat), " dat lines")
  cols = c(snp_col, b_col, se_col, freq_col, A1_col, A2_col,p_col, n_col,info_col)
  dat <- as.data.frame(dat)
  dat = dat[, names(dat) %in% cols]
  
  
  if(! snp_col %in% names(dat)){
    stop("SNP column not found")
  }
  
  if(! A1_col %in% names(dat) | (!A2_col %in% names(dat))){
    stop("A1/A2 columns not found")
  }
  
  if(!b_col %in% names(dat)){
    stop("signed statistics not found")
  }
  
  if(!n_col %in% names(dat)){
    stop("Information for sample size not found")
  }
  
  
  ## Remove NA
  names(dat)[which(names(dat) == snp_col)[1]] <- "SNP"
  dat$SNP <- tolower(dat$SNP)
  dat$SNP <- gsub("[[:space:]]", "", dat$SNP)
  dat <- subset(dat, !is.na(SNP))
  
  ## check info
  if(info_col %in% names(dat)){
    names(dat)[which(names(dat) == info_col)[1]] <- "info"
    dat[is.na(dat$info),"info"] = 1
    dat <- subset(dat, info > 0.9)
    message("Remove SNPs with imputation info less than 0.9 ...", ", remaining ", nrow(dat), " SNPs.")
  }
  
  ## Check effect_allele (A1)
  if(A1_col %in% names(dat)){
    
    names(dat)[which(names(dat) == A1_col)[1]] <- "A1"
    
    if(is.logical(dat$A1)){
      dat$A1 <- substr(as.character(dat$A1), 1, 1)
    }
    
    if(!is.character(dat$A1)){
      message("effect allele column is not character data. Coercing...")
      dat$A1 <- as.character(dat$A1)
    }
    
    dat$A1 <- toupper(dat$A1)
    
    index = !grepl("^[ACTG]+$", dat$A1)
    if(any(index)){
      dat = dat[-index,]
      message("effect allele column has some values that are not A/C/T/G. Remove these SNPs...", ", remaining ", nrow(dat), " SNPs.")
      
    }
    
  }
  
  
  ## Check other_allele (A2)
  if(A2_col %in% names(dat)){
    
    names(dat)[which(names(dat) == A2_col)[1]] <- "A2"
    
    if(is.logical(dat$A2))
    {
      dat$A2 <- substr(as.character(dat$A2), 1, 1)
    }
    if(!is.character(dat$A2))
    {
      message("the other allele column is not character data. Coercing...")
      dat$A2 <- as.character(dat$A2)
    }
    
    dat$A2 <- toupper(dat$A2)
    
    index = !grepl("^[ACTG]+$", dat$A2)
    if(any(index)){
      dat = dat[-index,]
      message("the other allele column has some values that are not A/C/T/G. Remove these SNPs...", ", remaining ", nrow(dat), " SNPs.")
    }
  }
  
  index = which((dat$A1=="A" & dat$A2=="T") |
                  (dat$A1=="T" & dat$A2=="A") |
                  (dat$A1=="C" & dat$A2=="G") |
                  (dat$A1=="G" & dat$A2=="C") |
                  (dat$A1=="A" & dat$A2=="A") |
                  (dat$A1=="T" & dat$A2=="T") |
                  (dat$A1=="C" & dat$A2=="C") |
                  (dat$A1=="G" & dat$A2=="G"))
  if(any(index)){
    dat = dat[-index,]
    message("Remove ambiguous SNPs ...", ", remaining ", nrow(dat), " SNPs.")
    rm(index)
  }
  
  
  ## remove MHC SNPs
  if(!is.null(snps.remove)){
    dat <- subset(dat, !SNP %in% snps.remove)
    message("Remove SNPs in MHC region ...", ", remaining ", nrow(dat), " SNPs.")
  }
  
  
  ## Duplicated SNPs
  dup_SNPs = dat$SNP[which(duplicated(dat$SNP))]
  dat = subset(dat, !SNP %in% dup_SNPs)
  message("Remove duplicated SNPs ...", ", remaining ", nrow(dat), " SNPs.")
  
  
  ## merge with hapmap3 SNPlists
  if(!is.null(snps.merge)){
    snps.merge = snps.merge[, c("SNP","A1","A2")]
    colnames(snps.merge) = c("SNP","ref.A1", "ref.A2")
    dat <- merge(dat, snps.merge, by="SNP")
    message("Merge SNPs with the hapmap3 snplist ...", ", remaining ", nrow(dat), " SNPs.")
    
    index = which(!((dat$A1 == dat$ref.A2 & dat$A2 == dat$ref.A1) |
                      (dat$A1 == dat$ref.A1 & dat$A2 == dat$ref.A2) |
                      (dat$A1 == comple(dat$ref.A2) & dat$A2 == comple(dat$ref.A1))|
                      (dat$A1 == comple(dat$ref.A1) & dat$A2 == comple(dat$ref.A2))))
    
    if(any(index)){
      dat = dat[-index,]
      message("Remove SNPs with alleles not matched with the hapmap3 snplist", ", remaining ", nrow(dat), " SNPs.")
      rm(index)
    }
  }
  
  
  ## Check effect size estimate (b) and se
  if(b_col %in% names(dat)){
    names(dat)[which(names(dat) == b_col)[1]] <- "b"
    if(!is.numeric(dat$b)){
      message("b column is not numeric. Coercing...")
      dat$b <- as.numeric(as.character(dat$b))
    }
    dat = dat[is.finite(dat$b),]
    # Check se
    if(se_col %in% names(dat)){
      names(dat)[which(names(dat) == se_col)[1]] <- "se"
      if(!is.numeric(dat$se)){
        message("se column is not numeric. Coercing...")
        dat$se <- as.numeric(as.character(dat$se))
      }
      dat = dat[is.finite(dat$se) & dat$se > 0,]
    }
  }
  
  
  ## Check freq
  if(freq_col %in% names(dat)){
    names(dat)[which(names(dat) == freq_col)[1]] <- "freq"
    if(!is.numeric(dat$freq))
    {
      message("freq column is not numeric. Coercing...")
      dat$freq <- as.numeric(as.character(dat$freq))
    }
    # remove SNP with allele freqcy less than min_freq
    dat[is.na(dat$freq),"freq"] = 0.5
    dat = dat[dat$freq > min_freq & dat$freq < (1-min_freq), ]
    message("Remove SNPs with MAF less than", min_freq, ", remaining ", nrow(dat), " SNPs.")
  }
  
  
  ## Check n
  if(n_col %in% names(dat)){
    names(dat)[which(names(dat) == n_col)[1]] <- "n"
    
    if(!is.numeric(dat$n)){
      message(samplesize_col, " column is not numeric")
      dat$n <- as.numeric(as.character(dat$n))
    }

  }
  
  ## Check pval
  if(p_col %in% names(dat)){
    names(dat)[which(names(dat) == p_col)[1]] <- "p"
    if(!is.numeric(dat$p)){
      message("p-value column is not numeric. Coercing...")
      dat$p <- as.numeric(dat$p)
    }
    dat = subset(dat, p>=0 & p <=1)
    message("Remove SNPs with p-value < 0 or p-value > 1", ", remaining ", nrow(dat), " SNPs.")
  }
  
  # calculate z if not contain z, but with b se
  if((! "z" %in% names(dat))  & ("b" %in% names(dat) & "se" %in% names(dat))){
    message("Inferring z score from b/se ...")
    dat$z = dat$b/dat$se
    dat$chi2 = dat$z^2
  }
  
  # calculate chi2 if no chi2
  if(!"chi2" %in% names(dat)) dat$chi2 = dat$z^2
  
  # calculate p if not contain p
  if(!"p" %in% names(dat)) dat$p = pchisq(dat$chi2, 1, lower.tail = F)
  
  if(! "z" %in% names(dat)){
    stop("Error: No information for z score ")
  }
  
  # check missing
  dat = dat[, c("SNP","A1","A2","z","n","chi2","p")]
  dat = na.omit(dat)
  message("Remove missing values", ", remaining ", nrow(dat), " SNPs.")
  
  n_min = mean(dat$n) - 5* sd(dat$n)
  n_max = mean(dat$n) + 5* sd(dat$n)
  dat = subset(dat, n >= n_min  & n <= n_max)
  message("Remove SNPs with sample size 5 standard deviations away from the mean", ", remaining ", nrow(dat), " SNPs.")
  
  chi2_max = max(c(80, median(dat$n)/1000))
  dat = subset(dat, chi2 < chi2_max)
  message("Remove SNPs with chi2 > chi2_max ... ", ", remaining ", nrow(dat), " SNPs.")
  
 # message("The formatted data has ", nrow(dat), " dat lines. \n")
  
  
  colnames(dat) = c("SNP","A1","A2","Z","N","chi2","P")
  
  dat %<>% dplyr::mutate_if(is.integer, as.numeric)
  
  return(dat)
  
}



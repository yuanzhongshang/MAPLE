#' @title A function estimate parameters of sample structure by LD score regression.
#' @description Format GWAS summary data for calculating Sample structure.
#' This function is adapted from the est_paras() function in YangLabHKUST/MR-APSS.
#' @param dat1 formatted GWAS summary-level data for exposure.
#' @param dat2 formatted GWAS summary-level data for outcome.
#' @param trait1.name specify the name of exposure, default `exposure`.
#' @param trait2.name specify the name of outcome, default `outcome`.
#' @param LDSC logical, whether to run LD score regression, default `TRUE`. If `FALSE`, the function will not give the parameter estimates.
#' @param h2.fix.intercept logical, whether to fix LD score regression intercept to 1, default `FALSE`.
#' @param ldscore.dir specify the path to the LD score files.
#' @param ld specify the data.frame for the ld scores of SNPs if available and setting ldscore.dir =NULL.
#' @param M specify the number of SNPs used for calculating LD score.
#'
#' @return List with the following elements:
#' \describe{
#' \item{Omega}{the estimated Omega matrix capturing the effects of sample structure}
#' \item{Omega.se}{the estimated matrix consists of the standard errors of the intercept estimates obtained from LD score regression}
#' }
#'
#' @export
#' @importFrom readr read_delim
#' @importFrom readr read_csv


est_SS <- function(dat1,
                      dat2,
                      trait1.name = "exposure",
                      trait2.name = "outcome",
                      LDSC = T,
                      h2.fix.intercept = F,
                      ldscore.dir = NULL,
                      ld=NULL,
                      M=NULL){
  dat1 = format_data(dat1,
                     snp_col = "SNP",
                     b_col = "b",
                     se_col = "se",
                     freq_col = "frq_A1",
                     A1_col = "A1",
                     A2_col = "A2",
                     p_col = "P",
                     n_col = "N")
  dat2 = format_data(dat2,
                     snp_col = "SNP",
                     b_col = "b",
                     se_col = "se",
                     freq_col = "frq_A1",
                     A1_col = "A1",
                     A2_col = "A2",
                     p_col = "P",
                     n_col = "N")
  dat1 %<>% dplyr::mutate_if(is.integer, as.numeric)
  dat2 %<>% dplyr::mutate_if(is.integer, as.numeric)
  dat1 %<>% dplyr::mutate_if(is.factor, as.character)
  dat2 %<>% dplyr::mutate_if(is.factor, as.character)
  
  #message("Merge dat1 and dat2 by SNP ...")
  dat = merge(dat1, dat2, by="SNP")
  
  
  flip.index = which((dat$A1.x == dat$A2.y & dat$A1.y == dat$A2.x) |
                       (dat$A1.x ==comple(dat$A2.y) & dat$A1.y == comple(dat$A2.x)))
  
  dat[,"A1.y"] = dat[,"A1.x"]
  dat[,"A2.y"] = dat[,"A2.x"]
  dat[flip.index ,"Z.y"] = -dat[flip.index ,"Z.y"]
  
  
  #message("Read in LD scores ... ")
  if(is.null(ldscore.dir)&is.null(ld)) stop("Please provide the information on LD scores")
  
  if(is.null(ld) & !is.null(ldscore.dir)){
    
    ld <- suppressMessages(readr::read_delim(paste0(ldscore.dir,"/1.l2.ldscore.gz"), "\t", escape_double = FALSE, trim_ws = TRUE,progress = F))
    
    for(i in 2:22){
      ld <- rbind(ld,suppressMessages(readr::read_delim(paste0(ldscore.dir, "/", i,".l2.ldscore.gz"), "\t", escape_double = FALSE, trim_ws = TRUE,progress = F)))
    }
  }
  
  if(is.null(M) & !is.null(ldscore.dir)){
    m.chr  <- suppressMessages(readr::read_csv(paste0(ldscore.dir,"/1.l2.M_5_50"),  col_names = FALSE))
    
    for(i in 2:22){
      m.chr <- rbind(m.chr,suppressMessages(readr::read_csv(paste0(ldscore.dir, "/", i,".l2.M_5_50"),  col_names = FALSE)))
    }
    
    M = sum(m.chr)  # the number of SNPs include in the LD score estimation
  }
  
  merged  = merge(dat, ld, by="SNP")
  C = NULL
  gcres12  = NULL
  
  if(LDSC){
    #message("Begin estimation of Omega using LDSC ...")
    gcres12 = ldsc_SS(merged,
                      trait1.name = "exposure",
                      trait2.name = "outcome",
                      Twostep = T,
                      M=M,
                      h2.fix.intercept = h2.fix.intercept)
    
    Omega = matrix(as.vector(gcres12$I), nrow=2, ncol=2)
    Omega.se = matrix(as.vector(gcres12$I.se), nrow=2, ncol=2)
  }
  
  return(list(ldsc_res=gcres12,
              Omega = Omega, Omega.se=Omega.se))
}

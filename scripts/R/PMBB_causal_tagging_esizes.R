################################################################################
######################## Causal-Tagging Simulations ############################
######################## Phenotype-specific simulations ########################
################################################################################
#' 
#' This script computes analytical parameters for the average causal effect under
#' the global model, given simulation parameters rho and r2. Numerical estimates 
#' are also computed, by simulating multiple realizations of beta_j's before
#' computing local model and global model hat(LAACor_j)'s.
#' 

#### Define Libraries and Directories ------------------------------------------
working_dir <- "/project/mathilab/aaw/admix_prs/"
eur_geno_backup_dir <-
  "/project/mathilab/aaw/admix_prs/PMBB/Synthetic_Data/Eur_Genomes/Causal_Tag_Data/backup/"
adm_geno_dir <- 
  "/project/mathilab/aaw/admix_prs/PMBB/v2/Causal_Tag_Data/"
result_dir <- "/project/mathilab/aaw/admix_prs/results/042825/"

library(dplyr)
library(bigsnpr)

# Read argument to select relevant combination row
args = commandArgs(trailingOnly=TRUE)

# Read argument to select relevant phenotype and R2_VEC
s <- as.numeric(args[1])

PHENOTYPE_VEC <- c("Standing-Height","Weight","Body-Mass-Index",
                   "Triglycerides","Neutrophil_Count","Platelet-counts")
PHENOTYPE <- PHENOTYPE_VEC[s]
RHO_VEC <- c(0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,
             0.7,0.75,0.8,0.85,0.9,
             0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1)

#### Helper functions ----------------------------------------------------------
## Concatenate digits for saving file easily ===================================
concatenate_digits <- function(decimal) {
  # Convert the decimal to a character string
  decimal_str <- as.character(decimal)
  
  # Remove the period and concatenate the digits
  result <- paste0(gsub("[^0-9]", "", decimal_str), collapse = "")
  
  return(result)
}

## Pool ADM cohort haplotype and LANC matrices =================================
getPooledHapLANC <- function(seed_haplotypes,
                             seed_afr_ancestries) {
  # Initialize hap and lanc matrices
  big_haps <- NULL
  big_afr_ancs <- NULL
  
  # For each quantile
  for (q in 1:4) {
    # Generate haplotype and African ancestry matrices
    haps <- cbind(seed_haplotypes[[paste0("Q",q,"H1")]], # n x 2p
                  seed_haplotypes[[paste0("Q",q,"H2")]])
    big_haps <- rbind(big_haps,haps)
    
    afr_ancs <- cbind(seed_afr_ancestries[[paste0("Q",q,"H1")]], # n x 2p
                      seed_afr_ancestries[[paste0("Q",q,"H2")]])
    big_afr_ancs <- rbind(big_afr_ancs,afr_ancs)
  }
  
  # Combine tag and causal lists
  to_return <- list(HAP=big_haps,
                    AFR_ANC=big_afr_ancs)
  # Return
  return(to_return)
}

## Check correlations in admixed data ==========================================
checkADMCor <- function(causal_hap, tag_hap) {
  # get no. SNPs
  p <- ncol(causal_hap)
  n_snps <- p/2 
  causal_geno <- causal_hap[,1:n_snps] + causal_hap[,(n_snps+1):(2*n_snps)]
  tag_geno <- tag_hap[,1:n_snps] + tag_hap[,(n_snps+1):(2*n_snps)]
  
  # Compute correlation between each causal/tagging at the genotype level
  to_return <- c()
  for (i in 1:n_snps) {
    to_return <- c(to_return, cor(causal_geno[,i],tag_geno[,i],
                                  use="complete.obs"))
  }
  
  # Return correlation vector
  return(to_return)
}

## Check correlations in unadmixed data ========================================
checkEURCor <- function(causal_geno, tag_geno) {
  # get no. SNPs
  n_snps <- ncol(causal_geno)
  
  # Compute correlation between each causal/tagging at the genotype level
  to_return <- c()
  for (i in 1:n_snps) {
    to_return <- c(to_return, cor(causal_geno[,i],tag_geno[,i],
                                  use="complete.obs"))
  }
  
  # Return correlation vector
  return(to_return)
}

## Get paired and marginal frequencies =========================================
#' Performs computation of allele frequencies for computing beta_t^E and beta_t^A 
#' in the admixed cohort. AFs for causal and tagging variants are computed. 
#' For causal variants, we also compute the local ancestry-agnostic AF 
#' (since we don't need local ancestry).
#' 
#' adm_cohort_list is a list(CAUSAL_LANC=causal_lanc_mat,
#'                           CAUSAL_HAP=causal_hap_mat,
#'                           TAG_LANC=tag_lanc_mat,
#'                           TAG_HAP=tag_hap_mat)
#' See: https://www.arslanzaidi.com/post/effect-size-of-a-tag-variant/
getfcft <- function(adm_cohort_list) {
  # extract the three data matrices
  causal_lanc_mat <- adm_cohort_list$CAUSAL_LANC
  causal_hap_mat <- adm_cohort_list$CAUSAL_HAP
  tag_lanc_mat <- adm_cohort_list$TAG_LANC
  tag_hap_mat <- adm_cohort_list$TAG_HAP
  
  # get the number of variants
  nvars <- ncol(tag_hap_mat)/2
  
  # f_t and f_c can be vectorized
  afr_causal_hap_mat <- causal_hap_mat*causal_lanc_mat # Afr causal
  unfolded_afr_causal_ad <- colSums(afr_causal_hap_mat,na.rm=TRUE) 
  unfolded_afr_causal_lanc_counts <- colSums(causal_lanc_mat,na.rm=TRUE) 
  afr_tag_hap_mat <- tag_hap_mat*tag_lanc_mat # Afr tag
  unfolded_afr_tag_ad <- colSums(afr_tag_hap_mat,na.rm=TRUE) 
  unfolded_afr_tag_lanc_counts <- colSums(tag_lanc_mat,na.rm=TRUE) 
  
  folded_afr_causal_ad <- unfolded_afr_causal_ad[1:nvars]+
    unfolded_afr_causal_ad[(nvars+1):(2*nvars)]
  folded_afr_tag_ad <- unfolded_afr_tag_ad[1:nvars]+
    unfolded_afr_tag_ad[(nvars+1):(2*nvars)]
  folded_afr_causal_lanc_counts <- unfolded_afr_causal_lanc_counts[1:nvars]+
    unfolded_afr_causal_lanc_counts[(nvars+1):(2*nvars)]
  folded_afr_tag_lanc_counts <- unfolded_afr_tag_lanc_counts[1:nvars]+
    unfolded_afr_tag_lanc_counts[(nvars+1):(2*nvars)]
  afr_causal_af <- folded_afr_causal_ad/folded_afr_causal_lanc_counts
  afr_tag_af <- folded_afr_tag_ad/folded_afr_tag_lanc_counts
  
  eur_causal_hap_mat <- causal_hap_mat*(1-causal_lanc_mat) # Eur causal
  unfolded_eur_causal_ad <- colSums(eur_causal_hap_mat,na.rm=TRUE)
  unfolded_eur_causal_lanc_counts <- colSums(1-causal_lanc_mat,na.rm=TRUE) 
  eur_tag_hap_mat <- tag_hap_mat*(1-tag_lanc_mat) # Eur tag
  unfolded_eur_tag_ad <- colSums(eur_tag_hap_mat,na.rm=TRUE) 
  unfolded_eur_tag_lanc_counts <- colSums(1-tag_lanc_mat,na.rm=TRUE) 
  
  folded_eur_causal_ad <- unfolded_eur_causal_ad[1:nvars]+
    unfolded_eur_causal_ad[(nvars+1):(2*nvars)]
  folded_eur_tag_ad <- unfolded_eur_tag_ad[1:nvars]+
    unfolded_eur_tag_ad[(nvars+1):(2*nvars)]
  folded_eur_causal_lanc_counts <- unfolded_eur_causal_lanc_counts[1:nvars]+
    unfolded_eur_causal_lanc_counts[(nvars+1):(2*nvars)]
  folded_eur_tag_lanc_counts <- unfolded_eur_tag_lanc_counts[1:nvars]+
    unfolded_eur_tag_lanc_counts[(nvars+1):(2*nvars)]
  eur_causal_af <- folded_eur_causal_ad/folded_eur_causal_lanc_counts
  eur_tag_af <- folded_eur_tag_ad/folded_eur_tag_lanc_counts
  
  # local ancestry-agnostic f_c
  unfolded_lanc_agnostic_causal_af <- colMeans(causal_hap_mat,na.rm=TRUE) 
  lanc_agnostic_causal_af <- 
    (unfolded_lanc_agnostic_causal_af[1:nvars]+
       unfolded_lanc_agnostic_causal_af[(nvars+1):(2*nvars)])/2
  
  # f_tc requires considering both tag_hap_mat and causal_hap_mat
  # as well as tag_lanc_mat and causal_lanc_mat
  afr_tagcausal_hap_mat <- 
    afr_causal_hap_mat*afr_tag_hap_mat # entry = 1 iff x_t = x_c = 1 
  eur_tagcausal_hap_mat <- 
    eur_causal_hap_mat*eur_tag_hap_mat # entry = 1 iff x_t = x_c = 1 
  afr_tagcausal_lanc_mat <- 
    causal_lanc_mat*tag_lanc_mat # entry = 1 iff a_i^t=a_i^c=Afr
  eur_tagcausal_lanc_mat <- 
    (1-causal_lanc_mat)*(1-tag_lanc_mat) # entry = 1 iff a_i^t=a_i^c=Eur
  
  unfolded_afr_tagcausal_ad <- colSums(afr_tagcausal_hap_mat,na.rm=TRUE) 
  unfolded_eur_tagcausal_ad <- colSums(eur_tagcausal_hap_mat,na.rm=TRUE) 
  unfolded_afr_tagcausal_lanc_counts <- colSums(afr_tagcausal_lanc_mat,na.rm=TRUE)
  unfolded_eur_tagcausal_lanc_counts <- colSums(eur_tagcausal_lanc_mat,na.rm=TRUE)
  
  folded_afr_tagcausal_ad <- unfolded_afr_tagcausal_ad[1:nvars] +
    unfolded_afr_tagcausal_ad[(nvars+1):(2*nvars)]
  folded_eur_tagcausal_ad <- unfolded_eur_tagcausal_ad[1:nvars] +
    unfolded_eur_tagcausal_ad[(nvars+1):(2*nvars)]
  folded_afr_tagcausal_lanc_counts <- unfolded_afr_tagcausal_lanc_counts[1:nvars]+
    unfolded_afr_tagcausal_lanc_counts[(nvars+1):(2*nvars)]
  folded_eur_tagcausal_lanc_counts <- unfolded_eur_tagcausal_lanc_counts[1:nvars]+
    unfolded_eur_tagcausal_lanc_counts[(nvars+1):(2*nvars)]
  
  afr_tagcausal_af <- folded_afr_tagcausal_ad/folded_afr_tagcausal_lanc_counts
  eur_tagcausal_af <- folded_eur_tagcausal_ad/folded_eur_tagcausal_lanc_counts
  
  eur_r_xc_xt <- (eur_tagcausal_af - eur_causal_af*eur_tag_af)/
    sqrt(eur_causal_af*(1-eur_causal_af)*eur_tag_af*(1-eur_tag_af))
  afr_r_xc_xt <- (afr_tagcausal_af - afr_causal_af*afr_tag_af)/
    sqrt(afr_causal_af*(1-afr_causal_af)*afr_tag_af*(1-afr_tag_af))
  
  # Return
  return(list(fc_lanc_agnostic=lanc_agnostic_causal_af,
              ftE = eur_tag_af,
              ftA = afr_tag_af,
              fcE = eur_causal_af,
              fcA = afr_causal_af,
              fctE = eur_tagcausal_af,
              fctA = afr_tagcausal_af,
              rctE = eur_r_xc_xt,
              rctA = afr_r_xc_xt))
}

## Prepare EUR cohort data =====================================================
#' Take in two obj.bigSNPs (one causal and one tagging) and two index vectors 
#' that re-order the columns (for ensuring causal and tagging are paired right).
#' 
#' Designed for genotypes, this function relies on cor(X_t,X_c) rather than 
#' computing the Wright-Robertson-Hill formula,
#' 
#' r = (f_AB - f_A*f_B)/sqrt((1-f_A)(1-f_B)f_Af_B)
#' 
#' There should be a small difference between the two quantities. It'd be ideal 
#' to run a HMM to phase and then compute r above, but this takes extra time so
#' we avoid it.
#' 
#' Computes the following quantities that are important for R2_PGS calculation:
#' - MAF vectors for tagging and causal variants (fc and ft)
#' - approximate LD between tagging and causal variants (rct)
#' - demeaned versions of the genotype matrices
#' - denominators to normalize PGS 
getunadmixedfcft <- function(bigsnpr_G_causal,bigsnpr_G_tag,
                             causal_index, tag_index) {
  # Get number of variants and number of individuals
  n_ind <- nrow(bigsnpr_G_causal)
  n_var <- ncol(bigsnpr_G_causal)
  
  # Get MAF (fc and ft)
  # [!] Variants are arranged in the same order as in ADM cohort 
  causal_maf <- colMeans(bigsnpr_G_causal[1:n_ind,causal_index],na.rm=T)/2
  tag_maf <- colMeans(bigsnpr_G_tag[1:n_ind,tag_index],na.rm=T)/2
  
  # Get demeaned genotypes
  # [!] Variants are arranged in the same order as in ADM cohort 
  copy_bigsnpr_G_causal <- big_copy(bigsnpr_G_causal,
                                    ind.col=causal_index, # rearrange variant index
                                    type="double") # convert to double for big_increment
  neg_colmeans <- t(replicate(nrow(bigsnpr_G_causal), 0-2*causal_maf))
  big_increment(copy_bigsnpr_G_causal,neg_colmeans)
  
  copy_bigsnpr_G_tag <- big_copy(bigsnpr_G_tag,
                                 ind.col=tag_index, # rearrange variant index
                                 type="double") # convert to double for big_increment
  neg_colmeans <- t(replicate(nrow(bigsnpr_G_tag), 0-2*tag_maf))
  big_increment(copy_bigsnpr_G_tag,neg_colmeans)
  
  # [!] Set NAs to 0 [!]
  big_apply(copy_bigsnpr_G_causal, function(X,ind) {
    X.sub <- X[,ind,drop=FALSE]
    ind_na <- which(is.na(X.sub), arr.ind=TRUE)
    ind_na[, 2] <- ind[ind_na[, 2]]
    X[ind_na] <- 0
    NULL
  }, a.combine='c',block.size=200)
  
  big_apply(copy_bigsnpr_G_tag, function(X,ind) {
    X.sub <- X[,ind,drop=FALSE]
    ind_na <- which(is.na(X.sub), arr.ind=TRUE)
    ind_na[, 2] <- ind[ind_na[, 2]]
    X[ind_na] <- 0
    NULL
  }, a.combine='c',block.size=200)
  
  
  # Compute numerator: Sum of element-wise product
  num <- colSums(copy_bigsnpr_G_causal[]*copy_bigsnpr_G_tag[],na.rm=TRUE)
  
  # Compute denominators: Square root of sum of squared deviations
  denom_causal <- sqrt(colSums(copy_bigsnpr_G_causal[]^2,na.rm=TRUE))
  denom_tag <- sqrt(colSums(copy_bigsnpr_G_tag[]^2,na.rm=TRUE))
  
  # Compute correlation
  r_xc_xt <- num / (denom_causal * denom_tag)
  
  # Return MAF, square-root of variances of MAFs and demeaned genotypes
  return(list(fc = causal_maf,
              ft = tag_maf,
              rct = r_xc_xt,
              G_demean_causal = copy_bigsnpr_G_causal,
              G_demean_tag = copy_bigsnpr_G_tag,
              causal_beta_denom = sqrt(sum(causal_maf*(1-causal_maf))),
              tag_beta_denom = sqrt(sum(tag_maf*(1-tag_maf)))))
}

## Simulate causal and tagging effects =========================================
#' Take in outputs of getfcft and getunadmixedfcft, and uses formula for tagging 
#' variant effect (as function of causal variant effect and LD and AF) to simulate 
#' both beta_causal and beta_tagging^Afr and beta_tagging^Eur in ADM. 
#' 
#' This function simulates tagging effects using AF and LD estimated from 
#' unadmixed cohort, as well as Eur-LANC-conditioned AF and LD estimated from 
#' admixed cohort. The former should be more accurate and reflective of real 
#' LD pattern differences, since there is a very small effective sample size for 
#' estimating LD and AF for Eur-LANC in the admixed cohort.
#' 
#' Causal effects are assumed identical across ancestries, but LD and AF differences
#' lead to differences in tagging effects. 
getEffectVecs <- function(r2, fcft_list,
                          cross_pop_rho=1,
                          unadmixed_fcft_list) {
  # Half the r2 so that the genotypes will have desired squared correlation
  eff_r2 <- r2/2
  
  # Get number of variants
  nvars <- length(fcft_list$ftE)
  
  # Generate (unscaled) causal effects
  causal_beta_vec <- rnorm(n=nvars,mean=0,sd=sqrt(eff_r2))
  eur_causal_beta_vec <- rnorm(n=nvars,mean=0,sd=sqrt(eff_r2))
  ind_gaussian <- rnorm(n=nvars,mean=0,sd=1)
  afr_causal_beta_vec <- eur_causal_beta_vec*cross_pop_rho + 
    ind_gaussian*sqrt(eff_r2*(1-cross_pop_rho^2))
  
  # Generate tagging effects for each ancestry
  admixed_eur_tag_beta_vec <- eur_causal_beta_vec*fcft_list$rctE*
    sqrt((fcft_list$fcE*(1-fcft_list$fcE))/(fcft_list$ftE*(1-fcft_list$ftE)))
  eur_tag_beta_vec <- eur_causal_beta_vec*unadmixed_fcft_list$rct* # difference in Eur AF estimated from external cohort and ADM in-sample!
    sqrt((unadmixed_fcft_list$fc*(1-unadmixed_fcft_list$fc))/
           (unadmixed_fcft_list$ft*(1-unadmixed_fcft_list$ft)))
  afr_tag_beta_vec <- afr_causal_beta_vec*fcft_list$rctA*
    sqrt((fcft_list$fcA*(1-fcft_list$fcA))/(fcft_list$ftA*(1-fcft_list$ftA)))
  
  # Set NAs to zero; NAs arise from non-segregating loci 
  afr_tag_beta_vec[which(is.na(afr_tag_beta_vec))] <- 0
  eur_tag_beta_vec[which(is.na(eur_tag_beta_vec))] <- 0
  admixed_eur_tag_beta_vec[which(is.na(admixed_eur_tag_beta_vec))] <- 0
  
  # Return
  return(list(EUR_CAUSAL=eur_causal_beta_vec,
              AFR_CAUSAL=afr_causal_beta_vec,
              EUR_TAG_ADMIXED=admixed_eur_tag_beta_vec,
              EUR_TAG=eur_tag_beta_vec,
              AFR_TAG=afr_tag_beta_vec))
}

## getOmegas ===================================================================
#'
#' Calculates omega_primes and omegas that appear in the formulae for local ancestry
#' average causal effect. NOTE: We use tagging variant matrix to approximate the 
#' global ancestry of the individual.  
#' 
getOmegas <- function(adm_cohort_list) {
  # Extract fixed data objects
  afr_A_prime <- adm_cohort_list$CAUSAL_LANC # n x 2p
  afr_A <- adm_cohort_list$TAG_LANC # n x 2p
  eur_A <- 1-afr_A
  eur_A_prime <- 1-afr_A_prime
  causal_hap_mat <- adm_cohort_list$CAUSAL_HAP # n x 2p
  tag_hap_mat <- adm_cohort_list$TAG_HAP # n x 2p 
  n_vars <- ncol(afr_A_prime)/2
  n_inds <- nrow(afr_A_prime)
  # Get haplotype-specific causal and tagging ancestry matrices 
  afr_A_prime_1 <- afr_A_prime[,1:n_vars]
  afr_A_prime_2 <- afr_A_prime[,(n_vars+1):(2*n_vars)]
  afr_A_1 <- afr_A[,1:n_vars]
  afr_A_2 <- afr_A[,(n_vars+1):(2*n_vars)]
  
  # Get global Afr ancestry matrix
  ind_global_anc_vec <- rowMeans(afr_A,na.rm=T) # global Afr ancestry based on causal loci 
  ind_global_anc_prime_vec <- rowMeans(afr_A_prime,na.rm=T) # global Afr ancestry based on causal loci 
  message(date(), 
          ": Correlation of global ancestries computed using causal vs tagging variants = ",
          cor(ind_global_anc_vec,ind_global_anc_prime_vec))
  message(date(), 
          ": Mean global ancestry computed using causal variants = ",
          mean(ind_global_anc_prime_vec))
  message(date(), 
          ": Mean global ancestry computed using tagging variants = ",
          mean(ind_global_anc_vec))
  
  ind_global_anc_mat <- matrix(rep(ind_global_anc_vec, times = n_vars), # n x p
                               nrow = length(ind_global_anc_vec), ncol = n_vars)
  
  # Calculate hat_a_j and hat_a_j_prime
  hat_a_j_prime_1 <- colMeans(afr_A_prime_1,na.rm=T) # p x 1
  hat_a_j_prime_2 <- colMeans(afr_A_prime_2,na.rm=T) # p x 1
  hat_a_j_1 <- colMeans(afr_A_1,na.rm=T) # p x 1
  hat_a_j_2 <- colMeans(afr_A_2,na.rm=T) # p x 1
  
  # Compute numerator and denominator of omega terms
  omega1_prime_num <- colSums((afr_A_prime_1+afr_A_prime_2)*(1-ind_global_anc_mat),
                              na.rm=T)
  omega2_prime_num <- colSums((afr_A_prime_1+afr_A_prime_2)*ind_global_anc_mat,
                              na.rm=T)
  omega3_prime_num <- colSums((2-afr_A_prime_1-afr_A_prime_2)*(1-ind_global_anc_mat),
                              na.rm=T)
  omega4_prime_num <- colSums((2-afr_A_prime_1-afr_A_prime_2)*ind_global_anc_mat,
                              na.rm=T)
  omega12_prime_denom <- n_inds*(hat_a_j_prime_1+hat_a_j_prime_2)
  omega34_prime_denom <- n_inds*(2-hat_a_j_prime_1-hat_a_j_prime_2)
  
  omega1_num <- colSums((afr_A_1+afr_A_2)*(1-ind_global_anc_mat),
                        na.rm=T)
  omega2_num <- colSums((afr_A_1+afr_A_2)*ind_global_anc_mat, 
                        na.rm=T)
  omega3_num <- colSums((2-afr_A_1-afr_A_2)*(1-ind_global_anc_mat),
                        na.rm=T)
  omega4_num <- colSums((2-afr_A_1-afr_A_2)*ind_global_anc_mat,
                        na.rm=T)
  omega12_denom <- n_inds*(hat_a_j_1+hat_a_j_2)
  omega34_denom <- n_inds*(2-hat_a_j_1-hat_a_j_2)
  
  # Return 
  to_return <- data.frame(OMEGA1_PRIME=omega1_prime_num/omega12_prime_denom,
                          OMEGA2_PRIME=omega2_prime_num/omega12_prime_denom,
                          OMEGA3_PRIME=omega3_prime_num/omega34_prime_denom,
                          OMEGA4_PRIME=omega4_prime_num/omega34_prime_denom,
                          OMEGA1=omega1_num/omega12_denom,
                          OMEGA2=omega2_num/omega12_denom,
                          OMEGA3=omega3_num/omega34_denom,
                          OMEGA4=omega4_num/omega34_denom)
}


#### Preamble ------------------------------------------------------------------
#source(paste0(working_dir,"PMBB_sim_master_script_geno_AFs_vers_030425.R"))
message(date(), ": Simulating and calculating genetic effects with ", PHENOTYPE,
        " haplotype/local ancestry matrices under local and global models")

#### Main Body -----------------------------------------------------------------
## 1. Load EUR and ADM data ====================================================
# European cohort causal and tagging matrices
obj.bigSNP <- snp_attach(paste0(eur_geno_backup_dir,
                                PHENOTYPE,"_causal.rds"))
pheno_EurAm_G <- obj.bigSNP$genotypes

obj.bigSNP.tag <- snp_attach(paste0(eur_geno_backup_dir,
                                    PHENOTYPE,"_tag.rds"))
pheno_EurAm_G_tag <- obj.bigSNP.tag$genotypes

# ADM cohort causal and tagging haplotypes / LANC matrices 
pheno_ADM_data <- readRDS(paste0(adm_geno_dir,PHENOTYPE,"_causal_tag_data.rds"))
causal_haps <- pheno_ADM_data[["CAUSAL_HAPLOTYPES"]]
causal_afr_ancs <- pheno_ADM_data[["CAUSAL_AFR_ANCESTRY"]]
tag_haps <- pheno_ADM_data[["TAG_HAPLOTYPES"]]
tag_afr_ancs <- pheno_ADM_data[["TAG_AFR_ANCESTRY"]]
var_ids_df <- pheno_ADM_data[["VAR_IDS"]]
rm(pheno_ADM_data)

causal_pooled_list <- getPooledHapLANC(seed_haplotypes=causal_haps,
                                       seed_afr_ancestries=causal_afr_ancs)
tag_pooled_list <- getPooledHapLANC(seed_haplotypes=tag_haps,
                                    seed_afr_ancestries=tag_afr_ancs)

# Reorder indices of EUR to match the ordering of variants in ADM haplotypes
# These can be used to rearrange the EUR cohort columns and also the allele freqs
causal_indices <- match(var_ids_df$CAUSAL,table=obj.bigSNP$map$marker.ID)
tag_indices <- match(var_ids_df$TAG,table=obj.bigSNP.tag$map$marker.ID)

# [!] Statistic Check: LD of tagging and causal variant
test_cor_vec <- checkADMCor(causal_hap=causal_pooled_list[["HAP"]],
                            tag_hap=tag_pooled_list[["HAP"]])
ADM_MEAN_SQ_LD <- mean(test_cor_vec^2,na.rm=T)
message("Mean squared LD / emp corr^2 between Tagging and Causal in ADM = ", 
        sprintf("%.4f", ADM_MEAN_SQ_LD))

# Precompute LD and allele frequencies (for use throughout simulations)
emp_af_list <- getfcft(adm_cohort_list=list(CAUSAL_LANC=causal_pooled_list[["AFR_ANC"]],
                                            CAUSAL_HAP=causal_pooled_list[["HAP"]],
                                            TAG_LANC=tag_pooled_list[["AFR_ANC"]],
                                            TAG_HAP=tag_pooled_list[["HAP"]]))
message("Mean squared LD / emp corr^2 between Tagging and Causal in AFR (using ADM Afr-LANC) = ", 
        sprintf("%.4f",mean(emp_af_list$rctA^2,na.rm=T)))
message("Mean squared LD / emp corr^2 between Tagging and Causal in EUR (using ADM Eur-LANC) = ", 
        sprintf("%.4f",mean(emp_af_list$rctE^2,na.rm=T)))

unadmixed_emp_af_list <- getunadmixedfcft(bigsnpr_G_causal=pheno_EurAm_G,
                                          bigsnpr_G_tag=pheno_EurAm_G_tag,
                                          causal_index=causal_indices, 
                                          tag_index=tag_indices)

eur_test_cor_vec <- checkEURCor(causal_geno=unadmixed_emp_af_list$G_demean_causal,
                                tag_geno=unadmixed_emp_af_list$G_demean_tag)
EUR_MEAN_SQ_LD <- mean(eur_test_cor_vec^2,na.rm=T)
message("Mean squared LD / emp corr^2 between Tagging and Causal in EUR (using Eur cohort) = ", 
        sprintf("%.4f",EUR_MEAN_SQ_LD))
message("Pearson corr between EUR causal allele frequencies estimated in-ADM-sample and in Eur cohort = ",
        sprintf("%.4f",cor(unadmixed_emp_af_list$fc,emp_af_list$fcE,use="complete.obs")))
message("Pearson corr between EUR tagging allele frequencies estimated in-ADM-sample and in Eur cohort = ",
        sprintf("%.4f",cor(unadmixed_emp_af_list$ft,emp_af_list$ftE,use="complete.obs")))

message(">>> Allele freq comparisons in EUR: Average f_c-f_t = ", 
        sprintf("%.4f",mean(unadmixed_emp_af_list$fc-unadmixed_emp_af_list$ft)))
message(">>> Allele freq comparisons in EUR: Var(f_c-f_t) = ", 
        sprintf("%.4f",var(unadmixed_emp_af_list$fc-unadmixed_emp_af_list$ft)))
message(">>> Allele freq comparisons in EUR: Average |f_c-f_t| = ", 
        sprintf("%.4f",mean(abs(unadmixed_emp_af_list$fc-unadmixed_emp_af_list$ft))))

# [!] Statistic check: global African ancestry for tagging and causal haplotypes
ind_tag_glo_AFR_anc <- rowMeans(tag_pooled_list[["AFR_ANC"]],na.rm=T)
ind_causal_glo_AFR_anc <- rowMeans(causal_pooled_list[["AFR_ANC"]],na.rm=T)
a_bar_tag <- mean(ind_tag_glo_AFR_anc)
a_bar_causal <- mean(ind_causal_glo_AFR_anc)
n_inds <- length(ind_tag_glo_AFR_anc)
n_vars <- ncol(tag_pooled_list[["AFR_ANC"]])/2

# Precompute omega and omega_prime vectors (for use throughout simulations)
omega_list <- getOmegas(adm_cohort_list=list(CAUSAL_LANC=causal_pooled_list[["AFR_ANC"]],
                                             CAUSAL_HAP=causal_pooled_list[["HAP"]],
                                             TAG_LANC=tag_pooled_list[["AFR_ANC"]],
                                             TAG_HAP=tag_pooled_list[["HAP"]]))

# Get haplotype-specific causal and tagging local ancestry matrices 
afr_A_prime_1 <- causal_pooled_list[["AFR_ANC"]][,1:n_vars]
afr_A_prime_2 <- causal_pooled_list[["AFR_ANC"]][,(n_vars+1):(2*n_vars)]
afr_A_1 <- tag_pooled_list[["AFR_ANC"]][,1:n_vars]
afr_A_2 <- tag_pooled_list[["AFR_ANC"]][,(n_vars+1):(2*n_vars)]

## We define R2_VEC to make sure the tagging variant R2 is approximately 
## what we observe in our analysis. 
if (PHENOTYPE=="Standing-Height") {
  R2_VEC <- round(c(0.15,0.2,0.25,0.3,0.35,
                    0.4,0.45,0.5,0.55,0.6)/EUR_MEAN_SQ_LD,
                  digits=3)
} else if (PHENOTYPE=="Weight") {
  R2_VEC <- round(c(0.05,0.055,0.06,0.065,0.07,
                    0.075,0.08,0.085,0.09,0.095)/EUR_MEAN_SQ_LD,
                  digits=3)
} else if (PHENOTYPE=="Body-Mass-Index") {
  R2_VEC <- round(c(0.035,0.04,0.045,0.05,0.055,
                    0.06,0.065,0.07,0.075,0.08)/EUR_MEAN_SQ_LD,
                  digits=3)
} else if (PHENOTYPE=="Triglycerides") {
  R2_VEC <- round(c(0.04,0.045,0.05,0.055,0.06,
                    0.065,0.07,0.075,0.08,0.085)/EUR_MEAN_SQ_LD,
                  digits=3)
} else if (PHENOTYPE=="Neutrophil_Count") {
  R2_VEC <- round(c(0.02,0.025,0.03,0.035,0.04,
                    0.045,0.05,0.055,0.06,0.065)/EUR_MEAN_SQ_LD,
                  digits=3)
} else if (PHENOTYPE=="Platelet-counts") {
  R2_VEC <- round(c(0.075,0.08,0.085,0.09,0.095,
                    0.1,0.105,0.11,0.115,0.12)/EUR_MEAN_SQ_LD,
                  digits=3)
} else {
  stop("PHENOTYPE not well-defined. Check trailing argument in Rscript.")
}

## 2. Run simulations and calculations =========================================
biggest_res_df <- NULL
for (R2 in R2_VEC) {
  for (RHO in RHO_VEC) {
    message("------------ ",date(), 
            ": Running simulations with causal variant r2 = ", 
            R2, " and Rho = ", RHO, " ------------")

    # Compute (tau', sigma'2_Eur, sigma'2_Afr) for (r2,rho) choice
    TAU_prime <- RHO*R2/(2*sqrt(
      sum(emp_af_list$fcE*(1-emp_af_list$fcE))*
        sum(emp_af_list$fcA*(1-emp_af_list$fcA)))) 
    SIGMA_prime2_Eur <- R2/(2*sum(emp_af_list$fcE*(1-emp_af_list$fcE)))
    SIGMA_prime2_Afr <- R2/(2*sum(emp_af_list$fcA*(1-emp_af_list$fcA)))
    
    # Compute analytical average CAUSAL effect correlations once per R2,RHO pair
    omega1_prime_vec <- omega_list$OMEGA1_PRIME
    omega2_prime_vec <- omega_list$OMEGA2_PRIME
    omega3_prime_vec <- omega_list$OMEGA3_PRIME
    omega4_prime_vec <- omega_list$OMEGA4_PRIME
    
    u_prime_vec <- SIGMA_prime2_Eur*omega1_prime_vec^2+
      2*TAU_prime*omega1_prime_vec*omega2_prime_vec+
      SIGMA_prime2_Afr*omega2_prime_vec^2
    v_prime_vec <- SIGMA_prime2_Eur*omega3_prime_vec^2+
      2*TAU_prime*omega3_prime_vec*omega4_prime_vec+
      SIGMA_prime2_Afr*omega4_prime_vec^2
    w_prime_vec <- SIGMA_prime2_Eur*omega1_prime_vec*omega3_prime_vec+
      TAU_prime*omega1_prime_vec*omega4_prime_vec+
      TAU_prime*omega2_prime_vec*omega3_prime_vec+
      SIGMA_prime2_Afr*omega2_prime_vec*omega4_prime_vec
    
    LAACor_glo_prime_vec <- w_prime_vec/sqrt(u_prime_vec*v_prime_vec) 
    
    # Compute analytical average TAGGING effect correlations once per R2,RHO pair
    omega1_vec <- omega_list$OMEGA1
    omega2_vec <- omega_list$OMEGA2
    omega3_vec <- omega_list$OMEGA3
    omega4_vec <- omega_list$OMEGA4
    
    theta_Afr_vec <- emp_af_list$rctA*
      sqrt(emp_af_list$fcA*(1-emp_af_list$fcA))/sqrt(emp_af_list$ftA*(1-emp_af_list$ftA))
    
    #[!] 5/4/25 -- Don't use eur_test_cor_vec; actually, should use it!
    # theta_Eur_vec <- eur_test_cor_vec*
    #   sqrt(emp_af_list$fcE*(1-emp_af_list$fcE))/sqrt(emp_af_list$ftE*(1-emp_af_list$ftE))
    #theta_Eur_vec <- emp_af_list$rctE*
    #  sqrt(emp_af_list$fcE*(1-emp_af_list$fcE))/sqrt(emp_af_list$ftE*(1-emp_af_list$ftE))
    theta_Eur_vec <- eur_test_cor_vec*
      sqrt(unadmixed_emp_af_list$fc*(1-unadmixed_emp_af_list$fc))/sqrt(unadmixed_emp_af_list$ft*(1-unadmixed_emp_af_list$ft))
    
    u_vec <- SIGMA_prime2_Eur*theta_Eur_vec^2*omega1_vec^2+
      2*TAU_prime*theta_Eur_vec*theta_Afr_vec*omega1_vec*omega2_vec+
      SIGMA_prime2_Afr*theta_Afr_vec^2*omega2_vec^2
    v_vec <- SIGMA_prime2_Eur*theta_Eur_vec^2*omega3_vec^2+
      2*TAU_prime*theta_Eur_vec*theta_Afr_vec*omega3_vec*omega4_vec+
      SIGMA_prime2_Afr*theta_Afr_vec^2*omega4_vec^2
    w_vec <- SIGMA_prime2_Eur*theta_Eur_vec^2*omega1_vec*omega3_vec+
      TAU_prime*theta_Eur_vec*theta_Afr_vec*omega1_vec*omega4_vec+
      TAU_prime*theta_Eur_vec*theta_Afr_vec*omega2_vec*omega3_vec+
      SIGMA_prime2_Afr*theta_Afr_vec^2*omega2_vec*omega4_vec
    
    LAACor_glo_vec <- w_vec/sqrt(u_vec*v_vec) 
    LAACor_loc_coeff <- sum(theta_Afr_vec*theta_Eur_vec,na.rm=T)/
      (sqrt(sum(theta_Eur_vec^2,na.rm=T))*sqrt(sum(theta_Afr_vec^2,na.rm=T)))
      
    big_res_df <- NULL
    for (rep in 1:100) {
      set.seed(rep)
      ## 2. Simulate causal and tagging effects 
      beta_list <- getEffectVecs(r2=R2, 
                                 fcft_list=emp_af_list,
                                 cross_pop_rho=RHO,
                                 unadmixed_fcft_list=unadmixed_emp_af_list)
      
      # Obtain variant-specific simulated effects
      spop_beta_vec = beta_list[["EUR_TAG"]]/
        sqrt(sum(emp_af_list$ftE*(1-emp_af_list$ftE)))
      tpop_beta_vec = beta_list[["AFR_TAG"]]/
        sqrt(sum(emp_af_list$ftA*(1-emp_af_list$ftA)))
      spop_causal_beta_vec = beta_list[["EUR_CAUSAL"]]/
        sqrt(sum(emp_af_list$fcE*(1-emp_af_list$fcE)))
      tpop_causal_beta_vec = beta_list[["AFR_CAUSAL"]]/
        sqrt(sum(emp_af_list$fcA*(1-emp_af_list$fcA)))
      
      ## 3. Calculate Global model average CAUSAL effect correlation 
      mixed_beta_prime_mat <- outer((1-ind_causal_glo_AFR_anc),spop_causal_beta_vec) + 
        outer(ind_causal_glo_AFR_anc,tpop_causal_beta_vec) # n x p
      num_LAAAfr_prime_hat <- colSums(mixed_beta_prime_mat*afr_A_prime_1,na.rm=T) +
        colSums(mixed_beta_prime_mat*afr_A_prime_2,na.rm=T) # p x 1 
      denom_LAAAfr_prime_hat <- colSums(afr_A_prime_1,na.rm=T) + colSums(afr_A_prime_2,na.rm=T) # p x 1 
      num_LAAEur_prime_hat <- colSums(mixed_beta_prime_mat*(1-afr_A_prime_1),na.rm=T) +
        colSums(mixed_beta_prime_mat*(1-afr_A_prime_2),na.rm=T) # p x 1 
      denom_LAAEur_prime_hat <- 2*n_inds-denom_LAAAfr_prime_hat # p x 1 
      LAAEur_prime_hat <- num_LAAEur_prime_hat/denom_LAAEur_prime_hat
      LAAAfr_prime_hat <- num_LAAAfr_prime_hat/denom_LAAAfr_prime_hat
      
      ## 3. Calculate Global model average TAGGING effect correlation 
      mixed_beta_mat <- outer((1-ind_tag_glo_AFR_anc),spop_beta_vec) + 
        outer(ind_tag_glo_AFR_anc,tpop_beta_vec) # n x p
      num_LAAAfr_hat <- colSums(mixed_beta_mat*afr_A_1,na.rm=T) +
        colSums(mixed_beta_mat*afr_A_2,na.rm=T) # p x 1 
      denom_LAAAfr_hat <- colSums(afr_A_1,na.rm=T) + colSums(afr_A_2,na.rm=T) # p x 1 
      num_LAAEur_hat <- colSums(mixed_beta_mat*(1-afr_A_1),na.rm=T) +
        colSums(mixed_beta_mat*(1-afr_A_2),na.rm=T) # p x 1 
      denom_LAAEur_hat <- 2*n_inds-denom_LAAAfr_hat # p x 1 
      LAAEur_hat <- num_LAAEur_hat/denom_LAAEur_hat
      LAAAfr_hat <- num_LAAAfr_hat/denom_LAAAfr_hat
      
      est_quantities <- data.frame(R2_Causal=R2,
                                   Rho=RHO,
                                   Global_LAACor_prime_hat=cor(LAAEur_prime_hat,LAAAfr_prime_hat),
                                   Mean_Global_LAACor_prime_j_True=sum(w_prime_vec,na.rm=T)/
                                     sqrt(sum(u_prime_vec,na.rm=T)*sum(v_prime_vec,na.rm=T)),
                                   Mean_Global_LAACor_prime_j_Approx=mean(LAACor_glo_prime_vec,na.rm=T),
                                   Local_LAACor_prime_hat=cor(spop_causal_beta_vec,tpop_causal_beta_vec),
                                   Mean_Local_LAACor_prime_j=RHO,
                                   Global_LAACor_hat=cor(LAAEur_hat,LAAAfr_hat),
                                   Mean_Global_LAACor_j_True=sum(w_vec,na.rm=T)/
                                     sqrt(sum(u_vec,na.rm=T)*sum(v_vec,na.rm=T)),
                                   Mean_Global_LAACor_j_Approx=mean(LAACor_glo_vec,na.rm=T),
                                   Local_LAACor_hat=cor(spop_beta_vec,tpop_beta_vec),
                                   Mean_Local_LAACor_j=RHO*LAACor_loc_coeff,
                                   a_bar_tag=a_bar_tag,
                                   a_bar_causal=a_bar_causal,
                                   nvars=length(emp_af_list$fcA))
      
      big_res_df <- rbind(big_res_df,est_quantities)
    }
    
    ## Save
    message(date(), ": Saving file......")
    readr::write_csv(big_res_df,
                     file=paste0(result_dir,PHENOTYPE,"_R2-",concatenate_digits(R2),
                                 "_Rho-",RHO,"_LAACor_results.csv"))
    biggest_res_df <- rbind(biggest_res_df,big_res_df)
  }
}

## 5. Collate everything into a large dataframe ================================
message("------------ ", date(), ": End of simulation ------------")
readr::write_csv(biggest_res_df,
                 file=paste0(result_dir,PHENOTYPE,
                             "_all_LAACor_results.csv"))
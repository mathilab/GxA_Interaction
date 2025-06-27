################################################################################
######################## Causal-Tagging Simulations ############################
######################## Phenotype-specific simulations ########################
################################################################################
#' 
#' This script computes expected values of Cor(ParPGS,Y)^2 and Cor(TotPGS,Y)^2 for 
#' the local model, given simulation parameters rho and r2. Multiple approximations
#' are considered to evaluate the relative accuracy of each approximation scheme. 
#' 
#' Compared to v8, this version also uses causal allele frequencies to simulate
#' phenotypes in simPGSandPheno, fixing a tiny inconsistency in the previous 
#' iteration, where tagging allele frequencies were used to demean BOTH causal
#' and tagging haplotype matrices.
#' 

#### Define Libraries and Directories ------------------------------------------
working_dir <- "/project/mathilab/aaw/admix_prs/"
eur_geno_backup_dir <-
  "/project/mathilab/aaw/admix_prs/PMBB/Synthetic_Data/Eur_Genomes/Causal_Tag_Data/backup/"
adm_geno_dir <- 
  "/project/mathilab/aaw/admix_prs/PMBB/v2/Causal_Tag_Data/"
result_dir <- "/project/mathilab/aaw/admix_prs/results/042325/"

library(dplyr)
library(bigsnpr)

# Read argument to select relevant combination row
args = commandArgs(trailingOnly=TRUE)

# Read argument to select relevant phenotype and R2_VEC
s <- as.numeric(args[1])
u <- as.numeric(args[2])

PHENOTYPE_VEC <- c("Standing-Height","Weight","Body-Mass-Index",
                   "Triglycerides","Neutrophil_Count","Platelet-counts")
PHENOTYPE <- PHENOTYPE_VEC[s]
MODEL <- ifelse(u==1,"local","global")

RHO_VEC <- c(0.92,0.94,0.96,0.98,1)

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
  eur_tag_beta_vec <- eur_causal_beta_vec*unadmixed_fcft_list$rct*
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

## Simulate phenotype and PGSs =================================================
#' Similar to simPGS2, but this function uses causal variants and causal effects
#' to generate the phenotype, while using tagging variants and their effects
#' to compute polygenic scores. Phenotypes are simulated without consideration
#' of local ancestries, owing to shared causal effects across all individuals. 
#' Genotype-phenotype r2 is controlled by normalizing effects by allele frequencies
#' similar to how it's done in our other simulations. Tagging effects are 
#' dependent on local ancestry, so we compute polygenic scores using local 
#' ancestry matrices. We demean all alleles, similar to how we do so for the 
#' rest of our study. 
simPGSandPheno <- function(r2,
                           admixed_cohort_list,
                           ind_target_pop_frac,
                           causal_af,
                           source_pop_af,
                           target_pop_af,
                           source_pop_causal_af,
                           target_pop_causal_af,
                           source_causal_beta_vec,
                           target_causal_beta_vec,
                           source_pop_beta_vec,
                           target_pop_beta_vec,
                           model="local") {
  
  # >>>>>>>>>>>>>>>>>>>>>>>>> Step 0: Extract objects <<<<<<<<<<<<<<<<<<<<<<<<<<
  causal_lanc_mat <- admixed_cohort_list[["CAUSAL_LANC"]] # n x 2p
  #ind_causal_lanc_frac <- rowMeans(causal_lanc_mat) # n 
  tag_lanc_mat <- admixed_cohort_list[["TAG_LANC"]] # n x 2p
  ind_tag_lanc_frac <- rowMeans(tag_lanc_mat,na.rm=TRUE) # n 
  causal_hap_mat <- admixed_cohort_list[["CAUSAL_HAP"]] # n x 2p
  tag_hap_mat <- admixed_cohort_list[["TAG_HAP"]] # n x 2p
  causal_allele_freq <- c(causal_af,causal_af) # p --> 2p
  source_pop_allele_freq <- c(source_pop_af,source_pop_af) # p --> 2p
  target_pop_allele_freq <- c(target_pop_af,target_pop_af) # p --> 2p
  source_pop_causal_allele_freq <- c(source_pop_causal_af,source_pop_causal_af) # p --> 2p
  target_pop_causal_allele_freq <- c(target_pop_causal_af,target_pop_causal_af) # p --> 2p
  big_source_causal_beta_vec <- c(source_causal_beta_vec,source_causal_beta_vec) # p --> 2p
  big_target_causal_beta_vec <- c(target_causal_beta_vec,target_causal_beta_vec) # p --> 2p
  big_source_pop_beta_vec <- c(source_pop_beta_vec,source_pop_beta_vec) # p --> 2p
  big_target_pop_beta_vec <- c(target_pop_beta_vec,target_pop_beta_vec) # p --> 2p
  
  # >>>>>>>>>>>>>>>>>>>>>>> Step 1: Simulate phenotypes <<<<<<<<<<<<<<<<<<<<<<<<
  # In this version, we use LANC-specific causal allele frequencies 
  #causal_hap_mat_demean <- sweep(causal_hap_mat,2,causal_allele_freq)
  #causal_X_demean_beta <- as.numeric(causal_hap_mat_demean %*% big_causal_beta_vec)
  causal_source_pop_masked_mat_demean <- (1-causal_lanc_mat)*sweep(causal_hap_mat,2,
                                                                   source_pop_causal_allele_freq)
  causal_target_pop_masked_mat_demean <- causal_lanc_mat*sweep(causal_hap_mat,2,
                                                               target_pop_causal_allele_freq)
  
  # [!] >>>>>>>> Convert NAs to zeros for purposes of computing PGS <<<<<<<< [!]
  causal_source_pop_masked_mat_demean[is.na(causal_source_pop_masked_mat_demean)] <- 0
  causal_target_pop_masked_mat_demean[is.na(causal_target_pop_masked_mat_demean)] <- 0
  
  if (model=="local") {
    causal_X_demean_beta <- 
      as.numeric(causal_source_pop_masked_mat_demean %*% big_source_causal_beta_vec) +
      as.numeric(causal_target_pop_masked_mat_demean %*% big_target_causal_beta_vec)
    adm_causal_emp_var <- var(causal_X_demean_beta,na.rm=T)
    if (adm_causal_emp_var >= 1) {
      message("ALERT: adm_causal_emp_var is not less than 1. Reset to 0.999")
      adm_causal_emp_var <- 0.999
    }
    
    causal_adm_noise_vec <- rnorm(n=nrow(causal_hap_mat),
                                  mean=0,
                                  sd=sqrt(1-adm_causal_emp_var)) 
    y <- causal_X_demean_beta+causal_adm_noise_vec # Y in ADM now has var(Y) ~ 1
  } else if (model=="global") {
    # [!] 4/4/2025 --- GLOBAL MODEL formula
    n_ind <- length(ind_target_pop_frac)
    mixed_beta_mat <- outer((1-ind_target_pop_frac),big_source_causal_beta_vec) + 
      outer(ind_target_pop_frac,big_target_causal_beta_vec)
    
    # (Demeaned X)*b is used to construct y
    causal_X_demean_beta_source_pop <- as.numeric(rowSums(
      causal_source_pop_masked_mat_demean*mixed_beta_mat, na.rm=TRUE)) 
    causal_X_demean_beta_target_pop <- as.numeric(rowSums(
      causal_target_pop_masked_mat_demean*mixed_beta_mat, na.rm=TRUE))
    causal_X_demean_beta <- causal_X_demean_beta_source_pop + 
      causal_X_demean_beta_target_pop
    
    # Compute noise vector for global model using analytical quantities
    adm_causal_emp_var <- var(causal_X_demean_beta,na.rm=T)
    if (adm_causal_emp_var >= 1) {
      message("ALERT: adm_causal_emp_var is not less than 1. Reset to 0.999")
      adm_causal_emp_var <- 0.999
    }
    
    causal_adm_noise_vec <- rnorm(n=nrow(causal_hap_mat),
                                  mean=0,
                                  sd=sqrt(1-adm_causal_emp_var)) 
    
    # Compute y
    y <- causal_X_demean_beta + causal_adm_noise_vec
  } else {
    stop("Please specify 'local' or 'global' for model parameter.")
  }
  
  # >>>>>>>>>>>>>>>>>>>> Step 2: Compute ParPGS and TotPGS <<<<<<<<<<<<<<<<<<<<<
  # Compute masked matrices
  #source_pop_masked_mat <- (1-tag_lanc_mat)*tag_hap_mat
  #target_pop_masked_mat <- tag_lanc_mat*tag_hap_mat
  
  # Demean matrices 
  source_pop_masked_mat_demean <- (1-tag_lanc_mat)*sweep(tag_hap_mat,2,
                                                         source_pop_allele_freq)
  target_pop_masked_mat_demean <- tag_lanc_mat*sweep(tag_hap_mat,2,
                                                     target_pop_allele_freq)
  
  # [!] >>>>>>>> Convert NAs to zeros for purposes of computing PGS <<<<<<<< [!]
  source_pop_masked_mat_demean[is.na(source_pop_masked_mat_demean)] <- 0
  target_pop_masked_mat_demean[is.na(target_pop_masked_mat_demean)] <- 0
  
  # Compute ParPGS and TotPGS using tagging variants 
  parPGS_demean <- as.numeric( # this is parPGS (demean)
    source_pop_masked_mat_demean%*%big_source_pop_beta_vec)
  totPGS_demean <- parPGS_demean +
    as.numeric(target_pop_masked_mat_demean%*%big_source_pop_beta_vec)
  admPGS_demean <- parPGS_demean +
    as.numeric(target_pop_masked_mat_demean%*%big_target_pop_beta_vec)
  
  # Return
  return(list(Y=y, # Xb + noise 
              G=causal_X_demean_beta, # true Xb
              totPGS_demean=totPGS_demean, # totPGS from tagging variants
              parPGS_demean=parPGS_demean, # parPGS from tagging variants
              admPGS_demean=admPGS_demean, # admPGS from tagging variants
              est_target_pop_frac=mean(ind_tag_lanc_frac))) # cohort a_bar (tagging)
} 

## Estimate r^2 using R2_PGS computed on European genotypes ====================
getR2PGS <- function(beta_vec,
                     r2,
                     unadmixed_fcft_list) {
  # Get scaled effect sizes
  causal_beta_vec_g <- beta_vec/unadmixed_fcft_list$causal_beta_denom
  
  # Compute causal variant genetic score vector
  # take hat(X)*beta
  euram_gene_score <- big_prodVec(X=unadmixed_fcft_list$G_demean_causal,
                                  y.col=causal_beta_vec_g)
  
  # Compute variance of causal variant genetic score vector
  emp_var <- var(euram_gene_score)
  
  # Add noise term to get phenotype
  epsilons_theory <- rnorm(length(euram_gene_score),mean=0,sd=sqrt(1-r2))
  epsilons_emp <- rnorm(length(euram_gene_score),mean=0,sd=sqrt(1-emp_var))
  euram_phenos_theory <- euram_gene_score + epsilons_theory
  euram_phenos_emp <- euram_gene_score + epsilons_emp
  
  # Generate tagging variant effect sizes 
  tag_beta_vec <- beta_vec*unadmixed_fcft_list$rct*
    sqrt((unadmixed_fcft_list$fc*(1-unadmixed_fcft_list$fc))/
           (unadmixed_fcft_list$ft*(1-unadmixed_fcft_list$ft)))
  # Set NA to 0. This happens when causal variant empirical allele freq is 0
  tag_beta_vec[which(is.na(tag_beta_vec))] <- 0
  tag_beta_vec_g <- tag_beta_vec/unadmixed_fcft_list$tag_beta_denom
  
  # Compute PGS from tagging variants
  euram_pgs <- big_prodVec(X=unadmixed_fcft_list$G_demean_tag,
                           y.col=tag_beta_vec_g)
  
  # Compute R2
  r2_pgs_theory <- summary(lm(euram_phenos_theory~euram_pgs))$r.squared
  r2_pgs_emp <- summary(lm(euram_phenos_emp~euram_pgs))$r.squared
  r2_pgs_causal_theory <- summary(lm(euram_phenos_theory~euram_gene_score))$r.squared
  r2_pgs_causal_emp <- summary(lm(euram_phenos_emp~euram_gene_score))$r.squared
  
  # Return
  # first entry is relevant; second is just for tracking
  return(list(R2VEC=c(r2_pgs_theory,
                      r2_pgs_emp,
                      r2_pgs_causal_theory,
                      r2_pgs_causal_emp),
              BETA_TAG_UNADMIXED_EUR=tag_beta_vec_g)) 
}

## getAdmixedHapList (Inherited from master script) ============================
#' Compute key quantities from seed data
#' 
#' Given seed data (both allelic dosages and ancestry calls matrices), this
#' function computes the following quantities: 
#' (1) population-specific allele frequencies per haplotype; 
#' (2) number of markers; 
#' (3) number of individuals; 
#' (4) ancestry-specific demeaned allelic dosage matrices; 
#' (5) global ancestry (in target population) fractions across individuals;
#' (6) global ancestry (in target population) fractions across markers; 
#' (7) global ancestry fraction across cohort
#' 
#' Additionally, the original admixed haplotype matrix and target population
#' ancestry call matrix are returned. 
#' 
#' This function was modified on 10/8/24 to compute genotype allele frequencies
#' instead of haplotype-specific allele frequencies.
getAdmixedHapList <- function(admixed_haps,
                              afr_ancs) {
  # Compute masked matrices
  target_pop_masked_mat <- admixed_haps*afr_ancs
  source_pop_masked_mat <- admixed_haps*(1-afr_ancs)
  
  # Get 2 x p so that can halve later
  two_p <- ncol(source_pop_masked_mat)
  
  # Compute pop-specific allele frequencies
  # [!] 10/8/2024 - use genotype allele frequencies rather than haplotypes
  source_pop_num <- colSums(source_pop_masked_mat,na.rm=TRUE)
  source_pop_num <- source_pop_num[1:(two_p/2)] + source_pop_num[(two_p/2+1):two_p]
  source_pop_denom <- colSums((1-afr_ancs),na.rm=TRUE) 
  source_pop_denom <- source_pop_denom[1:(two_p/2)] + 
    source_pop_denom[(two_p/2+1):two_p]
  source_pop_allele_freq <- source_pop_num/source_pop_denom
  source_pop_allele_freq <- c(source_pop_allele_freq,source_pop_allele_freq) 
  
  target_pop_num <- colSums(target_pop_masked_mat,na.rm=TRUE)
  target_pop_num <- target_pop_num[1:(two_p/2)] + target_pop_num[(two_p/2+1):two_p]
  target_pop_denom <- colSums(afr_ancs,na.rm=TRUE)
  target_pop_denom <- target_pop_denom[1:(two_p/2)] + 
    target_pop_denom[(two_p/2+1):two_p]
  target_pop_allele_freq <- target_pop_num/target_pop_denom
  target_pop_allele_freq <- c(target_pop_allele_freq,target_pop_allele_freq)
  
  n_markers <- two_p/2
  
  # Extract target population ancestry fraction markers
  marker_target_pop_frac <- colMeans(afr_ancs,na.rm=TRUE)
  ind_target_pop_frac <- rowMeans(afr_ancs,na.rm=TRUE)
  n_ind <- length(ind_target_pop_frac)
  
  # Estimate overall target population fraction of cohort
  cohort_est_target_pop_frac <- mean(marker_target_pop_frac)
  
  # Define source population ancestry matrix 
  eur_ancs <- 1-afr_ancs
  
  source_pop_masked_mat_demean <- eur_ancs* # n x 2p
    sweep(admixed_haps,2,source_pop_allele_freq)
  target_pop_masked_mat_demean <- afr_ancs* # n x 2p
    sweep(admixed_haps,2,target_pop_allele_freq) 
  
  # Define admixed_haps_list
  ADMIXED_HAPS_LIST <- list(HAP_MAT=admixed_haps, # n x 2p matrix
                            SOURCE_POP_AF=source_pop_allele_freq, # length = 2p
                            TARGET_POP_AF=target_pop_allele_freq, # length = 2p
                            SOURCE_POP_DEMEAN_HAP_MAT=source_pop_masked_mat_demean, # n x 2p matrix
                            TARGET_POP_DEMEAN_HAP_MAT=target_pop_masked_mat_demean, # n x 2p matrix
                            SOURCE_POP_ANC_MAT=eur_ancs, # n x 2p matrix
                            TARGET_POP_ANC_MAT=afr_ancs, # n x 2p matrix
                            MARKER_TARGET_POP_FRAC=marker_target_pop_frac, # length = 2p
                            IND_TARGET_POP_FRAC=ind_target_pop_frac, # length = n
                            N_IND = n_ind, # scalar
                            N_MARKERS = n_markers, # scalar
                            TPOP_FRAC=cohort_est_target_pop_frac) # scalar
  
  # Return
  return(ADMIXED_HAPS_LIST)
}

## getExpectedBetas ============================================================
#'
#' 
#' Calculates expected values of the following quantities:
#' - (beta_j^Eur)^2
#' - beta_j^Eur*beta_j'^Eur
#' - beta_j^Eur*beta_j'^Afr
#' 
getExpectedBetas <- function(eur_causal_af,
                             eur_tag_af,
                             afr_causal_af,
                             eur_causal_tag_ld,
                             r2,
                             rho) {
  # Compute terms that are used in all expressions
  eur_causal_var <- eur_causal_af*(1-eur_causal_af) # p x 1
  eur_tag_var <- eur_tag_af*(1-eur_tag_af) # p x 1
  afr_causal_var <- afr_causal_af*(1-afr_causal_af) # p x 1
  beta_Eur_terms_denom <- 2*sum(eur_causal_var) # scalar
  beta_Eur_beta_Afr_term_denom <- 2*sqrt(sum(afr_causal_var)*sum(eur_causal_var)) # scalar
  
  # Eur Tagging effect squared, (beta_j^Eur)^2
  beta_eur_sq <- (r2*eur_causal_tag_ld^2*eur_causal_var/eur_tag_var)/
    beta_Eur_terms_denom
  # Eur Causal effect * Tagging effect, beta_j^Eur*beta_j'^Eur
  beta_eur_beta_eur_prime <- (r2*eur_causal_tag_ld*sqrt(eur_causal_var/eur_tag_var))/
    beta_Eur_terms_denom
  # Afr Causal effect * Eur Tagging effect, beta_j^Eur*beta_j'^Afr
  beta_eur_beta_afr_prime <- (r2*rho*eur_causal_tag_ld^2*eur_causal_var/eur_tag_var)/
    beta_Eur_beta_Afr_term_denom
  
  # Return
  return(list(BETA_EUR_SQ=beta_eur_sq,
              BETA_EUR_x_BETA_EUR_PRIME=beta_eur_beta_eur_prime,
              BETA_EUR_x_BETA_AFR_PRIME=beta_eur_beta_afr_prime))
}

## getApproxVals ===============================================================
#' 
#' Calculates expected values of squared correlations from theory, using various 
#' approximation schemes. Depends on the functions checkEURCor and   
#' 
#' The following quantities are key and are annotated by the object that contains
#' it ([1] = adm_cohort_list, [2] = emp_af_list, [3] = unadmixed_emp_af_list)
#' - Fixed data objects
#'   - Causal local ancestry matrices (Afr & Eur) [1]
#'   - Tagging local ancestry matrices (Afr & Eur) [1]
#'   - Causal variant haplotype matrices [1]
#'   - Tagging variant haplotype matrices [1]
#'   - Causal ancestry-specific AF (Afr & Eur) [2]
#'   - Tagging ancestry-specific AF (Afr & Eur) [2]
#'   - Causal-Tagging LD (in admixed cohort, conditioned by local ancestry) [2]
#'   - Causal-Tagging LD (in unadmixed cohort) [3]
#' - Problem parameters
#'   - r2 
#'   - rho 
#'   
getApproxVals <- function(adm_cohort_list,
                          emp_af_list,
                          unadmixed_emp_af_list,
                          r2,
                          rho) {
  # Extract fixed data objects
  afr_A_prime <- adm_cohort_list$CAUSAL_LANC # n x 2p
  afr_A <- adm_cohort_list$TAG_LANC # n x 2p
  eur_A <- 1-afr_A
  eur_A_prime <- 1-afr_A_prime
  causal_hap_mat <- adm_cohort_list$CAUSAL_HAP # n x 2p
  tag_hap_mat <- adm_cohort_list$TAG_HAP # n x 2p 
  afr_f_prime <- emp_af_list$fcA # p x 1
  afr_f <- emp_af_list$ftA # p x 1
  eur_f_prime <- emp_af_list$fcE # p x 1
  eur_f <- emp_af_list$ftE # p x 1
  eur_rct_unadmixed <- checkEURCor(causal_geno=unadmixed_emp_af_list$G_demean_causal,
                                   tag_geno=unadmixed_emp_af_list$G_demean_tag) # p x 1
  afr_rct_admixed <- emp_af_list$rctA # p x 1
  eur_rct_admixed <- emp_af_list$rctE # p x 1
  
  causal_eur_masked_mat_demean <- (1-afr_A_prime)*sweep(causal_hap_mat,2,
                                                        c(eur_f_prime,eur_f_prime))
  causal_afr_masked_mat_demean <- afr_A_prime*sweep(causal_hap_mat,2,
                                                    c(afr_f_prime,afr_f_prime))
  tag_eur_masked_mat_demean <- (1-afr_A)*sweep(tag_hap_mat,2,
                                               c(eur_f,eur_f))
  tag_afr_masked_mat_demean <- afr_A*sweep(tag_hap_mat,2,
                                           c(afr_f,afr_f))
  # Get ancestry-specific AF-demeaned haplotypes
  X_prime <- causal_eur_masked_mat_demean+causal_afr_masked_mat_demean # demeaned, n x 2p
  X <- tag_eur_masked_mat_demean+tag_afr_masked_mat_demean # demeaned, n x 2p
  
  # Get number of variants
  n_vars <- dim(X)[2]/2
  
  # Make individual-level global ancestry matrix n x p
  ind_global_anc_vec <- rowMeans(afr_A,na.rm=T) # global Afr ancestry based on causal loci 
  ind_global_anc_mat <- matrix(rep(ind_global_anc_vec, times = n_vars), # n x p
                               nrow = length(ind_global_anc_vec), ncol = n_vars)
  
  # Compute expected values of effects, which appear throughout all expressions
  expected_betas_list <- getExpectedBetas(eur_causal_af=eur_f_prime,
                                          eur_tag_af=eur_f,
                                          afr_causal_af=afr_f_prime,
                                          eur_causal_tag_ld=eur_rct_unadmixed,
                                          r2=r2,
                                          rho=rho)
  E_beta_eur_sq <- expected_betas_list$BETA_EUR_SQ # p x 1
  E_beta_eur_beta_eur_prime <- expected_betas_list$BETA_EUR_x_BETA_EUR_PRIME # p x 1
  E_beta_eur_beta_afr_prime <- expected_betas_list$BETA_EUR_x_BETA_AFR_PRIME # p x 1
  
  # Terms 1 to 6 
  eur_tag_var <- eur_f*(1-eur_f) # p x 1
  afr_tag_var <- afr_f*(1-afr_f) # p x 1
  a_dot_js <- colMeans(afr_A,na.rm=T) # 2p x 1
  one_minus_a_dot_js <- colMeans(1-afr_A,na.rm=T) # 2p x 1
  cross_cov_tag_hap <- colMeans(X[,1:n_vars]*X[,(n_vars+1):(2*n_vars)],na.rm=T) # p x 1
  cross_cov_eur_lanc_tag_hap <- colMeans(eur_A[,1:n_vars]*
                                           eur_A[,(n_vars+1):(2*n_vars)]*
                                           X[,1:n_vars]*
                                           X[,(n_vars+1):(2*n_vars)],na.rm=T)
    
  term_1_y_2 <- c(E_beta_eur_sq,E_beta_eur_sq)*(
    one_minus_a_dot_js*c(eur_tag_var,eur_tag_var)+a_dot_js*c(afr_tag_var,afr_tag_var))
  term_3 <- 2*E_beta_eur_sq*cross_cov_tag_hap
  
  term_4_y_5 <- c(E_beta_eur_sq,E_beta_eur_sq)*(
    one_minus_a_dot_js*c(eur_tag_var,eur_tag_var))
  term_6 <- 2*E_beta_eur_sq*cross_cov_eur_lanc_tag_hap
  
  # Terms 7L - 10L
  # 7L and 10L can be computed together
  loc_eur_causal_tag_cov <- colMeans(eur_A_prime*X_prime*X,na.rm=T)
  loc_afr_causal_tag_cov <- colMeans(afr_A_prime*X_prime*X,na.rm=T)
  loc_afr_causal2_tag1_cov <- colMeans(afr_A_prime[,(n_vars+1):(2*n_vars)]*
                                         X_prime[,(n_vars+1):(2*n_vars)]*X[,1:n_vars],
                                       na.rm=T)
  loc_afr_causal1_tag2_cov <- colMeans(afr_A_prime[,1:n_vars]*
                                         X_prime[,1:n_vars]*X[,(n_vars+1):(2*n_vars)],
                                       na.rm=T)
  loc_eur_causal2_tag1_cov <- colMeans(eur_A_prime[,(n_vars+1):(2*n_vars)]*
                                         X_prime[,(n_vars+1):(2*n_vars)]*X[,1:n_vars],
                                       na.rm=T)
  loc_eur_causal1_tag2_cov <- colMeans(eur_A_prime[,1:n_vars]*
                                         X_prime[,1:n_vars]*X[,(n_vars+1):(2*n_vars)],
                                       na.rm=T)
  term_7L_y_10L <- c(E_beta_eur_beta_eur_prime,E_beta_eur_beta_eur_prime)*loc_eur_causal_tag_cov+
    c(E_beta_eur_beta_afr_prime,E_beta_eur_beta_afr_prime)*loc_afr_causal_tag_cov
  term_8L <- E_beta_eur_beta_afr_prime*loc_afr_causal2_tag1_cov+ 
    E_beta_eur_beta_eur_prime*loc_eur_causal2_tag1_cov
  term_9L <- E_beta_eur_beta_afr_prime*loc_afr_causal1_tag2_cov+
    E_beta_eur_beta_eur_prime*loc_eur_causal1_tag2_cov
  
  # Terms 11L - 14L
  # 11L and 14L can be computed together
  loc_eur_causal_eur_tag_cov <- colMeans(eur_A_prime*eur_A*X_prime*X,na.rm=T)
  loc_afr_causal_eur_tag_cov <- colMeans(afr_A_prime*eur_A*X_prime*X,na.rm=T)
  loc_afr_causal2_eur_tag1_cov <- colMeans(afr_A_prime[,(n_vars+1):(2*n_vars)]*
                                             X_prime[,(n_vars+1):(2*n_vars)]*
                                             X[,1:n_vars]*eur_A[,1:n_vars],
                                           na.rm=T)
  loc_afr_causal1_eur_tag2_cov <- colMeans(afr_A_prime[,1:n_vars]*
                                             X_prime[,1:n_vars]*
                                             X[,(n_vars+1):(2*n_vars)]*
                                             eur_A[,(n_vars+1):(2*n_vars)],
                                           na.rm=T)
  loc_eur_causal2_eur_tag1_cov <- colMeans(eur_A_prime[,(n_vars+1):(2*n_vars)]*
                                             X_prime[,(n_vars+1):(2*n_vars)]*
                                             X[,1:n_vars]*eur_A[,1:n_vars],
                                           na.rm=T)
  loc_eur_causal1_eur_tag2_cov <- colMeans(eur_A_prime[,1:n_vars]*
                                             X_prime[,1:n_vars]*
                                             X[,(n_vars+1):(2*n_vars)]*
                                             eur_A[,(n_vars+1):(2*n_vars)],
                                           na.rm=T)
  term_11L_y_14L <- c(E_beta_eur_beta_eur_prime,E_beta_eur_beta_eur_prime)*loc_eur_causal_eur_tag_cov+
    c(E_beta_eur_beta_afr_prime,E_beta_eur_beta_afr_prime)*loc_afr_causal_eur_tag_cov
  term_12L <- E_beta_eur_beta_afr_prime*loc_afr_causal2_eur_tag1_cov+ 
    E_beta_eur_beta_eur_prime*loc_eur_causal2_eur_tag1_cov
  term_13L <- E_beta_eur_beta_afr_prime*loc_afr_causal1_eur_tag2_cov+
    E_beta_eur_beta_eur_prime*loc_eur_causal1_eur_tag2_cov
  
  # Terms 7G - 10G
  big_ind_global_anc_mat <- cbind(ind_global_anc_mat,ind_global_anc_mat) # n x 2p 
  # 7G and 10G can be computed together
  glo_eur_anc_causal_tag_cov <- colMeans((1-big_ind_global_anc_mat)*X_prime*X,na.rm=T)
  glo_afr_anc_causal_tag_cov <- colMeans(big_ind_global_anc_mat*X_prime*X,na.rm=T)
  glo_eur_causal1_tag2_cov <- colMeans((1-ind_global_anc_mat)*
                                         X_prime[,1:n_vars]*
                                         X[,(n_vars+1):(2*n_vars)],
                                       na.rm=T)
  glo_eur_causal2_tag1_cov <- colMeans((1-ind_global_anc_mat)*
                                         X_prime[,(n_vars+1):(2*n_vars)]*
                                         X[,1:n_vars],
                                       na.rm=T)
  glo_afr_causal1_tag2_cov <- colMeans(ind_global_anc_mat*
                                         X_prime[,1:n_vars]*
                                         X[,(n_vars+1):(2*n_vars)],
                                       na.rm=T)
  glo_afr_causal2_tag1_cov <- colMeans(ind_global_anc_mat*
                                         X_prime[,(n_vars+1):(2*n_vars)]*
                                         X[,1:n_vars],
                                       na.rm=T)
  term_7G_y_10G <- c(E_beta_eur_beta_eur_prime,E_beta_eur_beta_eur_prime)*glo_eur_anc_causal_tag_cov+
    c(E_beta_eur_beta_afr_prime,E_beta_eur_beta_afr_prime)*glo_afr_anc_causal_tag_cov
  term_8G <- E_beta_eur_beta_afr_prime*glo_afr_causal2_tag1_cov+ 
    E_beta_eur_beta_eur_prime*glo_eur_causal2_tag1_cov
  term_9G <- E_beta_eur_beta_afr_prime*glo_afr_causal1_tag2_cov+
    E_beta_eur_beta_eur_prime*glo_eur_causal1_tag2_cov
  
  # Terms 11G - 14G
  # 11G and 14G can be computed together
  glo_eur_causal_eur_tag_cov <- colMeans((1-big_ind_global_anc_mat)*eur_A*X_prime*X,na.rm=T)
  glo_afr_causal_eur_tag_cov <- colMeans(big_ind_global_anc_mat*eur_A*X_prime*X,na.rm=T)
  glo_afr_causal2_eur_tag1_cov <- colMeans(ind_global_anc_mat*
                                             X_prime[,(n_vars+1):(2*n_vars)]*
                                             X[,1:n_vars]*eur_A[,1:n_vars],
                                           na.rm=T)
  glo_afr_causal1_eur_tag2_cov <- colMeans(ind_global_anc_mat*
                                             X_prime[,1:n_vars]*
                                             X[,(n_vars+1):(2*n_vars)]*
                                             eur_A[,(n_vars+1):(2*n_vars)],
                                           na.rm=T)
  glo_eur_causal2_eur_tag1_cov <- colMeans((1-ind_global_anc_mat)*
                                             X_prime[,(n_vars+1):(2*n_vars)]*
                                             X[,1:n_vars]*eur_A[,1:n_vars],
                                           na.rm=T)
  glo_eur_causal1_eur_tag2_cov <- colMeans((1-ind_global_anc_mat)*
                                             X_prime[,1:n_vars]*
                                             X[,(n_vars+1):(2*n_vars)]*
                                             eur_A[,(n_vars+1):(2*n_vars)],
                                           na.rm=T)
  term_11G_y_14G <- c(E_beta_eur_beta_eur_prime,E_beta_eur_beta_eur_prime)*glo_eur_causal_eur_tag_cov+
    c(E_beta_eur_beta_afr_prime,E_beta_eur_beta_afr_prime)*glo_afr_causal_eur_tag_cov
  term_12G <- E_beta_eur_beta_afr_prime*glo_afr_causal2_eur_tag1_cov+ 
    E_beta_eur_beta_eur_prime*glo_eur_causal2_eur_tag1_cov
  term_13G <- E_beta_eur_beta_afr_prime*glo_afr_causal1_eur_tag2_cov+
    E_beta_eur_beta_eur_prime*glo_eur_causal1_eur_tag2_cov
  
  # Return the following approximations
  # 1. Simplest approximation
  # 2. First order approximation (causal/tagging AF equal) [!]
  # 3. Second order approximation (average LD) [!]
  # 4. Average global ancestry approximation [!]
  sq_Cor_ParPGS_Y_Loc <- (sum(term_11L_y_14L)+sum(term_12L)+sum(term_13L))^2/
    (sum(term_4_y_5)+sum(term_6))
  sq_Cor_ParPGS_Y_Glo <- (sum(term_11G_y_14G)+sum(term_12G)+sum(term_13G))^2/
    (sum(term_4_y_5)+sum(term_6))
  sq_Cor_TotPGS_Y_Loc <- (sum(term_7L_y_10L)+sum(term_8L)+sum(term_9L))^2/
    (sum(term_1_y_2)+sum(term_3))
  sq_Cor_TotPGS_Y_Glo <- (sum(term_7G_y_10G)+sum(term_8G)+sum(term_9G))^2/
    (sum(term_1_y_2)+sum(term_3))
  
  return(data.frame(SQ_COR_PARPGS_LOC=sq_Cor_ParPGS_Y_Loc,
                    SQ_COR_PARPGS_GLO=sq_Cor_ParPGS_Y_Glo,
                    SQ_COR_TOTPGS_LOC=sq_Cor_TotPGS_Y_Loc,
                    SQ_COR_TOTPGS_GLO=sq_Cor_TotPGS_Y_Glo))
}

#### Preamble ------------------------------------------------------------------
#source(paste0(working_dir,"PMBB_sim_master_script_geno_AFs_vers_030425.R"))
message(date(), ": Simulating with ", PHENOTYPE,
        " haplotype/local ancestry matrices under ", MODEL, " model")

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
message("~~~ Allele freq comparisons in ADM Afr: Average f_c-f_t = ", 
        sprintf("%.4f",mean(emp_af_list$fcA-emp_af_list$ftA)))
message("~~~ Allele freq comparisons in ADM Afr: Var(f_c-f_t) = ", 
        sprintf("%.4f",var(emp_af_list$fcA-emp_af_list$ftA)))
message("~~~ Allele freq comparisons in ADM Afr: Average |f_c-f_t| = ", 
        sprintf("%.4f",mean(abs(emp_af_list$fcA-emp_af_list$ftA))))
message(">>> Allele freq comparisons in ADM Eur: Average f_c-f_t = ", 
        sprintf("%.4f",mean(emp_af_list$fcE-emp_af_list$ftE)))
message(">>> Allele freq comparisons in ADM Eur: Var(f_c-f_t) = ", 
        sprintf("%.4f",var(emp_af_list$fcE-emp_af_list$ftE)))
message(">>> Allele freq comparisons in ADM Eur: Average |f_c-f_t| = ", 
        sprintf("%.4f",mean(abs(emp_af_list$fcE-emp_af_list$ftE))))
message("### Tagging Variants Allele freq comparisons in ADM: Average f_tA-f_tE = ", 
        sprintf("%.4f",mean(emp_af_list$ftA-emp_af_list$ftE)))
message("### Tagging Variants Allele freq comparisons in ADM: Var(f_tA-f_tE) = ", 
        sprintf("%.4f",var(emp_af_list$ftA-emp_af_list$ftE)))
message("### Tagging Variants Allele freq comparisons in ADM: Average |f_tA-f_tE| = ", 
        sprintf("%.4f",mean(abs(emp_af_list$ftA-emp_af_list$ftE))))

unadmixed_emp_af_list <- getunadmixedfcft(bigsnpr_G_causal=pheno_EurAm_G,
                                          bigsnpr_G_tag=pheno_EurAm_G_tag,
                                          causal_index=causal_indices, 
                                          tag_index=tag_indices)

eur_test_cor_vec <- checkEURCor(causal_geno=unadmixed_emp_af_list$G_demean_causal,
                                tag_geno=unadmixed_emp_af_list$G_demean_tag)
EUR_MEAN_SQ_LD <- mean(eur_test_cor_vec^2,na.rm=T)
message("Mean squared LD / emp corr^2 between Tagging and Causal in EUR = ", 
        sprintf("%.4f",EUR_MEAN_SQ_LD))
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

## 2. Simulate causal and tagging effects ======================================
## 3. Simulate PGS and phenotypes ==============================================
## 4. Run estimation ===========================================================
biggest_res_df <- NULL
for (R2 in R2_VEC) {
  for (RHO in RHO_VEC) {
    message("------------ ",date(), 
            ": Running simulations with causal variant r2 = ", 
            R2, " and Rho = ", RHO, " ------------")
    # Initialize big_res_df
    big_res_df <- NULL
    
    # Compute analytical quantities once per R2,RHO pair
    E_sq_cor_list <- getApproxVals(adm_cohort_list=list(CAUSAL_LANC=causal_pooled_list[["AFR_ANC"]],
                                                        CAUSAL_HAP=causal_pooled_list[["HAP"]],
                                                        TAG_LANC=tag_pooled_list[["AFR_ANC"]],
                                                        TAG_HAP=tag_pooled_list[["HAP"]]),
                                   emp_af_list=emp_af_list,
                                   unadmixed_emp_af_list=unadmixed_emp_af_list,
                                   r2=R2,
                                   rho=RHO)
    
    for (rep in 1:200) {
      set.seed(rep)
      ## 2. Simulate causal and tagging effects 
      beta_list <- getEffectVecs(r2=R2, 
                                 fcft_list=emp_af_list,
                                 cross_pop_rho=RHO,
                                 unadmixed_fcft_list=unadmixed_emp_af_list)
      
      ## 3. Calculate PGS and Y in ADM cohort
      spop_beta_vec = beta_list[["EUR_TAG"]]/
        sqrt(sum(emp_af_list$ftE*(1-emp_af_list$ftE)))
      # spop_beta_vec = beta_list[["EUR_TAG"]]/
      #   sqrt(sum(unadmixed_emp_af_list$ft*(1-unadmixed_emp_af_list$ft)))
      tpop_beta_vec = beta_list[["AFR_TAG"]]/
        sqrt(sum(emp_af_list$ftA*(1-emp_af_list$ftA)))
      
      spop_causal_beta_vec = beta_list[["EUR_CAUSAL"]]/
        sqrt(sum(emp_af_list$fcE*(1-emp_af_list$fcE)))
      tpop_causal_beta_vec = beta_list[["AFR_CAUSAL"]]/
        sqrt(sum(emp_af_list$fcA*(1-emp_af_list$fcA)))
      
      sim_pgs <- simPGSandPheno(r2=R2,
                                admixed_cohort_list=list(CAUSAL_LANC=causal_pooled_list[["AFR_ANC"]],
                                                         CAUSAL_HAP=causal_pooled_list[["HAP"]],
                                                         TAG_LANC=tag_pooled_list[["AFR_ANC"]],
                                                         TAG_HAP=tag_pooled_list[["HAP"]]),
                                ind_target_pop_frac=ind_tag_glo_AFR_anc,
                                causal_af=emp_af_list$fc_lanc_agnostic,
                                source_pop_af=emp_af_list$ftE,
                                target_pop_af=emp_af_list$ftA,
                                source_pop_causal_af=emp_af_list$fcE,
                                target_pop_causal_af=emp_af_list$fcA,
                                source_causal_beta_vec=spop_causal_beta_vec, 
                                target_causal_beta_vec= tpop_causal_beta_vec,
                                source_pop_beta_vec=spop_beta_vec,
                                target_pop_beta_vec=tpop_beta_vec,
                                model=MODEL)
      message("[!] Cor(admPGS, Y)^2 in ADM cohort = ", 
              sprintf("%.4f",cor(sim_pgs$admPGS_demean,sim_pgs$Y)^2))
      
      ## 4. Calculate R2_PGS from EUR cohort 
      r2_pgs_list <- getR2PGS(beta_vec=beta_list[["EUR_CAUSAL"]],
                              r2=R2,
                              unadmixed_fcft_list=unadmixed_emp_af_list)
      r2_pgs_vec <- r2_pgs_list$R2VEC
      message("Cor(EUR Tagging in unadmixed, EUR Tagging in ADM) = ", 
              sprintf("%.4f",cor(r2_pgs_list$BETA_TAG_UNADMIXED_EUR, 
                                 beta_list$EUR_TAG)))
      message("[!] Tagging SNP R\u00B2_PGS computed for simulated phenotypes in EUR = ", 
              sprintf("%.4f",r2_pgs_vec[1]), 
              "; true r\u00B2 (if using causal SNPs) = ", R2, 
              " [estimated ", sprintf("%.4f",r2_pgs_vec[3]),"]")
      
      ## 5. Compute squared correlations
      # [!] Perform regression of Y against ParPGS and compute variance of residuals
      parPGS_sum_sq_resid <- sum(lm(sim_pgs$Y~sim_pgs$parPGS_demean)$residuals^2)
      totPGS_sum_sq_resid <- sum(lm(sim_pgs$Y~sim_pgs$totPGS_demean)$residuals^2)
      admPGS_sum_sq_resid <- sum(lm(sim_pgs$Y~sim_pgs$admPGS_demean)$residuals^2)
      G_sum_sq_resid <- sum(lm(sim_pgs$Y~sim_pgs$G)$residuals^2)
      
      if (u==1) {
        expected_cor2_y_totPGS_demean <- E_sq_cor_list$SQ_COR_TOTPGS_LOC
        expected_cor2_y_parPGS_demean <- E_sq_cor_list$SQ_COR_PARPGS_LOC
      } else {
        expected_cor2_y_totPGS_demean <- E_sq_cor_list$SQ_COR_TOTPGS_GLO
        expected_cor2_y_parPGS_demean <- E_sq_cor_list$SQ_COR_PARPGS_GLO
      }
      
      est_quantities <- data.frame(var_y=var(sim_pgs$Y),
                                   var_G=var(sim_pgs$G),
                                   var_totPGS_demean=var(sim_pgs$totPGS_demean),
                                   var_parPGS_demean=var(sim_pgs$parPGS_demean),
                                   cor2_y_G=cor(sim_pgs$Y,sim_pgs$G)^2,
                                   cor2_y_totPGS_demean=cor(sim_pgs$Y,sim_pgs$totPGS_demean)^2,
                                   cor2_y_parPGS_demean=cor(sim_pgs$Y,sim_pgs$parPGS_demean)^2,
                                   cor2_y_admPGS_demean=cor(sim_pgs$Y,sim_pgs$admPGS_demean)^2,
                                   cor2_G_admPGS_demean=cor(sim_pgs$G,sim_pgs$admPGS)^2,
                                   E_cor2_y_totPGS_demean=expected_cor2_y_totPGS_demean,
                                   E_cor2_y_parPGS_demean=expected_cor2_y_parPGS_demean,
                                   sumsq_resid_y_regress_G=G_sum_sq_resid,
                                   sumsq_resid_y_regress_totPGS=totPGS_sum_sq_resid,
                                   sumsq_resid_y_regress_parPGS=parPGS_sum_sq_resid,
                                   sumsq_resid_y_regress_admPGS=admPGS_sum_sq_resid,
                                   R2_Causal=R2,
                                   R2_PGS_EUR=r2_pgs_vec[1],
                                   Rho=RHO,
                                   mean_emp_sq_LD_ADM=ADM_MEAN_SQ_LD,
                                   mean_emp_sq_LD_EUR=EUR_MEAN_SQ_LD,
                                   mean_emp_sq_LD_AFR=mean(emp_af_list$rctA^2,na.rm=T),
                                   mean_emp_LD_ADM=mean(test_cor_vec,na.rm=T),
                                   mean_emp_LD_EUR=mean(eur_test_cor_vec,na.rm=T),
                                   mean_emp_LD_AFR=mean(emp_af_list$rctA,na.rm=T),
                                   a_bar_tag=a_bar_tag,
                                   a_bar_causal=a_bar_causal,
                                   cor_beta_causal_eur_tag=cor(beta_list$EUR_CAUSAL, 
                                                               beta_list$EUR_TAG),
                                   cor_beta_causal_afr_tag=cor(beta_list$AFR_CAUSAL, 
                                                               beta_list$AFR_TAG),
                                   cor_eur_tag_afr_tag=cor(beta_list$EUR_TAG, 
                                                           beta_list$AFR_TAG),
                                   cor_eur_tag_EUR_eur_tag=cor(r2_pgs_list$BETA_TAG_UNADMIXED_EUR, 
                                                               beta_list$EUR_TAG),
                                   mean_var_fcA_ADM=mean(emp_af_list$fcA*(1-emp_af_list$fcA),na.rm=T),
                                   mean_var_fcE_ADM=mean(emp_af_list$fcE*(1-emp_af_list$fcE),na.rm=T),
                                   mean_var_fcE_EUR=mean(unadmixed_emp_af_list$fc*
                                                           (1-unadmixed_emp_af_list$fc),na.rm=T),
                                   mean_var_ftA_ADM=mean(emp_af_list$ftA*(1-emp_af_list$ftA),na.rm=T),
                                   mean_var_ftE_ADM=mean(emp_af_list$ftE*(1-emp_af_list$ftE),na.rm=T),
                                   mean_var_ftE_EUR=mean(unadmixed_emp_af_list$ft*
                                                           (1-unadmixed_emp_af_list$ft),na.rm=T),
                                   nvars=length(emp_af_list$fcA))
      
      big_res_df <- rbind(big_res_df,est_quantities)
    }
    
    ## Save
    message(date(), ": Saving file......")
    readr::write_csv(big_res_df,
                     file=paste0(result_dir,PHENOTYPE,"_R2-",concatenate_digits(R2),
                                 "_Rho-",RHO,"_",MODEL,
                                 "_model_diagnostic_results.csv"))
    biggest_res_df <- rbind(biggest_res_df,big_res_df)
  }
}

## 5. Collate everything into a large dataframe ================================
message("------------ ", date(), ": End of simulation ------------")
readr::write_csv(biggest_res_df,
                 file=paste0(result_dir,PHENOTYPE,"_",MODEL,
                             "_model_all_diagnostic_results.csv"))
################################################################################
############################# Simulation Script ################################
################################# PMBB Version #################################
################################################################################

# ----------------------------------- NOTES ------------------------------------
# The reason for this script is that PMBB data has NAs so we have to deal with them.
# Also, we don't need to simulate admixed haplotype matrices since they are given.
# 
# (OLD) Mar 22, 2024: This version assumes haplotype-specific effects (beta_j's). 
# It also applies formulas that treat allelic dosages and ancestral assignments as fixed,
# rather than randomly drawn from some population with a mean population allele frequency.
#
# (NEW) Apr 18, 2024: This version allows user to specify the r^2 (var(G)/var(Y))
# in the European population, which would avoid having the estimate it using the
# PGS. There are multiple ways to estimate r^2: REML, LDSR, etc. See this article:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7330487/.
#
# (NEW) Jul 19, 2024: This version adds a new function that allows computation of 
# correlation of (average) effect size by local ancestry. See computeBetaCor on 
# line 393.
# 
# (NEW) Aug 21, 2024: This version adds new functions that process EurAm seed data 
# (29,410 x 1563 dataset), compute R^2_PGS as defined in our paper, and to plug 
# this quantity in as a substitute for r^2 for inference. This is to check that 
# the algorithm we run on PMBB phenotypes is also accurate in simulated phenotypes. 
# In LPC, this script is saved as PMBB_sim_master_script_082124.R (for purposes of
# tracking).
#
# (NEW) Oct 7, 2024: This version collapses the haplotype allele frequencies to 
# genotype allele frequencies. 
library(dplyr)
library(bigsnpr)
backup_dir <- "/project/mathilab/aaw/admix_prs/PMBB/Synthetic_Data/Eur_Genomes/backup/"

#' Simulate Unscaled Tagging Variant Effects (beta_j's)
#' 
#' Given the number of markers to generate effects,
#' expected squared correlation with phenotype, and the cross-population  
#' correlation of effect sizes between the two populations, this function
#' generates effect size vectors for both source and target populations. These
#' effects are not scaled down by allele frequency, to ensure same unscaled effects  
#' across the quantiles of PMBB cohort. Compare to simPopBetas in 
#' sim_master_script_012224.
#' 
# Example:
# 
# set.seed(2024)
# N_MARKERS <- 1000
# R2 <- 0.81
# RHO <- 0.4
# beta_vecs <-simPopBetas(n_markers=N_MARKERS,
#                         r2=R2,
#                         cross_pop_rho=RHO)
simPopBetas <- function(n_markers, 
                        r2,
                        cross_pop_rho) {
  # Half the r2 so that the genotypes will have desired squared correlation
  eff_r2 <- r2/2
  
  # Generate effect sizes for source and target populations
  beta_src <- rnorm(n=n_markers,mean=0,sd=sqrt(eff_r2))
  ind_gaussian <- rnorm(n=n_markers,mean=0,sd=1)
  beta_tar <- beta_src*cross_pop_rho + ind_gaussian*sqrt(eff_r2*(1-cross_pop_rho^2))
  
  # Return 
  return(list(SOURCE_POP_BETA_VEC=beta_src,
              TARGET_POP_BETA_VEC=beta_tar))
}

#' Generate EurAm demeaned genotypes and MAFs
#'
#' This is done once per seed, so that it can be used across all 
#' reps.
prepEurAmData <- function(bigsnpr_G) {
  # Get MAF
  maf <- colMeans(bigsnpr_G[1:29410,1:1563],na.rm=T)/2
  
  # Get demeaned genotypes
  copy_bigsnpr_G <- big_copy(bigsnpr_G,
                             type="double") # convert to double for big_increment
  neg_colmeans <- t(replicate(nrow(bigsnpr_G), 0-2*maf))
  big_increment(copy_bigsnpr_G,neg_colmeans)
  
  # Set NAs to 0
  big_apply(copy_bigsnpr_G, function(X,ind) {
    X.sub <- X[,ind,drop=FALSE]
    ind_na <- which(is.na(X.sub), arr.ind=TRUE)
    ind_na[, 2] <- ind[ind_na[, 2]]
    X[ind_na] <- 0
    NULL
  }, a.combine='c',block.size=200)
  
  # Return MAF, square-root of variances of MAFs and demeaned genotypes
  return(list(maf_vec = maf,
              G_demean = copy_bigsnpr_G,
              beta_denom = sqrt(sum(maf*(1-maf)))))
}

#' Compute PGS in EurAm cohort
#' 
#' Given betas and prepped_data (output of prepEurAmData), 
#' do the following:
#'  1. Construct EurAm beta_vec
#'  2. Construct PGS with EurAm beta_vec
#'  3. Simulate noise to obtain EurAm phenotypes 
#'  4. Compute R^2_PGS 
#'  5. To return: R^2_PGS (should be roughly squared correlation)
#'  
#'  Note to self: beta_vec is beta_vecs[["SOURCE_POP_BETA_VEC"]]
getEurAmPGS <- function(beta_vec,
                        r2,
                        prepped_data) {
  # Get scaled effect sizes
  spop_beta_vec_g <- beta_vec/prepped_data$beta_denom
  
  # Compute PGS 
  # take hat(X)*beta
  euram_pgs <- big_prodVec(X=prepped_data$G_demean,
                           y.col=spop_beta_vec_g)
  
  # Compute variance of PGS 
  emp_var <- var(euram_pgs)
  
  # Add noise term to get phenotype
  epsilons_theory <- rnorm(length(euram_pgs),mean=0,sd=sqrt(1-r2))
  epsilons_emp <- rnorm(length(euram_pgs),mean=0,sd=sqrt(1-emp_var))
  euram_phenos_theory <- euram_pgs + epsilons_theory
  euram_phenos_emp <- euram_pgs + epsilons_emp
  
  # Compute R^2_PGS
  r2_pgs_theory <- summary(lm(euram_phenos_theory~euram_pgs))$r.squared
  r2_pgs_emp <- summary(lm(euram_phenos_emp~euram_pgs))$r.squared
  
  # Return
  # first entry is relevant; second is just for tracking
  return(c(r2_pgs_theory,r2_pgs_emp)) 
} 

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
#' 
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

#' Compute admixed PGS (Version 2)
#' 
#' This function constructs the quantitative phenotype alongside various polygenic
#' scores (PGS). The phenotype is constructed based on the user specifying either 
#' a `model` (must be "global" or "local"; default setting is "local"), or both
#' a mixing parameter `lambda` and mixing model type `mix_type`, which together
#' specify the mixture model (between global and local). Effect sizes provided (
#' in `source_pop_beta_vec` and `target_pop_beta_vec`) and the haplotype matrices 
#' supplied (in `admixed_haps_list`) are used to construct the phenotype. 
#' Allelic dosages are demeaned by allele frequencies based on the 
#' labeled ancestry of the marker for the individual. Exogenous noise is simulated 
#' such that the total variance of the phenotype is roughly 1. (The parameters `r2` 
#' and `cross_pop_rho` are used to compute the genotypic variance.) 
#' 
#' PGSs are constructed using original haplotype matrices. We compute admPGS,
#' eurPGS and parPGS. These quantities do not depend on the model or the mixing
#' parameter and type. 
#' 
simPGS2 <- function(r2,
                    cross_pop_rho,
                    admixed_haps_list,
                    source_pop_beta_vec,
                    target_pop_beta_vec,
                    model,
                    lambda,
                    mix_type=NULL) {
  
  # Check for errors
  if (!is.null(model)) {
    assertthat::assert_that(model %in% c("global","local"),
                            msg="Model is neither global nor local.")
  } 
  if (!is.null(lambda)) {
    assertthat::assert_that(lambda < 1 & lambda > 0,
                            msg="Lambda is not strictly between 0 and 1")
    assertthat::assert_that(!is.null(mix_type),
                            msg="Mixture model type is not specified.")
  }
  if (!is.null(mix_type)) {
    assertthat::assert_that(mix_type %in% c("convex", "coinflip"),
                            msg="Mixture model type is neither convex nor coinflip.")
  }

  # Construct stacked effect size vectors to compute genotype-PGS
  big_source_pop_beta_vec <- source_pop_beta_vec
  big_target_pop_beta_vec <- target_pop_beta_vec
  
  # Compute masked matrices
  source_pop_masked_mat <- admixed_haps_list$SOURCE_POP_ANC_MAT*
    admixed_haps_list$HAP_MAT
  target_pop_masked_mat <- admixed_haps_list$TARGET_POP_ANC_MAT*
    admixed_haps_list$HAP_MAT
  
  # [!] >>>>>>>> Convert NAs to zeros for purposes of computing PGS <<<<<<<< [!]
  # Commented out on 10/9/2024
  #source_pop_masked_mat[is.na(source_pop_masked_mat)] <- 0
  #target_pop_masked_mat[is.na(target_pop_masked_mat)] <- 0
  
  # Extract target population ancestry fraction vectors
  marker_target_pop_frac <- admixed_haps_list$MARKER_TARGET_POP_FRAC
  ind_target_pop_frac <- admixed_haps_list$IND_TARGET_POP_FRAC
  n_ind <- length(ind_target_pop_frac)
  
  # Estimate the overall target population fraction of cohort
  cohort_est_target_pop_frac <- mean(marker_target_pop_frac)

  # Get 2 x p so that can halve later
  two_p <- ncol(source_pop_masked_mat)
  
  # Compute pop-specific allele frequencies 
  # [!] 10/8/2024 - use genotype allele frequencies rather than haplotypes
  source_pop_num <- colSums(source_pop_masked_mat,na.rm=TRUE)
  source_pop_num <- source_pop_num[1:(two_p/2)] + source_pop_num[(two_p/2+1):two_p]
  source_pop_denom <- colSums(admixed_haps_list$SOURCE_POP_ANC_MAT,na.rm=TRUE)
  source_pop_denom <- source_pop_denom[1:(two_p/2)] + 
    source_pop_denom[(two_p/2+1):two_p]
  source_pop_allele_freq <- source_pop_num/source_pop_denom
  source_pop_allele_freq <- c(source_pop_allele_freq,source_pop_allele_freq) 
  
  target_pop_num <- colSums(target_pop_masked_mat,na.rm=TRUE)
  target_pop_num <- target_pop_num[1:(two_p/2)] + target_pop_num[(two_p/2+1):two_p]
  target_pop_denom <- colSums(admixed_haps_list$TARGET_POP_ANC_MAT,na.rm=TRUE)
  target_pop_denom <- target_pop_denom[1:(two_p/2)] + 
    target_pop_denom[(two_p/2+1):two_p]
  target_pop_allele_freq <- target_pop_num/target_pop_denom
  target_pop_allele_freq <- c(target_pop_allele_freq,target_pop_allele_freq)
  
  # Demean masked matrices by previously computed allele frequencies
  source_pop_masked_mat_demean <- admixed_haps_list$SOURCE_POP_ANC_MAT*
    sweep(admixed_haps_list$HAP_MAT,2,source_pop_allele_freq)
  target_pop_masked_mat_demean <- admixed_haps_list$TARGET_POP_ANC_MAT*
    sweep(admixed_haps_list$HAP_MAT,2,target_pop_allele_freq)
  
  # [!] >>>>>>>> Convert NAs to zeros for purposes of computing PGS <<<<<<<< [!]
  source_pop_masked_mat_demean[is.na(source_pop_masked_mat_demean)] <- 0
  target_pop_masked_mat_demean[is.na(target_pop_masked_mat_demean)] <- 0
  
  # Construct PGSs (demean and original)
  # admPGS (original)
  X_beta_source_pop <- as.numeric( # this is parPGS (original)
    source_pop_masked_mat%*%big_source_pop_beta_vec) 
  X_beta_target_pop <- as.numeric(
    target_pop_masked_mat%*%big_target_pop_beta_vec)
  admPGS <- X_beta_source_pop + X_beta_target_pop
  
  # eurPGS (original)
  eurPGS <- X_beta_source_pop + as.numeric(
    target_pop_masked_mat%*%big_source_pop_beta_vec)
  
  # admPGS (demean)
  parPGS_demean <- as.numeric( # this is parPGS (demean)
    source_pop_masked_mat_demean%*%big_source_pop_beta_vec)
  admPGS_demean <- parPGS_demean + 
    as.numeric(target_pop_masked_mat_demean%*%big_target_pop_beta_vec)
  
  # eurPGS (demean)
  eurPGS_demean <- parPGS_demean +
    as.numeric(target_pop_masked_mat_demean%*%big_source_pop_beta_vec)
  
  if (!is.null(model)) {
    # Construct based on specified model
    if (model=="local") {
      # (Demeaned X)*b is used to construct y
      X_demean_beta_source_pop <- as.numeric(
        source_pop_masked_mat_demean%*%big_source_pop_beta_vec)
      X_demean_beta_target_pop <- as.numeric(
        target_pop_masked_mat_demean%*%big_target_pop_beta_vec)
      X_demean_beta <- X_demean_beta_source_pop + X_demean_beta_target_pop
      
      # Compute noise vector
      adm_noise_vec <- rnorm(n=nrow(admixed_haps_list$HAP_MAT),
                             mean=0,sd=sqrt(1-r2))
      # Compute y
      y <- X_demean_beta+adm_noise_vec
      
    } else {
      # Compute individual effect sizes based on global ancestry
      mixed_beta_mat <- outer((1-ind_target_pop_frac),big_source_pop_beta_vec) + 
        outer(ind_target_pop_frac,big_target_pop_beta_vec)
      
      # (Demeaned X)*b is used to construct y
      X_demean_beta_source_pop <- as.numeric(rowSums(
        source_pop_masked_mat_demean*mixed_beta_mat, na.rm=TRUE)) 
      X_demean_beta_target_pop <- as.numeric(rowSums(
        target_pop_masked_mat_demean*mixed_beta_mat, na.rm=TRUE))
      X_demean_beta <- X_demean_beta_source_pop + X_demean_beta_target_pop
      
      # Compute analytical quantities for genotypic variance under global model
      target_pop_anc_mat <- admixed_haps_list$TARGET_POP_ANC_MAT
      source_pop_anc_mat <- admixed_haps_list$SOURCE_POP_ANC_MAT
      
      source_pop_snp_var <- source_pop_allele_freq*(1-source_pop_allele_freq)
      target_pop_snp_var <- target_pop_allele_freq*(1-target_pop_allele_freq)
      
      inner_term_1 <- (1-ind_target_pop_frac)^2*r2/sum(source_pop_snp_var) # vector term 1
      inner_term_2 <- ind_target_pop_frac^2*r2/sum(target_pop_snp_var) # vector term 2
      inner_term_3 <- 2*cross_pop_rho*r2*
        (1-ind_target_pop_frac)*ind_target_pop_frac/
        sqrt(sum(source_pop_snp_var)*sum(target_pop_snp_var)) # vector term 3
      inner_term <- inner_term_1 + inner_term_2 + inner_term_3 # inner term (indexed by i)
      
      inner_term_mat <- matrix(inner_term,
                               nrow=length(inner_term),
                               ncol=length(source_pop_snp_var))
      source_pop_snp_var_mat <- matrix(source_pop_snp_var,
                                       nrow=n_ind,
                                       ncol=length(source_pop_snp_var),
                                       byrow=TRUE)
      target_pop_snp_var_mat <- matrix(target_pop_snp_var,
                                       nrow=n_ind,
                                       ncol=length(target_pop_snp_var),
                                       byrow=TRUE)
      
      total_var_1 <- colSums(source_pop_anc_mat*
                               source_pop_snp_var_mat*inner_term_mat, na.rm=TRUE)
      total_var_2 <- colSums(target_pop_anc_mat*
                               target_pop_snp_var_mat*inner_term_mat, na.rm=TRUE)
      total_var <- sum(total_var_1 + total_var_2)/n_ind
      
      # Compute noise vector for global model using analytical quantities
      adm_noise_vec <- rnorm(n=nrow(admixed_haps_list$HAP_MAT),
                             mean=0,sd=sqrt(1-total_var))
      # Compute y
      y <- X_demean_beta + adm_noise_vec
    }
  } else {
    # Construct based on mixture model
    message(date(),": Generating mixture model with lambda = ", lambda)
    target_pop_anc_mat <- admixed_haps_list$TARGET_POP_ANC_MAT
    source_pop_anc_mat <- admixed_haps_list$SOURCE_POP_ANC_MAT
    
    # Compute effects under the local and global model
    # Global effect size matrix
    global_beta_mat <- outer((1-ind_target_pop_frac),big_source_pop_beta_vec) + 
      outer(ind_target_pop_frac,big_target_pop_beta_vec)
    
    big_source_pop_beta_mat <- matrix(big_source_pop_beta_vec,
                                      nrow=n_ind,
                                      ncol=length(big_source_pop_beta_vec),
                                      byrow=TRUE)
    big_target_pop_beta_mat <- matrix(big_target_pop_beta_vec,
                                      nrow=n_ind,
                                      ncol=length(big_target_pop_beta_vec),
                                      byrow=TRUE)
    local_beta_mat <- source_pop_anc_mat*big_source_pop_beta_mat + 
      target_pop_anc_mat*big_target_pop_beta_mat
    
    # [!] >>>>>>> Convert NAs to zeros for purposes of computing PGS <<<<<<< [!]
    local_beta_mat[is.na(local_beta_mat)] <- 0
    
    if (mix_type=="convex") {
      message(date(),": Generating mixed effects using convex combination approach")
      combined_beta_mat <- lambda*local_beta_mat + (1-lambda)*global_beta_mat
      
    } else {
      message(date(),": Generating mixed effects using coin flip approach")
      local_model_markers <- rbinom(n=length(big_target_pop_beta_vec),
                                    size=1,prob=lambda)
      global_model_markers <- 1-local_model_markers
      local_model_marker_mat <- matrix(local_model_markers,
                                       nrow=n_ind,
                                       ncol=length(local_model_markers),
                                       byrow=TRUE)
      global_model_marker_mat <- matrix(global_model_markers,
                                        nrow=n_ind,
                                        ncol=length(global_model_markers),
                                        byrow=TRUE)
      combined_beta_mat <- local_beta_mat*local_model_marker_mat +
        global_beta_mat*global_model_marker_mat
      
      # [!] >>>>>> Convert NAs to zeros for purposes of computing PGS <<<<<< [!]
      combined_beta_mat[is.na(combined_beta_mat)] <- 0 
    }
    
    # Construct based on mixture model specified
    X_demean_beta_source_pop <- as.numeric(rowSums(source_pop_masked_mat_demean*
                                                     combined_beta_mat,
                                                   na.rm=TRUE)) 
    X_demean_beta_target_pop <- as.numeric(rowSums(target_pop_masked_mat_demean*
                                                     combined_beta_mat,
                                                   na.rm=TRUE))
    X_demean_beta <- X_demean_beta_source_pop + X_demean_beta_target_pop
    
    # Compute noise vector for mixture model
    total_var <- var(X_demean_beta)
    adm_noise_vec <- rnorm(n=nrow(admixed_haps_list$HAP_MAT),
                           mean=0,sd=sqrt(1-total_var))
    # Compute y
    y <- X_demean_beta+adm_noise_vec
  }
  
  # Return
  return(list(Y=y,
              admPGS_orig=admPGS,
              eurPGS_orig=eurPGS,
              parPGS_orig=X_beta_source_pop,
              admPGS_demean=admPGS_demean,
              eurPGS_demean=eurPGS_demean,
              parPGS_demean=parPGS_demean,
              est_target_pop_frac=cohort_est_target_pop_frac))
}

#' Compute analytical expectations
#' 
#' This function computes the analytical expectations of PGS metrics using the
#' problem parameters, which include (1) allele frequencies of each population 
#' (`source_pop_freq_vec` and `target_pop_freq_vec`); (2) squared correlation and
#' correlation of effect sizes (`r2` and `cross_pop_rho`); (3) target population
#' ancestry fractions at each marker (`marker_target_pop_frac`).
#' If under the global model or a mixture model (i.e., `lambda` is not 0 or 1), 
#' the target population ancestry fraction per individual (`ind_target_pop_frac`)
#' is also required, in addition to the number of individuals (`n_ind`), as well as
#' the ancestry label matrices (`source_pop_anc_mat` and `target_pop_anc_mat`). 
#' The following expectations of quantities are returned:
#' (1) squared correlation of parPGS and y
#' (2) squared correlation of eurPGS and y
#' (3) variance of parPGS
#' (4) variance of eurPGS
#'   
computeExpected <- function(r2,
                            cross_pop_rho,
                            admixed_haps_list,
                            model,
                            lambda=NULL) {
  # Check for errors
  if (!is.null(model)) {
    assertthat::assert_that(model %in% c("global","local"),
                            msg="Model is neither global nor local.")
  } 
  if (!is.null(lambda)) {
    assertthat::assert_that(lambda < 1 & lambda > 0,
                            msg="Lambda is not strictly between 0 and 1")
    assertthat::assert_that(is.null(model),
                            msg="Model should be set to NULL if using mixture model")
  }

  # Define key quantities required in rest of calculations
  # - demeaned masked matrices
  # - ancestry label matrices
  # - ancestry-specific allele frequencies
  # - marker-specific target pop frac, individual-specific global ancestry
  # - snp variation (f_j(1-f_j))
  
  source_pop_masked_demean_mat <- admixed_haps_list$SOURCE_POP_DEMEAN_HAP_MAT
  target_pop_masked_demean_mat <- admixed_haps_list$TARGET_POP_DEMEAN_HAP_MAT
  
  source_pop_anc_mat <- admixed_haps_list$SOURCE_POP_ANC_MAT
  target_pop_anc_mat <- admixed_haps_list$TARGET_POP_ANC_MAT
  
  marker_target_pop_frac <- admixed_haps_list$MARKER_TARGET_POP_FRAC
  ind_target_pop_frac <- admixed_haps_list$IND_TARGET_POP_FRAC
  
  source_pop_allele_freq <- admixed_haps_list$SOURCE_POP_AF
  target_pop_allele_freq <- admixed_haps_list$TARGET_POP_AF
  
  source_pop_snp_var <- source_pop_allele_freq*(1-source_pop_allele_freq) # f_j^s(1-f_j^s), assumed length 2p 
  target_pop_snp_var <- target_pop_allele_freq*(1-target_pop_allele_freq) # f_j^t(1-f_j^t), assumed length 2p 
  
  # Split all quantities into two haplotypes
  n_ind <- admixed_haps_list[["N_IND"]]
  n_markers <- length(source_pop_snp_var)/2
  
  # Check that n_markers computed above is the same as N_MARKERS computed from
  # admixed_haps_list
  assertthat::assert_that(n_markers==admixed_haps_list$N_MARKERS,
                          msg="n_markers computed from source_pop_snp_var disagrees with N_MARKERS.")
  
  h1_source_pop_snp_var <- source_pop_snp_var[1:n_markers]
  h2_source_pop_snp_var <- source_pop_snp_var[(n_markers+1):(2*n_markers)]
  h1_target_pop_snp_var <- target_pop_snp_var[1:n_markers]
  h2_target_pop_snp_var <- target_pop_snp_var[(n_markers+1):(2*n_markers)]
  
  h1_marker_target_pop_frac <- marker_target_pop_frac[1:n_markers]
  h2_marker_target_pop_frac <- marker_target_pop_frac[(n_markers+1):(2*n_markers)]
  
  h1_source_pop_masked_demean_mat <- source_pop_masked_demean_mat[,1:n_markers]
  h2_source_pop_masked_demean_mat <- source_pop_masked_demean_mat[,(n_markers+1):(2*n_markers)]
  h1_target_pop_masked_demean_mat <- target_pop_masked_demean_mat[,1:n_markers]
  h2_target_pop_masked_demean_mat <- target_pop_masked_demean_mat[,(n_markers+1):(2*n_markers)]
  
  # E[Var(parPGS)]
  E_var_parPGS_term_1 <- 0.5*
    sum((1-h1_marker_target_pop_frac)*h1_source_pop_snp_var)/
    sum(h1_source_pop_snp_var)
  E_var_parPGS_term_2 <- 0.5*
    sum((1-h2_marker_target_pop_frac)*h2_source_pop_snp_var)/
    sum(h2_source_pop_snp_var)
  E_var_parPGS_term_3 <- sum(colMeans(h1_source_pop_masked_demean_mat*
                                        h2_source_pop_masked_demean_mat))/
    sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
  E_var_parPGS <- r2*(E_var_parPGS_term_1+E_var_parPGS_term_2+E_var_parPGS_term_3)
  
  # E[Var(eurPGS)]
  E_var_eurPGS_term_1 <- E_var_parPGS_term_1 +
    0.5*sum(h1_marker_target_pop_frac*h1_target_pop_snp_var)/
    sum(h1_source_pop_snp_var)
  E_var_eurPGS_term_2 <- E_var_parPGS_term_2 +
    0.5*sum(h2_marker_target_pop_frac*h2_target_pop_snp_var)/
    sum(h2_source_pop_snp_var)
  h1_demean_mat <- h1_source_pop_masked_demean_mat+h1_target_pop_masked_demean_mat
  h2_demean_mat <- h2_source_pop_masked_demean_mat+h2_target_pop_masked_demean_mat
  E_var_eurPGS_term_3 <- sum(colMeans(h1_demean_mat*h2_demean_mat))/
    sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
  E_var_eurPGS <- r2*(E_var_eurPGS_term_1+E_var_eurPGS_term_2+E_var_eurPGS_term_3)
  
  # [!] >>>>>> 3/22/24 <<<<<<
  # Calculate expected quantities
  if (identical(model,"local")) {
    # Compute local model expressions
    # E[Cov(eurPGS,y)]
    same_anc_terms <- E_var_parPGS_term_1+
      E_var_parPGS_term_2+
      0.5*(sum(colMeans(h1_demean_mat*h2_source_pop_masked_demean_mat))+ 
             sum(colMeans(h2_demean_mat*h1_source_pop_masked_demean_mat)))/
      sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
    
    cross_anc_term_1 <- 0.5*sum(h1_marker_target_pop_frac*h1_target_pop_snp_var)/
      sqrt(sum(h1_source_pop_snp_var)*sum(h1_target_pop_snp_var))
    cross_anc_term_2 <- 0.5*sum(h2_marker_target_pop_frac*h2_target_pop_snp_var)/
      sqrt(sum(h2_source_pop_snp_var)*sum(h2_target_pop_snp_var))
    cross_anc_term_3 <- 0.5*sum(colMeans(h2_target_pop_masked_demean_mat*
                                           h1_demean_mat))/
      sqrt(sum(h1_source_pop_snp_var)*sum(h2_target_pop_snp_var))
    cross_anc_term_4 <- 0.5*sum(colMeans(h1_target_pop_masked_demean_mat*
                                           h2_demean_mat))/
      sqrt(sum(h1_target_pop_snp_var)*sum(h2_source_pop_snp_var))
    cross_anc_terms <- cross_anc_term_1+cross_anc_term_2+
      cross_anc_term_3+cross_anc_term_4
    E_cov_eurPGS_y <- r2*same_anc_terms+r2*cross_pop_rho*cross_anc_terms
    E_cor2_eurPGS_y <- E_cov_eurPGS_y^2/E_var_eurPGS
    
    # E[Cov(parPGS,y)]
    same_anc_terms <- E_var_parPGS_term_1+
      E_var_parPGS_term_2+
      sum(colMeans(h1_source_pop_masked_demean_mat*h2_source_pop_masked_demean_mat))/
      sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
    cross_anc_term_1 <- 0.5*sum(colMeans(h1_source_pop_masked_demean_mat*h2_target_pop_masked_demean_mat))/
      sqrt(sum(h1_source_pop_snp_var)*sum(h2_target_pop_snp_var))
    cross_anc_term_2 <- 0.5*sum(colMeans(h2_source_pop_masked_demean_mat*h1_target_pop_masked_demean_mat))/
      sqrt(sum(h2_source_pop_snp_var)*sum(h1_target_pop_snp_var))
    cross_anc_terms <- cross_anc_term_1+cross_anc_term_2
    E_cov_parPGS_y <- r2*same_anc_terms+r2*cross_pop_rho*cross_anc_terms
    E_cor2_parPGS_y <- E_cov_parPGS_y^2/E_var_parPGS  
    
  } else {
    # Compute global model expressions
    # Calculate helpful matrices -- individual-level global ancestry matrices 
    # These matrices are n x p rather than n x 2p, to facilitate their use 
    # in computing haplotype-specific quantities (since they don't differ)
    source_pop_global_anc_mat <- matrix(1-ind_target_pop_frac,
                                        nrow=n_ind,
                                        ncol=n_markers,
                                        byrow=FALSE)
    target_pop_global_anc_mat <- matrix(ind_target_pop_frac,
                                        nrow=n_ind,
                                        ncol=n_markers,
                                        byrow=FALSE)
    # eurPGS
    same_anc_term_1 <- 0.5*sum(colMeans(source_pop_global_anc_mat*h1_demean_mat*h1_demean_mat))/
      sum(h1_source_pop_snp_var)
    same_anc_term_2 <- 0.5*sum(colMeans(source_pop_global_anc_mat*h2_demean_mat*h2_demean_mat))/
      sum(h2_source_pop_snp_var)
    same_anc_term_3 <- sum(colMeans(source_pop_global_anc_mat*h1_demean_mat*h2_demean_mat))/
      sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
    same_anc_terms <- same_anc_term_1+same_anc_term_2+same_anc_term_3
    
    cross_anc_term_1 <- 0.5*sum(colMeans(target_pop_global_anc_mat*h1_demean_mat*h1_demean_mat))/
      sqrt(sum(h1_source_pop_snp_var)*sum(h1_target_pop_snp_var))  
    cross_anc_term_2 <- 0.5*sum(colMeans(target_pop_global_anc_mat*h2_demean_mat*h2_demean_mat))/
      sqrt(sum(h2_source_pop_snp_var)*sum(h2_target_pop_snp_var))  
    cross_anc_term_3 <- 0.5*sum(colMeans(target_pop_global_anc_mat*h1_demean_mat*h2_demean_mat))/
      sqrt(sum(h2_source_pop_snp_var)*sum(h1_target_pop_snp_var))
    cross_anc_term_4 <- 0.5*sum(colMeans(target_pop_global_anc_mat*h1_demean_mat*h2_demean_mat))/
      sqrt(sum(h1_source_pop_snp_var)*sum(h2_target_pop_snp_var))
    cross_anc_terms <- cross_anc_term_1+cross_anc_term_2+cross_anc_term_3+cross_anc_term_4
    
    E_cov_eurPGS_y <- r2*same_anc_terms+r2*cross_pop_rho*cross_anc_terms
    E_cor2_eurPGS_y <- E_cov_eurPGS_y^2/E_var_eurPGS
    
    # parPGS
    same_anc_term_1 <- 0.5*sum(colMeans(source_pop_global_anc_mat*
                                      h1_source_pop_masked_demean_mat*h2_demean_mat))/
      sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
    same_anc_term_2 <- 0.5*sum(colMeans(source_pop_global_anc_mat*
                                         h2_source_pop_masked_demean_mat*h1_demean_mat))/
      sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
    same_anc_term_3 <- 0.5*sum(colMeans(source_pop_global_anc_mat*
                                           h1_source_pop_masked_demean_mat*h1_demean_mat))/
      sum(h1_source_pop_snp_var)
    same_anc_term_4 <- 0.5*sum(colMeans(source_pop_global_anc_mat*
                                          h2_source_pop_masked_demean_mat*h2_demean_mat))/
      sum(h2_source_pop_snp_var)
    same_anc_terms <- same_anc_term_1+same_anc_term_2+same_anc_term_3+same_anc_term_4 
    
    cross_anc_term_1 <- 0.5*sum(colMeans(target_pop_global_anc_mat*
                                           h1_source_pop_masked_demean_mat*h1_demean_mat))/
      sqrt(sum(h1_source_pop_snp_var)*sum(h1_target_pop_snp_var))
    cross_anc_term_2 <- 0.5*sum(colMeans(target_pop_global_anc_mat*
                                           h2_source_pop_masked_demean_mat*h2_demean_mat))/
      sqrt(sum(h2_source_pop_snp_var)*sum(h2_target_pop_snp_var))
    cross_anc_term_3 <- 0.5*sum(colMeans(target_pop_global_anc_mat*
                                           h2_source_pop_masked_demean_mat*h1_demean_mat))/
      sqrt(sum(h1_target_pop_snp_var)*sum(h2_source_pop_snp_var))
    cross_anc_term_4 <- 0.5*sum(colMeans(target_pop_global_anc_mat*
                                           h1_source_pop_masked_demean_mat*h2_demean_mat))/
      sqrt(sum(h2_target_pop_snp_var)*sum(h1_source_pop_snp_var))
    cross_anc_terms <- cross_anc_term_1+cross_anc_term_2+cross_anc_term_3+cross_anc_term_4
    
    E_cov_parPGS_y <- r2*same_anc_terms+r2*cross_pop_rho*cross_anc_terms
    E_cor2_parPGS_y <- E_cov_parPGS_y^2/E_var_parPGS
    
    if (is.null(model)) {
      # Compute local model expressions (to combine with global model)
      # E[Cov(eurPGS,y)]
      same_anc_terms <- E_var_parPGS_term_1+
        E_var_parPGS_term_2+
        0.5*(sum(colMeans(h1_demean_mat*h2_source_pop_masked_demean_mat))+ 
               sum(colMeans(h2_demean_mat*h1_source_pop_masked_demean_mat)))/
        sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
      
      cross_anc_term_1 <- 0.5*sum(h1_marker_target_pop_frac*h1_target_pop_snp_var)/
        sqrt(sum(h1_source_pop_snp_var)*sum(h1_target_pop_snp_var))
      cross_anc_term_2 <- 0.5*sum(h2_marker_target_pop_frac*h2_target_pop_snp_var)/
        sqrt(sum(h2_source_pop_snp_var)*sum(h2_target_pop_snp_var))
      cross_anc_term_3 <- 0.5*sum(colMeans(h2_target_pop_masked_demean_mat*
                                             h1_demean_mat))/
        sqrt(sum(h1_source_pop_snp_var)*sum(h2_target_pop_snp_var))
      cross_anc_term_4 <- 0.5*sum(colMeans(h1_target_pop_masked_demean_mat*
                                             h2_demean_mat))/
        sqrt(sum(h1_target_pop_snp_var)*sum(h2_source_pop_snp_var))
      cross_anc_terms <- cross_anc_term_1+cross_anc_term_2+
        cross_anc_term_3+cross_anc_term_4
      E_cov_eurPGS_y <- r2*same_anc_terms+r2*cross_pop_rho*cross_anc_terms
      E_cor2_eurPGS_y_loc <- E_cov_eurPGS_y^2/E_var_eurPGS
      
      # E[Cov(parPGS,y)]
      same_anc_terms <- E_var_parPGS_term_1+
        E_var_parPGS_term_2+
        sum(colMeans(h1_source_pop_masked_demean_mat*h2_source_pop_masked_demean_mat))/
        sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
      cross_anc_term_1 <- 0.5*sum(colMeans(h1_source_pop_masked_demean_mat*h2_target_pop_masked_demean_mat))/
        sqrt(sum(h1_source_pop_snp_var)*sum(h2_target_pop_snp_var))
      cross_anc_term_2 <- 0.5*sum(colMeans(h2_source_pop_masked_demean_mat*h1_target_pop_masked_demean_mat))/
        sqrt(sum(h2_source_pop_snp_var)*sum(h1_target_pop_snp_var))
      cross_anc_terms <- cross_anc_term_1+cross_anc_term_2
      E_cov_parPGS_y <- r2*same_anc_terms+r2*cross_pop_rho*cross_anc_terms
      E_cor2_parPGS_y_loc <- E_cov_parPGS_y^2/E_var_parPGS  
      
      
      # Take convex combination
      E_cor2_eurPGS_y <- (lambda*sqrt(E_cor2_eurPGS_y_loc)+
                          (1-lambda)*sqrt(E_cor2_eurPGS_y))^2
      E_cor2_parPGS_y <- (lambda*sqrt(E_cor2_parPGS_y_loc)+
                          (1-lambda)*sqrt(E_cor2_parPGS_y))^2
    }
  }
  
  # Return
  return(data.frame(cor2_y_eurPGS=E_cor2_eurPGS_y,
                    cor2_y_parPGS=E_cor2_parPGS_y,
                    var_eurPGS=E_var_eurPGS,
                    var_parPGS=E_var_parPGS))
}

#' Inference of model parameters
#' 
#' Given empirical squared correlations and variances, as well as estimated
#' target population ancestry fractions and population-specific allele frequencies,
#' this functions performs a sequential inference procedure. First it estimates 
#' the unsigned correlation between the admPGS and phenotype (`r`) as well as the
#' global correlation of effects (`rho`), using the empirical correlation of 
#' phenotype with the eurPGS. Next, it estimates the parPGS correlations under 
#' Global and Local models. These are subsequently used to estimate the mixture
#' parameter (`lambda`).
#' 
inferParams <- function(emp_cor2_y_eurPGS,
                        emp_cor2_y_parPGS,
                        emp_var_eurPGS,
                        emp_var_parPGS,
                        admixed_haps_list,
                        snp_herit=0.8,
                        ignore_var=TRUE,
                        R2_PGS=c(NA,NA),
                        use_R2_PGS=FALSE) {
  # Define number of individuals and number of markers first
  n_ind <- admixed_haps_list[["N_IND"]]
  n_markers <- admixed_haps_list[["N_MARKERS"]]
  
  # Step 1: Infer r and r*rho
  message("Step 1: Inferring r")
  
  # Extract quantities from admixed_haps_list
  source_pop_masked_demean_mat <- admixed_haps_list$SOURCE_POP_DEMEAN_HAP_MAT
  target_pop_masked_demean_mat <- admixed_haps_list$TARGET_POP_DEMEAN_HAP_MAT
  source_pop_anc_mat <- admixed_haps_list$SOURCE_POP_ANC_MAT
  target_pop_anc_mat <- admixed_haps_list$TARGET_POP_ANC_MAT
  marker_target_pop_frac <- admixed_haps_list$MARKER_TARGET_POP_FRAC
  ind_target_pop_frac <- admixed_haps_list$IND_TARGET_POP_FRAC
  source_pop_allele_freq <- admixed_haps_list$SOURCE_POP_AF
  target_pop_allele_freq <- admixed_haps_list$TARGET_POP_AF
  source_pop_snp_var <- source_pop_allele_freq*(1-source_pop_allele_freq) 
  target_pop_snp_var <- target_pop_allele_freq*(1-target_pop_allele_freq) 
  
  # Extract h1 and h2 versions from the stacked version
  h1_source_pop_snp_var <- source_pop_snp_var[1:n_markers]
  h2_source_pop_snp_var <- source_pop_snp_var[(n_markers+1):(2*n_markers)]
  h1_target_pop_snp_var <- target_pop_snp_var[1:n_markers]
  h2_target_pop_snp_var <- target_pop_snp_var[(n_markers+1):(2*n_markers)]
  
  h1_marker_target_pop_frac <- marker_target_pop_frac[1:n_markers]
  h2_marker_target_pop_frac <- marker_target_pop_frac[(n_markers+1):(2*n_markers)]
  
  h1_source_pop_masked_demean_mat <- source_pop_masked_demean_mat[,1:n_markers]
  h2_source_pop_masked_demean_mat <- source_pop_masked_demean_mat[,(n_markers+1):(2*n_markers)]
  h1_target_pop_masked_demean_mat <- target_pop_masked_demean_mat[,1:n_markers]
  h2_target_pop_masked_demean_mat <- target_pop_masked_demean_mat[,(n_markers+1):(2*n_markers)]
  
  # Estimate hat_r2_par from E[Var(parPGS)] equation
  E_var_parPGS_term_1 <- 0.5*
    sum((1-h1_marker_target_pop_frac)*h1_source_pop_snp_var)/
    sum(h1_source_pop_snp_var)
  E_var_parPGS_term_2 <- 0.5*
    sum((1-h2_marker_target_pop_frac)*h2_source_pop_snp_var)/
    sum(h2_source_pop_snp_var)
  E_var_parPGS_term_3 <- sum(colMeans(h1_source_pop_masked_demean_mat*
                                        h2_source_pop_masked_demean_mat))/
    sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
  hat_r2_par <- emp_var_parPGS/
    (E_var_parPGS_term_1+E_var_parPGS_term_2+E_var_parPGS_term_3)
  
  # Estimate hat_r2_eur from E[Var(eurPGS)] equation
  E_var_eurPGS_term_1 <- E_var_parPGS_term_1 +
    0.5*sum(h1_marker_target_pop_frac*h1_target_pop_snp_var)/
    sum(h1_source_pop_snp_var)
  E_var_eurPGS_term_2 <- E_var_parPGS_term_2 +
    0.5*sum(h2_marker_target_pop_frac*h2_target_pop_snp_var)/
    sum(h2_source_pop_snp_var)
  h1_demean_mat <- h1_source_pop_masked_demean_mat+h1_target_pop_masked_demean_mat
  h2_demean_mat <- h2_source_pop_masked_demean_mat+h2_target_pop_masked_demean_mat
  E_var_eurPGS_term_3 <- sum(colMeans(h1_demean_mat*h2_demean_mat))/
    sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
  hat_r2_eur <- emp_var_eurPGS/
    (E_var_eurPGS_term_1+E_var_eurPGS_term_2+E_var_eurPGS_term_3)
  
  # Compute r_hat by taking root mean square 
  if (ignore_var) {
    message("Ignoring var(parPGS) and var(eurPGS), using SNP heritability")
    if (use_R2_PGS) {
      message("SNP heritability computed on EurAm simulated phenotype")
      r_hat <- sqrt(R2_PGS[1])
    } else {
      message("SNP heritability pre-computed (typically from an actual phenotype)")
      r_hat <- sqrt(snp_herit)
    }
  } else {
    r_hat <- sqrt(median(c(hat_r2_par,hat_r2_eur)))
  }
  
  # r2 estimates are probably unbiased, but sqrt(r2) will be biased downward, by Jensen's inequality 
  message("hat_r2_par = ", hat_r2_par)
  message("hat_r2_eur = ", hat_r2_eur)
  message("r_hat = ", r_hat)
  
  # Step 2: Infer rho
  message("Step 2: Inferring rho")
  
  # Compute terms from E[corr_eurPGS_y] under Local Model
  # c0 numerator
  c0_num_term <- E_var_parPGS_term_1+
    E_var_parPGS_term_2+
    0.5*(sum(colMeans(h1_demean_mat*h2_source_pop_masked_demean_mat))+ 
           sum(colMeans(h2_demean_mat*h1_source_pop_masked_demean_mat)))/
    sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
  
  cross_anc_term_1 <- 0.5*sum(h1_marker_target_pop_frac*h1_target_pop_snp_var)/
    sqrt(sum(h1_source_pop_snp_var)*sum(h1_target_pop_snp_var))
  cross_anc_term_2 <- 0.5*sum(h2_marker_target_pop_frac*h2_target_pop_snp_var)/
    sqrt(sum(h2_source_pop_snp_var)*sum(h2_target_pop_snp_var))
  cross_anc_term_3 <- 0.5*sum(colMeans(h2_target_pop_masked_demean_mat*
                                         h1_demean_mat))/
    sqrt(sum(h1_source_pop_snp_var)*sum(h2_target_pop_snp_var))
  cross_anc_term_4 <- 0.5*sum(colMeans(h1_target_pop_masked_demean_mat*
                                         h2_demean_mat))/
    sqrt(sum(h1_target_pop_snp_var)*sum(h2_source_pop_snp_var))
  # c1 numerator
  c1_num_term <- cross_anc_term_1+cross_anc_term_2+
    cross_anc_term_3+cross_anc_term_4
  
  denom_term <- sqrt(E_var_eurPGS_term_1+E_var_eurPGS_term_2+E_var_eurPGS_term_3)
  c0 <- c0_num_term/denom_term
  c1 <- c1_num_term/denom_term
  
  # Compute rho_hat
  rho_hat <- (sqrt(emp_cor2_y_eurPGS)/r_hat - c0)/c1
  message("rho_hat = ", rho_hat)
  
  # Step 3: Infer r_G and r_L
  message("Step 3: Inferring r_G and r_L")
  
  # Compute local model expected correlation (E^local[cor_parPGS_y])
  # E[Cov(parPGS,y)]
  same_anc_terms <- E_var_parPGS_term_1+
    E_var_parPGS_term_2+
    sum(colMeans(h1_source_pop_masked_demean_mat*h2_source_pop_masked_demean_mat))/
    sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
  cross_anc_term_1 <- 0.5*sum(colMeans(h1_source_pop_masked_demean_mat*h2_target_pop_masked_demean_mat))/
    sqrt(sum(h1_source_pop_snp_var)*sum(h2_target_pop_snp_var))
  cross_anc_term_2 <- 0.5*sum(colMeans(h2_source_pop_masked_demean_mat*h1_target_pop_masked_demean_mat))/
    sqrt(sum(h2_source_pop_snp_var)*sum(h1_target_pop_snp_var))
  cross_anc_terms <- cross_anc_term_1+cross_anc_term_2
  E_cov_parPGS_y <- r_hat^2*same_anc_terms+r_hat^2*rho_hat*cross_anc_terms
  rL_hat <- E_cov_parPGS_y/
    sqrt(r_hat^2*(E_var_parPGS_term_1+E_var_parPGS_term_2+E_var_parPGS_term_3))
  
  # Compute global model expected correlation (E^global[cor_parPGS_y])
  source_pop_global_anc_mat <- matrix(1-ind_target_pop_frac,
                                      nrow=n_ind,
                                      ncol=n_markers,
                                      byrow=FALSE)
  target_pop_global_anc_mat <- matrix(ind_target_pop_frac,
                                      nrow=n_ind,
                                      ncol=n_markers,
                                      byrow=FALSE)
  # E[Cov(parPGS,y)]
  same_anc_term_1 <- 0.5*sum(colMeans(source_pop_global_anc_mat*
                                        h1_source_pop_masked_demean_mat*h2_demean_mat))/
    sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
  same_anc_term_2 <- 0.5*sum(colMeans(source_pop_global_anc_mat*
                                        h2_source_pop_masked_demean_mat*h1_demean_mat))/
    sqrt(sum(h1_source_pop_snp_var)*sum(h2_source_pop_snp_var))
  same_anc_term_3 <- 0.5*sum(colMeans(source_pop_global_anc_mat*
                                        h1_source_pop_masked_demean_mat*h1_demean_mat))/
    sum(h1_source_pop_snp_var)
  same_anc_term_4 <- 0.5*sum(colMeans(source_pop_global_anc_mat*
                                        h2_source_pop_masked_demean_mat*h2_demean_mat))/
    sum(h2_source_pop_snp_var)
  same_anc_terms <- same_anc_term_1+same_anc_term_2+same_anc_term_3+same_anc_term_4 
  
  cross_anc_term_1 <- 0.5*sum(colMeans(target_pop_global_anc_mat*
                                         h1_source_pop_masked_demean_mat*h1_demean_mat))/
    sqrt(sum(h1_source_pop_snp_var)*sum(h1_target_pop_snp_var))
  cross_anc_term_2 <- 0.5*sum(colMeans(target_pop_global_anc_mat*
                                         h2_source_pop_masked_demean_mat*h2_demean_mat))/
    sqrt(sum(h2_source_pop_snp_var)*sum(h2_target_pop_snp_var))
  cross_anc_term_3 <- 0.5*sum(colMeans(target_pop_global_anc_mat*
                                         h2_source_pop_masked_demean_mat*h1_demean_mat))/
    sqrt(sum(h1_target_pop_snp_var)*sum(h2_source_pop_snp_var))
  cross_anc_term_4 <- 0.5*sum(colMeans(target_pop_global_anc_mat*
                                         h1_source_pop_masked_demean_mat*h2_demean_mat))/
    sqrt(sum(h2_target_pop_snp_var)*sum(h1_source_pop_snp_var))
  cross_anc_terms <- cross_anc_term_1+cross_anc_term_2+cross_anc_term_3+cross_anc_term_4
  
  E_cov_parPGS_y <- r_hat^2*same_anc_terms+r_hat^2*rho_hat*cross_anc_terms
  rG_hat <- E_cov_parPGS_y/
    sqrt(r_hat^2*(E_var_parPGS_term_1+E_var_parPGS_term_2+E_var_parPGS_term_3))
  
  # Step 4: Infer lambda
  message("Step 4: Inferring mixture fraction lambda")
  lambda_hat <- (sqrt(emp_cor2_y_parPGS)-rG_hat)/(rL_hat-rG_hat)
  
  # Return all inferred parameters
  return(data.frame(R_eurPGS_HAT=hat_r2_eur,
                    R_parPGS_HAT=hat_r2_par,
                    R2_PGS_Theory=R2_PGS[1],
                    R2_PGS_Emp=R2_PGS[2],
                    R_HAT=r_hat,
                    RHO_HAT=rho_hat,
                    RG_HAT=rG_hat,
                    RL_HAT=rL_hat,
                    EMP_COR_Y_EURPGS=sqrt(emp_cor2_y_eurPGS),
                    EMP_COR_Y_PARPGS=sqrt(emp_cor2_y_parPGS),
                    LAMBDA_HAT=lambda_hat))
}

#' Generate admixed haplotype lists for all quantiles of PMBB
#' 
#' This function calls getAdmixedHapList to generate admixed haplotype quantities
#' for each quantile of PMBB. This function is run once so that the large matrix
#' operations involved will not be repeatedly performed per run. 
getPMBBAdmixedHapLists <- function(seed_haplotypes,
                                   seed_afr_ancestries) {
  # Create list of lists to return
  to_return <- vector("list",length=4)
  
  # For each quantile
  for (q in 1:4) {
    # Generate haplotype and African ancestry matrices
    haps <- cbind(seed_haplotypes[[paste0("Q",q,"H1")]], # n x 2p
                  seed_haplotypes[[paste0("Q",q,"H2")]])
    afr_ancs <- cbind(seed_afr_ancestries[[paste0("Q",q,"H1")]], # n x 2p
                      seed_afr_ancestries[[paste0("Q",q,"H2")]])
    
    # Run getAdmixedHapList on that quantile and save output
    to_return[[q]] <- getAdmixedHapList(haps,
                                        afr_ancs)
  }
  
  # Return
  return(to_return)
}

#' Generate PGSs and phenotype from scratch 
#' 
#' If `theory` is set to TRUE, this function works out the analytical expectations
#' of various PGS metrics, including squared correlations between eurPGS and parPGS
#' with the phenotype.
#'
#' Note that if INFERENCE=TRUE, then the following parameters should be specified:
#' - S (to allow loading of EurAm genotype files in a seed-specific manner)
#' - use_R2_PGS (to specify if, during the inference step, R2_PGS should be used in place of r_hat)
#'
# Example
#
# set.seed(2024)
# spop_freq_vec <- runif(n=5e3,min=0.3,max=0.5)
# tpop_freq_vec <- spop_freq_vec
doRRuns <- function(S,
                    seed_haplotypes,
                    seed_afr_ancestries,
                    n_markers,
                    r2,
                    cross_pop_rho,
                    model,
                    lambda=NULL,
                    mix_type=NULL,
                    n_reps=100,
                    use_R2_PGS=FALSE,
                    theory=FALSE,
                    inference=FALSE,
                    print_msg=FALSE) {
  # Check for errors
  if (!is.null(model)) {
    assertthat::assert_that(model %in% c("global","local"),
                            msg="Model is neither global nor local.")
  } 
  if (!is.null(lambda)) {
    assertthat::assert_that(lambda < 1 & lambda > 0,
                            msg="Lambda is not strictly between 0 and 1")
    assertthat::assert_that(is.null(model),
                            msg="Model should be set to NULL if using mixture model")
    assertthat::assert_that(!is.null(mix_type),
                            msg="Mixing type is not specified although lambda is.")
  }
  
  # Compute admixed_haps_list for each PMBB quantile (done just once!)
  admixed_haps_lists <- getPMBBAdmixedHapLists(seed_haplotypes,
                                               seed_afr_ancestries)
  
  # Create big dataframes to return
  big_simulation_df <- data.frame(var_y=numeric(),
                                  var_admPGS_orig=numeric(),
                                  var_eurPGS_orig=numeric(),
                                  var_parPGS_orig=numeric(),
                                  var_admPGS_demean=numeric(),
                                  var_eurPGS_demean=numeric(),
                                  var_parPGS_demean=numeric(),
                                  cor2_y_admPGS_orig=numeric(),
                                  cor2_y_eurPGS_orig=numeric(),
                                  cor2_y_parPGS_orig=numeric(),
                                  cor2_y_admPGS_demean=numeric(),
                                  cor2_y_eurPGS_demean=numeric(),
                                  cor2_y_parPGS_demean=numeric(),
                                  Quantile=character(),
                                  Rep=numeric())
  if (theory) {
    big_theory_df <- data.frame(cor2_y_eurPGS=numeric(),
                                cor2_y_parPGS=numeric(),
                                var_eurPGS=numeric(),
                                var_parPGS=numeric(),
                                MODEL=character(),
                                LAMBDA=numeric(),
                                A_BAR=numeric(),
                                Quantile=character(),
                                Rep=numeric())
  } else {
    big_theory_df <- NULL
  }
  
  if (inference) {
    big_infer_df <- data.frame(R_eurPGS_HAT=numeric(),
                               R_parPGS_HAT=numeric(),
                               R2_PGS_Theory=numeric(),
                               R2_PGS_Emp=numeric(),
                               R_HAT=numeric(),
                               RHO_HAT=numeric(),
                               RG_HAT=numeric(),
                               RL_HAT=numeric(),
                               LAMBDA_HAT=numeric(),
                               Quantile=character(),
                               Rep=numeric())
  } else {
    big_infer_df <- NULL
  }
  
  # If using R^2_PGS
  if (use_R2_PGS) {
    if (print_msg) {
      message(date(),": Computing prepped_data on PMBB European Americans")
    }
    
    # Compute prepped_data on EurAm (for computing R^2_PGS later)
    obj.bigSNP <- snp_attach(paste0(backup_dir,"PMBB_eur_QC_seed_",S,".rds"))
    seed_EurAm_G <- obj.bigSNP$genotypes
    PREPPED_DATA <- prepEurAmData(seed_EurAm_G)
  }
  
  
  if (print_msg) {
    start_time <- Sys.time()
    message(date(), ": Simulating ", n_reps, " phenotypes from PMBB genotypes ++++++++++++++++++++++++")
  }

  # For each rep
  for (rep in 1:n_reps) {
    if (print_msg) {
      message(date(),": Working on Simulation Rep #", rep , " ----------------------------------------")
    }
    
    # Create lists 
    simulation_list <- vector("list",4)
    theory_list <- vector("list",4)
    infer_list <- vector("list",4)
    
    # Simulate for each quantile
    for (q in 1:4) {
      # Simulate unscaled population effect sizes
      if (print_msg) {
        message(date(),": Simulating unscaled effect sizes for source and target populations in Q", q)
        message("Parameters: r2 = ", r2, ", rho = ", cross_pop_rho)
      }
      beta_vecs <- simPopBetas(n_markers=n_markers,
                               r2=r2,
                               cross_pop_rho=cross_pop_rho)
      
      # If using R^2_PGS 
      if (use_R2_PGS) {
        r2_pgs_vec <- getEurAmPGS(beta_vec=beta_vecs[["SOURCE_POP_BETA_VEC"]],
                                  r2=r2,
                                  prepped_data=PREPPED_DATA)
      }
      
      # Extract admixed hap list
      ADMIXED_HAPS_LIST <- admixed_haps_lists[[q]]
      
      # Simulate scaled effect sizes and phenotypes using specific rho and r2 and lambda
      SPOP_FREQ_VEC <- ADMIXED_HAPS_LIST[["SOURCE_POP_AF"]]
      TPOP_FREQ_VEC <- ADMIXED_HAPS_LIST[["TARGET_POP_AF"]]
      N_MARKERS <- ADMIXED_HAPS_LIST[["N_MARKERS"]]
      N_IND <- ADMIXED_HAPS_LIST[["N_IND"]]
      
      TPOP_FRAC <- ADMIXED_HAPS_LIST[["TPOP_FRAC"]] # scalar
      
      # Check that n_markers used to generate effect size vector
      # and N_MARKERS computed from getAdmixedHapList match
      assertthat::assert_that(N_MARKERS==n_markers,
                              msg="Mismatch between N_MARKERS and n_markers.")
      
      # Create length 2p beta vectors for simPGS2
      spop_beta_vec_h1 <- beta_vecs[["SOURCE_POP_BETA_VEC"]]/
        sqrt(sum(SPOP_FREQ_VEC[1:N_MARKERS]*(1-SPOP_FREQ_VEC[1:N_MARKERS])))
      spop_beta_vec_h2 <- beta_vecs[["SOURCE_POP_BETA_VEC"]]/
        sqrt(sum(SPOP_FREQ_VEC[(N_MARKERS+1):(2*N_MARKERS)]*(1-SPOP_FREQ_VEC[(N_MARKERS+1):(2*N_MARKERS)])))
      spop_beta_vec <- c(spop_beta_vec_h1,spop_beta_vec_h2)
      
      tpop_beta_vec_h1 <- beta_vecs[["TARGET_POP_BETA_VEC"]]/
        sqrt(sum(TPOP_FREQ_VEC[1:N_MARKERS]*(1-TPOP_FREQ_VEC[1:N_MARKERS])))
      tpop_beta_vec_h2 <- beta_vecs[["TARGET_POP_BETA_VEC"]]/
        sqrt(sum(TPOP_FREQ_VEC[(N_MARKERS+1):(2*N_MARKERS)]*(1-TPOP_FREQ_VEC[(N_MARKERS+1):(2*N_MARKERS)])))
      tpop_beta_vec <- c(tpop_beta_vec_h1,tpop_beta_vec_h2)
      
      # Simulate admixed PGS
      if (print_msg) {
        message(date(),": Simulating phenotype and PGSs for subsample (Q", q, ")")
        if (!is.null(model)) {
          message("Parameters: model = ", model)
        } else {
          message("Parameters: lambda = ", lambda, ", mixing type = ", mix_type)
        }
      }
      sim_pgs <- simPGS2(r2=r2,
                         cross_pop_rho=cross_pop_rho,
                         admixed_haps_list=ADMIXED_HAPS_LIST,
                         source_pop_beta_vec=spop_beta_vec,
                         target_pop_beta_vec=tpop_beta_vec,
                         model=model,
                         lambda=lambda,
                         mix_type=mix_type)
      
      # Compute quantities to be estimated
      if (print_msg) {
        message(date(),": Estimating variance and sq. correlation quantities")
      }
      est_quantities <- data.frame(var_y=var(sim_pgs$Y),
                                   var_admPGS_orig=var(sim_pgs$admPGS_orig),
                                   var_eurPGS_orig=var(sim_pgs$eurPGS_orig),
                                   var_parPGS_orig=var(sim_pgs$parPGS_orig),
                                   var_admPGS_demean=var(sim_pgs$admPGS_demean),
                                   var_eurPGS_demean=var(sim_pgs$eurPGS_demean),
                                   var_parPGS_demean=var(sim_pgs$parPGS_demean),
                                   cor2_y_admPGS_orig=cor(sim_pgs$Y,sim_pgs$admPGS_orig)^2,
                                   cor2_y_eurPGS_orig=cor(sim_pgs$Y,sim_pgs$eurPGS_orig)^2,
                                   cor2_y_parPGS_orig=cor(sim_pgs$Y,sim_pgs$parPGS_orig)^2,
                                   cor2_y_admPGS_demean=cor(sim_pgs$Y,sim_pgs$admPGS_demean)^2,
                                   cor2_y_eurPGS_demean=cor(sim_pgs$Y,sim_pgs$eurPGS_demean)^2,
                                   cor2_y_parPGS_demean=cor(sim_pgs$Y,sim_pgs$parPGS_demean)^2)
      
      # Compute analytical expectations
      # In theory, this can be done just once since it's not dependent on random beta
      # Am putting it here just because I want to check that the results are consistent
      # across the reps. 
      if (theory) {
        if (print_msg) {
          message(date(),": Computing analytical expectations under generative model")
        }
        expected_quantities <- computeExpected(r2=r2,
                                               cross_pop_rho=cross_pop_rho,
                                               admixed_haps_list=ADMIXED_HAPS_LIST,
                                               model=model,
                                               lambda=lambda)
        expected_quantities$MODEL <- model
        expected_quantities$LAMBDA <- lambda
        expected_quantities$A_BAR <- TPOP_FRAC
      } else {
        expected_quantities <- NULL
      }
      
      # Infer parameters
      if (inference) {
        if (print_msg) {
          message(date(), ": Performing inference on simulated data")
        }
        
        if (use_R2_PGS) {
          inferred_params <- inferParams(emp_cor2_y_eurPGS=est_quantities$cor2_y_eurPGS_demean,
                                         emp_cor2_y_parPGS=est_quantities$cor2_y_parPGS_demean,
                                         emp_var_eurPGS=est_quantities$var_eurPGS_demean,
                                         emp_var_parPGS=est_quantities$var_parPGS_demean,
                                         admixed_haps_list=ADMIXED_HAPS_LIST,
                                         snp_herit=NA,
                                         ignore_var=TRUE,
                                         R2_PGS=r2_pgs_vec,
                                         use_R2_PGS=TRUE)
        } else {
          inferred_params <- inferParams(emp_cor2_y_eurPGS=est_quantities$cor2_y_eurPGS_demean,
                                         emp_cor2_y_parPGS=est_quantities$cor2_y_parPGS_demean,
                                         emp_var_eurPGS=est_quantities$var_eurPGS_demean,
                                         emp_var_parPGS=est_quantities$var_parPGS_demean,
                                         admixed_haps_list=ADMIXED_HAPS_LIST,
                                         snp_herit=NA,
                                         ignore_var=FALSE,
                                         R2_PGS=c(NA,NA),
                                         use_R2_PGS=FALSE)
        }
        
      } else {
        inferred_params <- NULL
      }
      
      # Append quantile and rep information
      message(date(), ": Appending Quantile and Rep info for Q",q)
      est_quantities$Quantile <- paste0("Q",q)
      est_quantities$Rep <- rep
      
      if (theory) {
        expected_quantities$Quantile <- paste0("Q",q)
        expected_quantities$Rep <- rep
      }
      if (inference) {
        inferred_params$Quantile <- paste0("Q",q)
        inferred_params$Rep <- rep
      }
      
      simulation_list[[q]] <- est_quantities
      theory_list[[q]] <- expected_quantities
      infer_list[[q]] <- inferred_params
    }
    
    # Combine everything
    simulation_unlist <- do.call(rbind,simulation_list)
    theory_unlist <- do.call(rbind,theory_list)
    infer_unlist <- do.call(rbind,infer_list)
    
    if (print_msg) {
      message(date(),": Finished Simulation Rep #", rep , " ----------------------------------------")
    }
    
    # Append to large dataframe that will be returned
    big_simulation_df <- rbind(big_simulation_df,simulation_unlist)
    if (theory) {
      big_theory_df <- rbind(big_theory_df,theory_unlist)
    }
    if (inference) {
      big_infer_df <- rbind(big_infer_df,infer_unlist)
    }
  }
  if (print_msg) {
    end_time <- Sys.time()
    message("Total time elapsed = ", 
            round(difftime(end_time,start_time,units="mins"),digits=3),
            " mins")
  }
  
  # Return
  return(list(SIMULATION=big_simulation_df,
              THEORY=big_theory_df,
              INFERENCE=big_infer_df))
}
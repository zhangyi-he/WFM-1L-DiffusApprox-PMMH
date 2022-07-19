#' @title Estimating temporally variable selection intensity from ancient DNA data
#' @author Zhangyi He, Xiaoyang Dai, Wenyang Lyu, Mark Beaumont, Feng Yu

#' version 1.1 (PMMH)
#' Phenotypes controlled by a single gene
#' Constant natural selection and non-constant demographic histories

#' Input: called genotypes
#' Output: posteriors for the selection coefficient

#' R functions

#install.packages("gtools")
library("gtools")

#install.packages("purrr")
library("purrr")

#install.packages("MASS")
library("MASS")

#install.packages("coda")
library("coda")

#install.packages("inline")
library("inline")
#install.packages("Rcpp")
library("Rcpp")
#install.packages("RcppArmadillo")
library("RcppArmadillo")

#install.packages("compiler")
library("compiler")
#enableJIT(1)

# call C++ functions
sourceCpp("./Code v1.1/CFUN.cpp")

################################################################################

#' Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_frq the initial mutant allele frequency of the population
#' @param int_gen the generation that the simulated mutant allele frequency trajectory started
#' @param lst_gen the generation that the simulated mutant allele frequency trajectory ended

#' Standard version
simulateWFM <- function(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen) {
  fts_mat <- calculateFitnessMat_arma(sel_cof, dom_par)

  WFM <- simulateWFM_arma(fts_mat, pop_siz, int_frq, int_gen, lst_gen)
  mut_frq_pth <- as.vector(WFM$mut_frq_pth)
  gen_frq_pth <- as.matrix(WFM$gen_frq_pth)

  return(list(mut_frq_pth = mut_frq_pth,
              gen_frq_pth = gen_frq_pth))
}
#' Compiled version
cmpsimulateWFM <- cmpfun(simulateWFM)

########################################

#' Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param int_frq the initial mutant allele frequency of the population
#' @param int_gen the generation that the simulated mutant allele frequency trajectory started
#' @param lst_gen the generation that the simulated mutant allele frequency trajectory ended
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

#' Standard version
simulateWFD <- function(sel_cof, dom_par, pop_siz, ref_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE) {
  mut_frq_pth <- simulateWFD_arma(sel_cof, dom_par, pop_siz, ref_siz, int_frq, int_gen, lst_gen, ptn_num)
  mut_frq_pth <- as.vector(mut_frq_pth)
  
  if (dat_aug == FALSE) {
    return(mut_frq_pth[(0:(lst_gen - int_gen)) * ptn_num + 1])
  } else {
    return(mut_frq_pth)
  }  
}
#' Compiled version
cmpsimulateWFD <- cmpfun(simulateWFD)

########################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_con the initial mutant allele frequency of the population
#' @param smp_lab the identifier of the sample assigned
#' @param smp_gen the generation of the sample drawn
#' @param smp_qua the quality of the sample tested
#' @param thr_val the threshold for genotype calling
#' @param ref_siz the reference size of the horse population
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Standard version
simulateHMM <- function(model, sel_cof, dom_par, pop_siz, int_con, smp_lab, smp_gen, smp_qua, thr_val, ...) {
  int_gen <- min(smp_gen)
  lst_gen <- max(smp_gen)
  
  # generate the population genotype frequency trajectories
  if (model == "WFM") {
    WFM <- cmpsimulateWFM(sel_cof, dom_par, pop_siz, int_con, int_gen, lst_gen)
    mut_frq <- as.vector(WFM$mut_frq_pth)
    gen_frq <- as.matrix(WFM$gen_frq_pth)
  }
  if (model == "WFD") {
    mut_frq <- cmpsimulateWFD(sel_cof, dom_par, pop_siz, ref_siz, int_con, int_gen, lst_gen, ptn_num, dat_aug = FALSE)
    mut_frq <- as.vector(mut_frq)
    gen_frq <- matrix(NA, nrow = 3, ncol = length(mut_frq))
    fts_mat <- calculateFitnessMat_arma(sel_cof, dom_par)
    for (k in 1:length(mut_frq)) {
      pop_ale_frq <- c(1 - mut_frq[k], mut_frq[k])
      pop_gen_frq <- fts_mat * (pop_ale_frq %*% t(pop_ale_frq)) / sum(fts_mat * (pop_ale_frq %*% t(pop_ale_frq)))
      pop_gen_frq[lower.tri(pop_gen_frq, diag = FALSE)] <- NA
      gen_frq[, k] <- discard(as.vector(2 * pop_gen_frq - diag(diag(pop_gen_frq), nrow = 2, ncol = 2)), is.na)
    }
    gen_frq <- as.matrix(gen_frq)
  }
  
  # generate the sample individual genotypes at all sampling time points
  tru_smp <- NULL
  for (k in 1:length(smp_gen)) {
    tru_smp <- cbind(tru_smp, rbind(rep(smp_gen[k], times = 1), rmultinom(1, size = 1, prob = gen_frq[, smp_gen[k] - int_gen + 1])))
  }
  tru_smp <- as.data.frame(t(tru_smp))
  rownames(tru_smp) <- smp_lab
  colnames(tru_smp) <- c("generation", "A0A0", "A0A1", "A1A1")
  tru_smp <- tru_smp[order(tru_smp$generation), ]

  raw_smp <- tru_smp
  sel_idx <- which(tru_smp$'A0A0' == 1)
  raw_smp[sel_idx, 2:4] <- rdirichlet(length(sel_idx), alpha = c(smp_qua[1], (1 - smp_qua[1]) / 2, (1 - smp_qua[1]) / 2) * smp_qua[2])
  sel_idx <- which(tru_smp$'A0A1' == 1)
  raw_smp[sel_idx, 2:4] <- rdirichlet(length(sel_idx), alpha = c((1 - smp_qua[1]) / 2, smp_qua[1], (1 - smp_qua[1]) / 2) * smp_qua[2])
  sel_idx <- which(tru_smp$'A1A1' == 1)
  raw_smp[sel_idx, 2:4] <- rdirichlet(length(sel_idx), alpha = c((1 - smp_qua[1]) / 2, (1 - smp_qua[1]) / 2, smp_qua[1]) * smp_qua[2])

  cal_smp <- raw_smp
  mis_idx <- 1:nrow(raw_smp)
  sel_idx <- which(raw_smp$'A0A0' / raw_smp$'A0A1' >= thr_val & raw_smp$'A0A0' / raw_smp$'A1A1' >= thr_val)
  if (length(sel_idx) > 1) {
    cal_smp[sel_idx, 2:4] <- matrix(rep(c(1, 0, 0), times = length(sel_idx)), nrow = length(sel_idx), ncol = 3, byrow = TRUE)
    mis_idx <- mis_idx[-which(mis_idx %in% sel_idx)]
  }
  sel_idx <- which(raw_smp$'A0A1' / raw_smp$'A0A0' >= thr_val & raw_smp$'A0A1' / raw_smp$'A1A1' >= thr_val)
  if (length(sel_idx) > 1) {
    cal_smp[sel_idx, 2:4] <- matrix(rep(c(0, 1, 0), times = length(sel_idx)), nrow = length(sel_idx), ncol = 3, byrow = TRUE)
    mis_idx <- mis_idx[-which(mis_idx %in% sel_idx)]
  }
  sel_idx <- which(raw_smp$'A1A1' / raw_smp$'A0A0' >= thr_val & raw_smp$'A1A1' / raw_smp$'A0A1' >= thr_val)
  if (length(sel_idx) > 1) {
    cal_smp[sel_idx, 2:4] <- matrix(rep(c(0, 0, 1), times = length(sel_idx)), nrow = length(sel_idx), ncol = 3, byrow = TRUE)
    mis_idx <- mis_idx[-which(mis_idx %in% sel_idx)]
  }
  if (length(mis_idx) > 1) {
    cal_smp[mis_idx, 2:4] <- matrix(rep(c(NA, NA, NA), times = length(mis_idx)), nrow = length(mis_idx), ncol = 3, byrow = TRUE)
  }

  fil_smp <- cal_smp[complete.cases(cal_smp), ]

  cal_rat <- nrow(fil_smp) / nrow(cal_smp)
  # cal_rat
  err_rat <- sum(rowSums(cal_smp[, -1] == tru_smp[, -1]) == 1, na.rm = TRUE) / nrow(fil_smp)
  # err_rat

  return(list(tru_smp = tru_smp,
              raw_smp = raw_smp,
              cal_smp = cal_smp,
              fil_smp = fil_smp,
              cal_rat = cal_rat,
              err_rat = err_rat,
              gen_frq = gen_frq,
              mut_frq = mut_frq))
}
#' Compiled version
cmpsimulateHMM <- cmpfun(simulateHMM)

########################################

#' Run the bootstrap particle filter (BPF) with the single-locus Wright-Fisher diffusion with selection
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

#' Standard version
runBPF <- function(sel_cof, dom_par, pop_siz, ref_siz, raw_smp, ptn_num, pcl_num) {
  # preprocess the raw sample
  raw_smp <- as.data.frame(raw_smp)
  colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  raw_smp <- raw_smp[complete.cases(raw_smp), ]
  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  # run the BPF
  BPF <- runBPF_arma(sel_cof, dom_par, pop_siz, ref_siz, raw_smp, ptn_num, pcl_num)
  
  return(list(lik = BPF$lik, 
              wght = BPF$wght, 
              mut_frq_pre_resmp = BPF$mut_frq_pre_resmp, 
              mut_frq_pst_resmp = BPF$mut_frq_pst_resmp,
              gen_frq_pre_resmp = BPF$gen_frq_pre_resmp,
              gen_frq_pst_resmp = BPF$gen_frq_pst_resmp))
}
#' Compiled version
cmprunBPF <- cmpfun(runBPF)

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse populations
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

#' Standard version
calculateOptimalParticleNum <- function(sel_cof, dom_par, pop_siz, ref_siz, raw_smp, ptn_num, pcl_num, gap_num) {
  # preprocess the raw sample
  raw_smp <- as.data.frame(raw_smp)
  colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  raw_smp <- raw_smp[complete.cases(raw_smp), ]
  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  # calculate the optimal particle number
  OptNum <- calculateOptimalParticleNum_arma(sel_cof, dom_par, pop_siz, ref_siz, raw_smp, ptn_num, pcl_num, gap_num)
  
  return(list(opt_pcl_num = as.vector(OptNum$opt_pcl_num), 
              log_lik_sdv = as.vector(OptNum$log_lik_sdv)))
}
#' Compiled version
cmpcalculateOptimalParticleNum <- cmpfun(calculateOptimalParticleNum)

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH

#' Standard version
runPMMH <- function(sel_cof, dom_par, pop_siz, ref_siz, raw_smp, ptn_num, pcl_num, itn_num) {
  # preprocess the raw sample
  raw_smp <- as.data.frame(raw_smp)
  colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  raw_smp <- raw_smp[complete.cases(raw_smp), ]
  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))
  
  # run the PMMH
  PMMH <- runPMMH_arma(sel_cof, dom_par, pop_siz, ref_siz, raw_smp, ptn_num, pcl_num, itn_num)
  sel_cof_chn <- as.vector(PMMH)
  
  return(sel_cof_chn)
}
#' Compiled version
cmprunPMMH <- cmpfun(runPMMH)

########################################

#' Run the Bayesian Procedure for the inference of selection
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning

#' Standard version
runBayesianProcedure <- function(sel_cof, dom_par, pop_siz, ref_siz, raw_smp, ptn_num, pcl_num, itn_num, brn_num, thn_num) {
  # preprocess the raw sample
  raw_smp <- as.data.frame(raw_smp)
  colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  raw_smp <- raw_smp[complete.cases(raw_smp), ]
  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  # run the PMMH
  PMMH <- runPMMH_arma(sel_cof, dom_par, pop_siz, ref_siz, raw_smp, ptn_num, pcl_num, itn_num)
  sel_cof_chn <- as.vector(PMMH)

  # burn-in and thinning
  sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
  sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]
  
  # MMSE estimate for selection coefficient
  sel_cof_est <- mean(sel_cof_chn)
  
  # 95% HPD interval for selection coefficient
  sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
  
  return(list(sel_cof_est = sel_cof_est, 
              sel_cof_hpd = sel_cof_hpd, 
              sel_cof_chn = sel_cof_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################

#' @title Estimating temporally variable selection intensity from ancient DNA data
#' @author Zhangyi He, Xiaoyang Dai, Wenyang Lyu, Mark Beaumont, Feng Yu

#' version 1.3 (PMMH)
#' Phenotypes controlled by a single gene
#' Non-constant natural selection and non-constant demographic histories

#' Input: called genotypes
#' Output: posteriors for the selection coefficient and the genotype frequency trajectories of the population

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

#install.packages("emdbook")
library("emdbook")

# call R functions
source("./RFUN.R")

################################################################################

#' Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_frq the initial mutant allele frequency of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation that the simulated mutant allele frequency trajectory started
#' @param lst_gen the generation that the simulated mutant allele frequency trajectory ended

sel_cof <- c(0e+00, 1e-02)
dom_par <- 0e+00
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_frq <- 2e-01
evt_gen <- 240
int_gen <- 0
lst_gen <- 500

frq_pth <- cmpsimulateWFM(sel_cof, dom_par, pop_siz, int_frq, evt_gen, int_gen, lst_gen)$mut_frq_pth

k <- int_gen:lst_gen
plot(k, frq_pth, type = 'l', lwd = 1.5,
     xlab = "Generation", ylab = "Allele frequency",
     main = "WFM: the mutant allele frequency trajectory")

########################################

#' Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param int_frq the initial mutant allele frequency of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation that the simulated mutant allele frequency trajectory started
#' @param lst_gen the generation that the simulated mutant allele frequency trajectory ended
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

sel_cof <- c(0e+00, 1e-02)
dom_par <- 0e+00
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
ref_siz <- 1e+04
int_frq <- 2e-01
evt_gen <- 240
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

frq_pth <- cmpsimulateWFD(sel_cof, dom_par, pop_siz, ref_siz, int_frq, evt_gen, int_gen, lst_gen, ptn_num, dat_aug = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / ref_siz
plot(t, frq_pth, type = 'l', lwd = 1.5,
     xlab = "Generation", ylab = "Allele frequency",
     main = "WFD: the mutant allele frequency trajectory")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof <- c(0e+00, 1e-02)
dom_par <- 0e+00
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
ref_siz <- 1e+04
int_frq <- 2e-01
evt_gen <- 240
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00
sim_num <- 1e+06

smp_WFM <- numeric(sim_num)
smp_WFD <- numeric(sim_num)
for (i in 1:sim_num) {
  print(i)
  smp_WFM[i] <- cmpsimulateWFM(sel_cof, dom_par, pop_siz, int_frq, evt_gen, int_gen, lst_gen)$mut_frq[(lst_gen - int_gen) + 1]
  smp_WFD[i] <- cmpsimulateWFD(sel_cof, dom_par, pop_siz, ref_siz, int_frq, evt_gen, int_gen, lst_gen, ptn_num, dat_aug = FALSE)[(lst_gen - int_gen) + 1]
}

hist(smp_WFM, breaks = seq(min(smp_WFM, smp_WFD), max(smp_WFM, smp_WFD), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_WFM, smp_WFD), max(smp_WFM, smp_WFD)),
     xlab = "Allele frequency", main = paste("Histograms of the mutant allele frequency at generation", lst_gen))
hist(smp_WFD, breaks = seq(min(smp_WFM, smp_WFD), max(smp_WFM, smp_WFD), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

################################################################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_con the initial mutant allele frequency of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_lab the identifier of the sample assigned
#' @param smp_gen the generation of the sample drawn
#' @param smp_qua the quality of the sample tested
#' @param thr_val the threshold for genotype calling
#' @param ref_siz the reference size of the horse population
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Simulate the dataset under the Wright-Fisher model
model <- "WFM"
sel_cof <- c(0e+00, 1e-02)
dom_par <- 0e+00
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_con <- 2e-01
evt_gen <- 240
smp_lab <- 1:220
smp_gen <- rep((0:10) * 50, each = 20)
smp_qua <- c(0.95, 1)
thr_val <- 10

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, smp_qua, thr_val)
tru_smp <- aggregate(. ~ generation, data = sim_HMM_WFM$tru_smp, sum)
smp_gen <- tru_smp[, 1]
smp_siz <- rowSums(tru_smp[, -1])
smp_cnt <- t(as.matrix(tru_smp[, 2:4]))
smp_frq <- smp_cnt %*% diag(1 / smp_siz)
pop_frq <- sim_HMM_WFM$gen_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[1, ], pop_frq[1, ]), max(smp_frq[1, ], pop_frq[1, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the A0A0 genotype")
points(smp_gen, smp_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[2, ], pop_frq[2, ]), max(smp_frq[2, ], pop_frq[2, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the A0A1 genotype")
points(smp_gen, smp_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[3, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[3, ], pop_frq[3, ]), max(smp_frq[3, ], pop_frq[3, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the A1A1 genotype")
points(smp_gen, smp_frq[3, ], col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion
model <- "WFD"
sel_cof <- c(0e+00, 1e-02)
dom_par <- 0e+00
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_con <- 2e-01
evt_gen <- 240
smp_lab <- 1:220
smp_gen <- rep((0:10) * 50, each = 20)
smp_qua <- c(0.95, 1)
thr_val <- 10
ref_siz <- 1e+04
ptn_num <- 5e+00

sim_HMM_WFD <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, smp_qua, thr_val, ref_siz, ptn_num)
tru_smp <- aggregate(. ~ generation, data = sim_HMM_WFD$tru_smp, sum)
smp_gen <- tru_smp[, 1]
smp_siz <- rowSums(tru_smp[, -1])
smp_cnt <- t(as.matrix(tru_smp[, 2:4]))
smp_frq <- smp_cnt %*% diag(1 / smp_siz)
pop_frq <- sim_HMM_WFD$gen_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[1, ], pop_frq[1, ]), max(smp_frq[1, ], pop_frq[1, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFD-HMM: the A0A0 genotype")
points(smp_gen, smp_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[2, ], pop_frq[2, ]), max(smp_frq[2, ], pop_frq[2, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFD-HMM: the A0A1 genotype")
points(smp_gen, smp_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[3, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[3, ], pop_frq[3, ]), max(smp_frq[3, ], pop_frq[3, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFD-HMM: the A1A1 genotype")
points(smp_gen, smp_frq[3, ], col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset under the Wright-Fisher model
test_seed <- 2
set.seed(test_seed)

model <- "WFM"
sel_cof <- c(0e+00, 1e-02)
dom_par <- 0e+00
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_con <- 2e-01
evt_gen <- 240
smp_lab <- 1:220
smp_gen <- rep((0:10) * 50, each = 20)
smp_qua <- c(0.95, 1)
thr_val <- 10

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, smp_qua, thr_val)

save(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, smp_qua, thr_val, sim_HMM_WFM,
     file = "./TEST_SimData.rda")

load("./TEST_SimData.rda")

tru_smp <- aggregate(. ~ generation, data = sim_HMM_WFM$tru_smp, sum)
smp_gen <- tru_smp[, 1]
smp_siz <- rowSums(tru_smp[, -1])
smp_cnt <- t(as.matrix(tru_smp[, 2:4]))
smp_frq <- smp_cnt %*% diag(1 / smp_siz)
pop_frq <- sim_HMM_WFM$gen_frq

pdf(file = "./TEST_SimData.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
k <- min(smp_gen):max(smp_gen)
plot(k, pop_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[1, ], pop_frq[1, ]), max(smp_frq[1, ], pop_frq[1, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the A0A0 genotype")
points(smp_gen, smp_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[2, ], pop_frq[2, ]), max(smp_frq[2, ], pop_frq[2, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the A0A1 genotype")
points(smp_gen, smp_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_frq[3, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[3, ], pop_frq[3, ]), max(smp_frq[3, ], pop_frq[3, ])),
     xlab = "Generation", ylab = "Genotype frequency",
     main = "WFM-HMM: the A1A1 genotype")
points(smp_gen, smp_frq[3, ], col = 'red', pch = 17, cex = 1)
dev.off()

########################################

#' Run the bootstrap particle filter (BPF) with the single-locus Wright-Fisher diffusion with selection
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

load("./TEST_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
pop_siz
ref_siz <- 1e+04
evt_gen <- 240
raw_smp <- sim_HMM_WFM$tru_smp
ptn_num <- 5e+00
pcl_num <- 1e+05

system.time(BPF <- cmprunBPF(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, BPF,
     file = "./TEST_BPF.rda")

load("./TEST_BPF.rda")

tru_smp <- aggregate(. ~ generation, data = sim_HMM_WFM$tru_smp, sum)
smp_gen <- tru_smp[, 1]
smp_siz <- rowSums(tru_smp[, -1])
smp_cnt <- t(as.matrix(tru_smp[, 2:4]))
smp_frq <- smp_cnt %*% diag(1 / smp_siz)

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}

pdf(file = "./TEST_BPF_Likelihood.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:pcl_num, log(lik), type = 'l',
     xlab = "Number of particles", ylab = "Log likelihood",
     main = "Log likelihood through the bootstrap particle filter")
dev.off()

pop_frq_pre_resmp <- BPF$gen_frq_pre_resmp
pop_frq_pst_resmp <- BPF$gen_frq_pst_resmp

pdf(file = "./TEST_BPF_Particle.pdf", width = 24, height = 66)
par(mfrow = c(11, 3), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
for (k in 1:length(smp_gen)) {
  hist_pst_resmp <- hist(pop_frq_pst_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k], smp_frq[1, k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k], smp_frq[1, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype A0A0 at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[1, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_frq_pst_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k], smp_frq[2, k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k], smp_frq[2, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype A0A1 at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[2, k], col = 'red', lty = 2, lwd = 2)

  hist_pst_resmp <- hist(pop_frq_pst_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k], smp_frq[3, k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k], smp_frq[3, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Genotype frequency",
       main = paste("Genotype A1A1 at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[3, k], col = 'red', lty = 2, lwd = 2)
}
dev.off()

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

load("./TEST_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
pop_siz
ref_siz <- 1e+04
evt_gen <- 240
raw_smp <- sim_HMM_WFM$tru_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
gap_num <- 1e+02

system.time(OptNum <- calculateOptimalParticleNum(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, gap_num))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, gap_num, OptNum,
     file = "./TEST_OptNum.rda")

load("./TEST_OptNum.rda")

opt_pcl_num <- OptNum$opt_pcl_num
log_lik_sdv <- OptNum$log_lik_sdv

pdf(file = "./TEST_OptNum.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(opt_pcl_num, log_lik_sdv, type = 'b', lwd = 2,
     xlab = "Particle number", ylab = "Log-likelihood standard deviation",
     main = "Optimal particle number in the PMMH")
abline(h = 1.7, col = 'red', lty = 2, lwd = 2)
abline(h = 1.0, col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH

load("./TEST_SimData.rda")

set.seed(test_seed)

sel_cof <- c(0e+00, 0e+00)
dom_par
pop_siz
ref_siz <- 1e+04
evt_gen <- 240
raw_smp <- sim_HMM_WFM$tru_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num))

load("./TEST_SimData.rda")

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, PMMH,
     file = "./TEST_PMMH.rda")

load("./TEST_PMMH.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./TEST_PMMH_Traceplot.pdf", width = 16, height = 6)
par(mfrow = c(1, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before event")
abline(h = sel_cof[1], col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after event")
abline(h = sel_cof[2], col = 'red', lty = 2, lwd = 2)
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./TEST_PMMH_Posterior_SelCoeff.pdf", width = 16, height = 6)
par(mfrow = c(1, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before event")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after event")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./TEST_PMMH_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = sel_cof[2] - sel_cof[1], col = 'red', lty = 2, lwd = 2)
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

frq_pth_est <- colMeans(frq_pth_chn)

frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
  frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./TEST_PMMH_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = 'n', xlim = c(min(raw_smp[, 1]), max(raw_smp[, 1])), ylim = c(min(frq_pth_chn, sim_HMM_WFM$mut_frq), max(frq_pth_chn, sim_HMM_WFM$mut_frq)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")
for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(raw_smp[, 1]):max(raw_smp[, 1]), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp[, 1]):max(raw_smp[, 1]), sim_HMM_WFM$mut_frq, col = 'red', lty = 1, lwd = 2)
lines(min(raw_smp[, 1]):max(raw_smp[, 1]), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp[, 1]):max(raw_smp[, 1]), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp[, 1]):max(raw_smp[, 1]), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the adaptive particle marginal Metropolis-Hastings (AdaptPMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

load("./TEST_SimData.rda")

set.seed(test_seed)

sel_cof <- c(0e+00, 0e+00)
dom_par
pop_siz
ref_siz <- 1e+04
evt_gen <- 240
raw_smp <- sim_HMM_WFM$tru_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

load("./TEST_SimData.rda")

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./TEST_AdaptPMMH.rda")

load("./TEST_AdaptPMMH.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./TEST_AdaptPMMH_Traceplot.pdf", width = 16, height = 6)
par(mfrow = c(1, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before event")
abline(h = sel_cof[1], col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after event")
abline(h = sel_cof[2], col = 'red', lty = 2, lwd = 2)
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./TEST_AdaptPMMH_Posterior_SelCoeff.pdf", width = 16, height = 6)
par(mfrow = c(1, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before event")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after event")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./TEST_AdaptPMMH_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = sel_cof[2] - sel_cof[1], col = 'red', lty = 2, lwd = 2)
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

frq_pth_est <- colMeans(frq_pth_chn)

frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
  frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./TEST_AdaptPMMH_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = 'n', xlim = c(min(raw_smp[, 1]), max(raw_smp[, 1])), ylim = c(min(frq_pth_chn, sim_HMM_WFM$mut_frq), max(frq_pth_chn, sim_HMM_WFM$mut_frq)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")
for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(raw_smp[, 1]):max(raw_smp[, 1]), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp[, 1]):max(raw_smp[, 1]), sim_HMM_WFM$mut_frq, col = 'red', lty = 1, lwd = 2)
lines(min(raw_smp[, 1]):max(raw_smp[, 1]), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp[, 1]):max(raw_smp[, 1]), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp[, 1]):max(raw_smp[, 1]), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the Bayesian Procedure for the inference of selection
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param adp_set = TRUE/FALSE (return the result with the adaptive setting or not)
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

load("./TEST_SimData.rda")

set.seed(test_seed)

sel_cof <- c(0e+00, 0e+00)
dom_par
pop_siz
ref_siz <- 1e+04
evt_gen <- 240
raw_smp <- sim_HMM_WFM$tru_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
brn_num <- 1e+04
thn_num <- 5e+00
adp_set <- TRUE
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(BayesianProcedure <- cmprunBayesianProcedure(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, brn_num, thn_num, adp_set, stp_siz, apt_rto))

load("./TEST_SimData.rda")

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, brn_num, thn_num, adp_set, stp_siz, apt_rto, BayesianProcedure,
     file = "./TEST_BayesProc.rda")

load("./TEST_BayesProc.rda")

sel_cof_chn <- BayesianProcedure$sel_cof_chn
sel_cof_est <- BayesianProcedure$sel_cof_est
sel_cof_hpd <- BayesianProcedure$sel_cof_hpd

pdf(file = "./TEST_BayesProc_Posterior.pdf", width = 16, height = 6)
par(mfrow = c(1, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before event")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof[1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after event")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof[2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- BayesianProcedure$dif_sel_chn
dif_sel_est <- BayesianProcedure$dif_sel_est
dif_sel_hpd <- BayesianProcedure$dif_sel_hpd

pdf(file = "./TEST_BayesProc_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = sel_cof[2] - sel_cof[1], col = 'red', lty = 2, lwd = 2)
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

frq_pth_chn <- BayesianProcedure$frq_pth_chn
frq_pth_est <- BayesianProcedure$frq_pth_est
frq_pth_hpd <- BayesianProcedure$frq_pth_hpd

pdf(file = "./TEST_BayesProc_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = 'n', xlim = c(min(raw_smp[, 1]), max(raw_smp[, 1])), ylim = c(min(frq_pth_chn, sim_HMM_WFM$mut_frq), max(frq_pth_chn, sim_HMM_WFM$mut_frq)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")
for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(raw_smp[, 1]):max(raw_smp[, 1]), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(raw_smp[, 1]):max(raw_smp[, 1]), sim_HMM_WFM$mut_frq, col = 'red', lty = 1, lwd = 2)
lines(min(raw_smp[, 1]):max(raw_smp[, 1]), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(raw_smp[, 1]):max(raw_smp[, 1]), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(raw_smp[, 1]):max(raw_smp[, 1]), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################

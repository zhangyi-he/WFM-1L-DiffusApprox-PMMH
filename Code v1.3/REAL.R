#' @title Estimating temporally variable selection intensity from ancient DNA data
#' @author Zhangyi He, Xiaoyang Dai, Wenyang Lyu, Mark Beaumont, Feng Yu

#' version 1.3 (PMMH)
#' Phenotypes controlled by a single gene
#' Non-constant natural selection and non-constant demographic histories

#' Input: called genotypes
#' Output: posteriors for the selection coefficient and the genotype frequency trajectories of the population

#' Horse base coat colours (ASIP & MC1R) and white spotting patterns (KIT13 & KIT16 & TRMP1)

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
source("./Code v1.3/RFUN.R")

################################################################################

#' ASIP

#' Raw data of Wutke et al. (2016) from 12500 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- ASIP
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0e+00
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_COL_ASIP_1.rda")

load("./REAL_COL_ASIP_1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_COL_ASIP_1_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient from domestication")
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

pdf(file = "./REAL_COL_ASIP_1_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient from domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient from domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./REAL_COL_ASIP_1_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
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

pdf(file = "./REAL_COL_ASIP_1_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 9322 BC (Holocene 9700 BC)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- ASIP
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp <- raw_smp[which(raw_smp$age_mean <= 9700 + 2000), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0e+00
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_COL_ASIP_2.rda")

load("./REAL_COL_ASIP_2.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_COL_ASIP_2_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient from domestication")
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

pdf(file = "./REAL_COL_ASIP_2_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient from domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient from domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./REAL_COL_ASIP_2_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
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

pdf(file = "./REAL_COL_ASIP_2_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

############################################################

#' MC1R

#' Raw data of Wutke et al. (2016) from 4300 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- MC1R
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0e+00
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_COL_MC1R_1.rda")

load("./REAL_COL_MC1R_1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_COL_MC1R_1_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient from domestication")
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

pdf(file = "./REAL_COL_MC1R_1_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient from domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient from domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./REAL_COL_MC1R_1_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
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

pdf(file = "./REAL_COL_MC1R_1_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 9322 BC (Holocene 9700 BC)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- MC1R
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp <- raw_smp[which(raw_smp$age_mean <= 9700 + 2000), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0e+00
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_COL_MC1R_2.rda")

load("./REAL_COL_MC1R_2.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_COL_MC1R_2_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient from domestication")
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

pdf(file = "./REAL_COL_MC1R_2_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient from domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient from domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./REAL_COL_MC1R_2_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
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

pdf(file = "./REAL_COL_MC1R_2_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

############################################################

#' KIT13

#' Raw data of Wutke et al. (2016) from 3645 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT13
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1e+00
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_PTN_KIT13_1.rda")

load("./REAL_PTN_KIT13_1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_PTN_KIT13_1_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient from the Middle Ages")
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

pdf(file = "./REAL_PTN_KIT13_1_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient from the Middle Ages",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the Middle Ages")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient from the Middle Ages")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./REAL_PTN_KIT13_1_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
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

pdf(file = "./REAL_PTN_KIT13_1_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 3500 BC (Domestication)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT13
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp <- raw_smp[which(raw_smp$age_mean <= 3500 + 2000), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1e+00
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_PTN_KIT13_2.rda")

load("./REAL_PTN_KIT13_2.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_PTN_KIT13_2_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient from the Middle Ages")
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

pdf(file = "./REAL_PTN_KIT13_2_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient from the Middle Ages",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the Middle Ages")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient from the Middle Ages")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./REAL_PTN_KIT13_2_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
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

pdf(file = "./REAL_PTN_KIT13_2_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

############################################################

#' KIT16

#' Raw data of Wutke et al. (2016) from 3500 BC
load("./Data/REAL.rda")

set.seed(2)
raw_smp <- KIT16
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 5e-01
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_PTN_KIT16_1.rda")

load("./REAL_PTN_KIT16_1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_PTN_KIT16_1_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient from the Middle Ages")
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

pdf(file = "./REAL_PTN_KIT16_1_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient from the Middle Ages",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the Middle Ages")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient from the Middle Ages")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./REAL_PTN_KIT16_1_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
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

pdf(file = "./REAL_PTN_KIT16_1_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 3500 BC (Domestication)
load("./Data/REAL.rda")

set.seed(2)
raw_smp <- KIT16
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp <- raw_smp[which(raw_smp$age_mean <= 3500 + 2000), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 5e-01
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_PTN_KIT16_2.rda")

load("./REAL_PTN_KIT16_2.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_PTN_KIT16_2_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient from the Middle Ages")
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

pdf(file = "./REAL_PTN_KIT16_2_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient from the Middle Ages",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the Middle Ages")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient from the Middle Ages")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./REAL_PTN_KIT16_2_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
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

pdf(file = "./REAL_PTN_KIT16_2_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

############################################################

#' TRPM1

#' Raw data of Wutke et al. (2016) from 14500 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- TRPM1
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 5e-01
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_PTN_TRPM1_1.rda")

load("./REAL_PTN_TRPM1_1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_PTN_TRPM1_1_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient from domestication")
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

pdf(file = "./REAL_PTN_TRPM1_1_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient from domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient from domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./REAL_PTN_TRPM1_1_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
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

pdf(file = "./REAL_PTN_TRPM1_1_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

#######################################

#' Raw data of Wutke et al. (2016) from 9322 BC (Holocene 9700 BC)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- TRPM1
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp <- raw_smp[which(raw_smp$age_mean <= 9700 + 2000), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 5e-01
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./REAL_PTN_TRPM1_2.rda")

load("./REAL_PTN_TRPM1_2.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn

pdf(file = "./REAL_PTN_TRPM1_2_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient from domestication")
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

pdf(file = "./REAL_PTN_TRPM1_2_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient from domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient from domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./REAL_PTN_TRPM1_2_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
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

pdf(file = "./REAL_PTN_TRPM1_2_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################

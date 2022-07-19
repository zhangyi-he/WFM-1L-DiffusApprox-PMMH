#' @title Estimating temporally variable selection intensity from ancient DNA data
#' @author Zhangyi He, Xiaoyang Dai, Wenyang Lyu, Mark Beaumont, Feng Yu

#' version 1.0 (PMMH)
#' Phenotypes controlled by a single gene
#' Constant natural selection and constant demographic histories

#' Input: called genotypes
#' Output: posteriors for the selection coefficient

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
source("./Code v1.0/RFUN.R")

################################################################################

#' ASIP

#' Raw data of Wutke et al. (2016) from 12500 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- ASIP
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- 0e+00
dom_par <- 0e+00
pop_siz <- tail(pop_siz, n = 1)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_COL_ASIP_1.rda")

load("./REAL_COL_ASIP_1.rda")

sel_cof_chn <- PMMH

pdf(file = "./REAL_COL_ASIP_1_Traceplot_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./REAL_COL_ASIP_1_Posterior_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 9322 BC (Holocene 9700 BC)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- ASIP
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
raw_smp <- raw_smp[which(raw_smp$age_mean <= 9700 + 2000), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- 0e+00
dom_par <- 0e+00
pop_siz <- tail(pop_siz, n = 1)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_COL_ASIP_2.rda")

load("./REAL_COL_ASIP_2.rda")

sel_cof_chn <- PMMH

pdf(file = "./REAL_COL_ASIP_2_Traceplot_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./REAL_COL_ASIP_2_Posterior_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

############################################################

#' MC1R

#' Raw data of Wutke et al. (2016) from 4300 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- MC1R
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- 0e+00
dom_par <- 0e+00
pop_siz <- tail(pop_siz, n = 1)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_COL_MC1R_1.rda")

load("./REAL_COL_MC1R_1.rda")

sel_cof_chn <- PMMH

pdf(file = "./REAL_COL_MC1R_1_Traceplot_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./REAL_COL_MC1R_1_Posterior_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 9322 BC (Holocene 9700 BC)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- MC1R
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
raw_smp <- raw_smp[which(raw_smp$age_mean <= 9700 + 2000), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- 0e+00
dom_par <- 0e+00
pop_siz <- tail(pop_siz, n = 1)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_COL_MC1R_2.rda")

load("./REAL_COL_MC1R_2.rda")

sel_cof_chn <- PMMH

pdf(file = "./REAL_COL_MC1R_2_Traceplot_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./REAL_COL_MC1R_2_Posterior_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

############################################################

#' KIT13

#' Raw data of Wutke et al. (2016) from 3645 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT13
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- 0e+00
dom_par <- 1e+00
pop_siz <- tail(pop_siz, n = 1)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_PTN_KIT13_1.rda")

load("./REAL_PTN_KIT13_1.rda")

sel_cof_chn <- PMMH

pdf(file = "./REAL_PTN_KIT13_1_Traceplot_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./REAL_PTN_KIT13_1_Posterior_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 3500 BC (Domestication)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT13
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
raw_smp <- raw_smp[which(raw_smp$age_mean <= 3500 + 2000), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- 0e+00
dom_par <- 1e+00
pop_siz <- tail(pop_siz, n = 1)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_PTN_KIT13_2.rda")

load("./REAL_PTN_KIT13_2.rda")

sel_cof_chn <- PMMH

pdf(file = "./REAL_PTN_KIT13_2_Traceplot_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./REAL_PTN_KIT13_2_Posterior_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

############################################################

#' KIT16

#' Raw data of Wutke et al. (2016) from 3500 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT16
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- 0e+00
dom_par <- 5e-01
pop_siz <- tail(pop_siz, n = 1)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_PTN_KIT16_1.rda")

load("./REAL_PTN_KIT16_1.rda")

sel_cof_chn <- PMMH

pdf(file = "./REAL_PTN_KIT16_1_Traceplot_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./REAL_PTN_KIT16_1_Posterior_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 3500 BC (Domestication)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT16
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
raw_smp <- raw_smp[which(raw_smp$age_mean <= 3500 + 2000), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- 0e+00
dom_par <- 5e-01
pop_siz <- tail(pop_siz, n = 1)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_PTN_KIT16_2.rda")

load("./REAL_PTN_KIT16_2.rda")

sel_cof_chn <- PMMH

pdf(file = "./REAL_PTN_KIT16_2_Traceplot_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./REAL_PTN_KIT16_2_Posterior_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

############################################################

#' TRPM1

#' Raw data of Wutke et al. (2016) from 14500 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- TRPM1
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- 0e+00
dom_par <- 5e-01
pop_siz <- tail(pop_siz, n = 1)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_PTN_TRPM1_1.rda")

load("./REAL_PTN_TRPM1_1.rda")

sel_cof_chn <- PMMH

pdf(file = "./REAL_PTN_TRPM1_1_Traceplot_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./REAL_PTN_TRPM1_1_Posterior_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 9322 BC (Holocene 9700 BC)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- TRPM1
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
raw_smp <- raw_smp[which(raw_smp$age_mean <= 9700 + 2000), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]
raw_smp <- raw_smp[which(rowSums(raw_smp[, c(2, 3, 4)]) != 0), ]

sel_cof <- 0e+00
dom_par <- 5e-01
pop_siz <- tail(pop_siz, n = 1)
raw_smp
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, raw_smp, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_PTN_TRPM1_2.rda")

load("./REAL_PTN_TRPM1_2.rda")

sel_cof_chn <- PMMH

pdf(file = "./REAL_PTN_TRPM1_2_Traceplot_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./REAL_PTN_TRPM1_2_Posterior_SelCoeff.pdf", width = 12, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################

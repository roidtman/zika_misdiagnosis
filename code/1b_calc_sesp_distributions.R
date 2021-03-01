# Author: Rachel Oidtman

# This script is to calculate the distributions from with the sensitivity and 
# specificity values for confirmed and suspected cases are drawn.


#=============================================================================#
# load data and libraries
#=============================================================================#
rm(list = ls())
require(mvtnorm)
set.seed(1)

# assorted values --> se & sp for confirmed cases
se_c = c(0.87, 0.82)
sp_c = c(0.95, 0.96)

# assorted values --> se & sp for confirmed cases
se_igg = c(0.82, 0.87, 0.49, 0.69, 0.94, 0.9999, 0.942, 0.902, 0.951, 0.714, 0.32, 0.54, 0.57, 0.65, 
           0.9999, 0.9999, 0.14, 0.82, 0.85, 0.29, 0.74)
sp_igg = c(0.69, 0.62, 0.99, 0.96, 0.309, 0.925, 0.993, 0.959, 0.982, 0.233, 0.97, 0.97, 0.97, 0.54,
           0.74, 0.11, 0.9999, 0.72, 0.56, 0.9999, 0.86)

# braga values --> se & sp for suspected cases
se_s = c(0.813, 1, 0.583, 0.809, 0.756, 0.286)
se_s[which(se_s == 1)] = 0.9999
sp_s = c(0.109, 0.014, 0.519, 0.580, 0.635, 0.973)


#=============================================================================#
# calculate distributions for se_c & sp_c and se_s & sp_s
# need to make sure we have num_samples samples from the posterior that 
# are ALL between 0-1
#=============================================================================#


# set number of samples in the distributions
num_samples = 1500


# samples from distribution for se_c and sp_c
diag_dist_c = rmvnorm(num_samples, c(mean(se_c), mean(sp_c)), sigma = cov(matrix(c(se_c, sp_c), ncol = 2)))
colnames(diag_dist_c) = c('se', 'sp')


# samples from distribution for se_s and sp_s
diag_dist_s = rmvnorm(num_samples, c(mean(se_s), mean(sp_s)), sigma = cov(matrix(c(se_s, sp_s), ncol = 2)))
colnames(diag_dist_s) = c('se', 'sp')


# samples from distribution for se_igg and sp_igg
diag_dist_igg = rmvnorm(3000, c(mean(se_igg), mean(sp_igg)), sigma = cov(matrix(c(se_igg, sp_igg), ncol = 2)))
colnames(diag_dist_igg) = c('se', 'sp')


# only keep random draws that are between 0-1
# diag_dist_c_keep = diag_dist_c[-which(diag_dist_c[,1] < 0),]
# diag_dist_c_keep = diag_dist_c_keep[-which(diag_dist_c[,1]>1),]
# diag_dist_c_keep = diag_dist_c[-which(diag_dist_c[,2]>1),]
# diag_dist_c_keep = diag_dist_c_keep[-which(diag_dist_c_keep[,2] < 0),]
# dim(diag_dist_c)
# diag_dist_c_keep = diag_dist_c_keep[1:1000,]
# diag_dist_c = diag_dist_c_keep
diag_dist_c = diag_dist_c[1:1000,]


diag_dist_s_keep = diag_dist_s[-which(diag_dist_s[,1] < 0),]
diag_dist_s_keep = diag_dist_s_keep[-which(diag_dist_s_keep[,1] >1 ),]
diag_dist_s_keep = diag_dist_s_keep[-which(diag_dist_s_keep[,2] > 1),]
diag_dist_s_keep = diag_dist_s_keep[-which(diag_dist_s_keep[,2] < 0), ]
dim(diag_dist_s_keep)
diag_dist_s_keep = diag_dist_s_keep[1:1000,]
diag_dist_s = diag_dist_s_keep



diag_dist_igg_keep = diag_dist_igg[-which(diag_dist_igg[,1] > 1),]
diag_dist_igg_keep = diag_dist_igg_keep[-which(diag_dist_igg[,1] < 0), ]
diag_dist_igg_keep = diag_dist_igg_keep[-which(diag_dist_igg_keep[,2] > 1),]
diag_dist_igg_keep = diag_dist_igg_keep[-which(diag_dist_igg_keep[,2] < 0), ]
diag_dist_igg_keep = diag_dist_igg_keep[1:1000,]
diag_dist_igg = diag_dist_igg_keep


#=============================================================================#
# output distributions
#=============================================================================#

save(diag_dist_c, diag_dist_s, file = '~/Dropbox/misdiagnosis/output/diagnostic_distributions.RData')
save(diag_dist_igg, file = '~/Dropbox/misdiagnosis/output/diagnostic_igg_distribution.RData')

#=============================================================================#
# visualizations
#=============================================================================#
load('~/Dropbox/misdiagnosis/output/diagnostic_distributions.RData')

library(viridis)

# prob surface for confirmed cases 
se.sp.expand.c = expand.grid(se = seq(0,1,by=0.01), sp = seq(0,1,by=0.01))
diag.mean.c = c(mean(se_c), mean(sp_c))
sigma.c = sigma = cov(matrix(c(se_c, sp_c), ncol = 2))
se.keep = numeric()
sp.keep = numeric()
prob = numeric()
for(kk in 1:nrow(se.sp.expand.c)){
  se.sp = c(se.sp.expand.c$se[kk], se.sp.expand.c$sp[kk])
  prob.kk = dmvnorm(x = se.sp, mean = diag.mean.c, sigma = sigma.c, log = TRUE)
  se.keep = c(se.keep, se.sp[1])
  sp.keep = c(sp.keep, se.sp[2])
  prob = c(prob, prob.kk)
} 

df.out.c = data.frame(se = se.keep,
                    sp = sp.keep,
                    prob = prob)

# store list.out in matrices 
pdf('../output/results_20190508//se_sp_c.pdf', height = 5, width = 6)
mat.prob.c = matrix(df.out.c$prob, nrow = sqrt(nrow(df.out.c)), ncol = sqrt(nrow(df.out.c)))
layout(1)
par(mar = c(3.5,4.2,1,1))
image(mat.prob.c, col=viridis(25), yaxt = 'n', cex.axis = 1.5)
axis(side = 2, las = 2, cex.axis = 1.5)
contour(mat.prob.c, add = TRUE, col = rgb(0,0,0,1), lwd = 1.5)
points(diag_dist_c, pch = 16)
points(se_c, sp_c, cex = 2, col = 'red', pch = 16)
mtext(side = 1, 'sensitivity', line = 2)
mtext(side = 2, 'specificity', line = 3)
dev.off()


# prob surface for suspected cases 
se.sp.expand.s = expand.grid(se = seq(0,1,by=0.01), sp = seq(0,1,by=0.01))
diag.mean.s = c(mean(se_s), mean(sp_s))
sigma.s = sigma = cov(matrix(c(se_s, sp_s), ncol = 2))
se.keep = numeric()
sp.keep = numeric()
prob = numeric()
for(kk in 1:nrow(se.sp.expand.s)){
  se.sp = c(se.sp.expand.s$se[kk], se.sp.expand.s$sp[kk])
  prob.kk = dmvnorm(x = se.sp, mean = diag.mean.s, sigma = sigma.s, log = TRUE)
  se.keep = c(se.keep, se.sp[1])
  sp.keep = c(sp.keep, se.sp[2])
  prob = c(prob, prob.kk)
} 

df.out.s = data.frame(se = se.keep,
                    sp = sp.keep,
                    prob = prob)

# store list.out in matrices 
pdf('../output/results_20190508/se_sp_s.pdf', height = 5, width = 6)
mat.prob.s = matrix(df.out.s$prob, nrow = sqrt(nrow(df.out.s)), ncol = sqrt(nrow(df.out.s)))
layout(1)
par(mar = c(3.5,4.2,1,1))
image(mat.prob.s, col=viridis(25), yaxt = 'n', cex.axis = 1.5)
axis(side = 2, las = 2, cex.axis = 1.5)
contour(mat.prob.s, add = TRUE, col = rgb(0,0,0,1), lwd = 1.5)
points(diag_dist_s, pch = 16)
points(se_s, sp_s, cex = 2, col = 'red', pch = 16)
mtext(side = 1, 'sensitivity', line = 2)
mtext(side = 2, 'specificity', line = 2)
dev.off()


# prob surface for igg-confirmed cases 
se.sp.expand.igg = expand.grid(se = seq(0,1,by=0.01), sp = seq(0,1,by=0.01))
diag.mean.igg = c(mean(se_igg), mean(sp_igg))
sigma.igg = sigma = cov(matrix(c(se_igg, sp_igg), ncol = 2))
se.keep = numeric()
sp.keep = numeric()
prob = numeric()
for(kk in 1:nrow(se.sp.expand.igg)){
  se.sp = c(se.sp.expand.igg$se[kk], se.sp.expand.igg$sp[kk])
  prob.kk = dmvnorm(x = se.sp, mean = diag.mean.igg, sigma = sigma.igg, log = TRUE)
  se.keep = c(se.keep, se.sp[1])
  sp.keep = c(sp.keep, se.sp[2])
  prob = c(prob, prob.kk)
} 

df.out.igg = data.frame(se = se.keep,
                      sp = sp.keep,
                      prob = prob)

# store list.out in matrices 
# pdf('../output/results_20190508/se_sp_s.pdf', height = 5, width = 6)
mat.prob.igg = matrix(df.out.igg$prob, nrow = sqrt(nrow(df.out.igg)), ncol = sqrt(nrow(df.out.igg)))
layout(1)
par(mar = c(3.5,4.2,1,1))
image(mat.prob.igg, col=viridis(25), yaxt = 'n', cex.axis = 1.5)
axis(side = 2, las = 2, cex.axis = 1.5)
contour(mat.prob.igg, add = TRUE, col = rgb(0,0,0,1), lwd = 1.5)
points(diag_dist_igg, pch = 16)
points(se_igg, sp_igg, cex = 2, col = 'red', pch = 16)
mtext(side = 1, 'sensitivity', line = 2)
mtext(side = 2, 'specificity', line = 2)
# dev.off()

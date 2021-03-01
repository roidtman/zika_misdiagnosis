# Author: Rachel Oidtman

# In this script we calculated the weighted average of p' and p for confirmed and 
# suspected cases through time for every country.



#=============================================================================#
# load data and libraries
#=============================================================================#
rm(list = ls())
library('vioplot')
library(plyr)
library(MASS)
set.seed(1)
# library(dplyr)

load('../output/3f_revised_p_spatial_monthly_aggregated_confirmed-igm.RData')
load('../data/processed/processed_arbo_americas.RData')
load('../output/diagnostic_distributions.RData')

'%!in%' <- function(x,y)!('%in%'(x,y))

zik_c = cbind(zik_c[,1:2], rowSums(zik_c[,3:ncol(zik_c)]))
zik_s = cbind(zik_s[,1:2], rowSums(zik_s[,3:ncol(zik_s)]))

other_c = deng_c
other_c[, 3:ncol(deng_c)] = deng_c[, 3:ncol(deng_c)] + chik_c[, 3:ncol(chik_c)]
other_c = cbind(other_c[,1:2], rowSums(other_c[,3:ncol(other_c)]))
other_s = deng_s
other_s[, 3:ncol(deng_s)] = deng_s[, 3:ncol(deng_s)] + chik_s[, 3:ncol(chik_s)]
other_s = cbind(other_s[,1:2], rowSums(other_s[,3:ncol(other_s)]))

#=============================================================================#
# aggregate data monthly 
#=============================================================================#

library(lubridate)

years_dat = other_c$Year[which(other_s$Week == 39 & other_s$Year == 2015):nrow(zik_c)]
weeks_dat = other_c$Week[which(other_s$Week == 39 & other_s$Year == 2015):nrow(zik_c)]
months_dat = lubridate::month(as.Date(paste0(years_dat, '-', weeks_dat, '-', 10), format = "%Y-%U-%u"))

other_c_month = matrix(NA, nrow = length(c(9:12, 1:12, 1:8)), ncol = ncol(other_c))
colnames(other_c_month) = colnames(other_c)
colnames(other_c_month)[2] = 'Month'
other_c_month[,2] = c(9:12, 1:12, 1:8)
other_c_month[,1] = c(rep(2015, 4), rep(2016, 12), rep(2017, 8))

other_s_month = other_c_month
zik_c_month = other_c_month
zik_s_month = other_c_month


other_c = other_c[which(other_c$Week == 39 & other_c$Year == 2015):nrow(other_c),]
other_s = other_s[which(other_s$Week == 39 & other_s$Year == 2015):nrow(other_s),]
zik_c = zik_c[which(zik_c$Week == 39 & zik_c$Year == 2015):nrow(zik_c),]
zik_s = zik_s[which(zik_s$Week == 39 & zik_s$Year == 2015):nrow(zik_s),]


for(mm in 1:nrow(other_s_month)){
  
  other_s_month[mm,3] = sum(other_s[which(other_s[,1] == other_s_month[mm,1] & months_dat == other_s_month[mm,2]),3])
  
  other_c_month[mm,3] = sum(other_c[which(other_s[,1] == other_s_month[mm,1] & months_dat == other_s_month[mm,2]),3])
  
  zik_s_month[mm,3] = sum(zik_s[which(other_s[,1] == other_s_month[mm,1] & months_dat == other_s_month[mm,2]),3])
  
  zik_c_month[mm,3] = sum(zik_c[which(other_s[,1] == other_s_month[mm,1] & months_dat == other_s_month[mm,2]),3])
}



zik_c = zik_c_month
zik_s = zik_s_month
other_c = other_c_month
other_s = other_s_month


#=============================================================================#
# functions
#=============================================================================#

# function to calculate the number of missed Zika infections
# if output is positive --> there were Zika cases missed
# if output is negative --> there were dengue/chik cases missed
calc_missed_zika = function(c_Z_in, c_O_in, p_z_in, browse = F){
  if(browse) browser()
  cases_missed_out = ((c_Z_in + c_O_in) * p_z_in) - c_Z_in
  return(cases_missed_out)
}


# function to write 95% CI
write_CI = function(param_posterior, lwr = 0.025, upr = 0.975){
  
  mat = matrix(0, ncol = ncol(param_posterior), nrow = 3)
  
  for(ii in 1:ncol(param_posterior)){
    mat[1, ii] = quantile(param_posterior[, ii], c(lwr, upr), na.rm = T)[1]
    mat[2, ii] = quantile(param_posterior[, ii], c(lwr, upr), na.rm = T)[2]  
    mat[3, ii] = median(param_posterior[, ii], na.rm = T)
  }
  
  return(mat)
  
}


# moment matching beta function
beta_mm = function(mu, var){   
  a = mu*((mu* (1 - mu) / var) - 1)
  b = a*(1 - mu) / mu
  return(data.frame(a = a, b = b))
}


# calculate prob density fxn of beta input variable
calc_dist = function(p_in, browse = F){
  out = tryCatch({
    params = fitdistr(p_in, densfun = 'beta', list(shape1 = 0.5, shape2 = 0.5))
  }, error = function(err){
    if(sum(grep('optimization failed', err$message) > 0)){
      out = NA
    }else{
      stop(err)
    }
  })
  
  if(!is.na(out)){
    d1 = dbeta(seq(0.00001, 0.99999,length.out = 1000), shape1 = out$estimate[1], shape2 = out$estimate[2])
  }else{
    out_new = beta_mm(mu = mean(p_in, na.rm = T),
                      var = var(p_in, na.rm = T))
    d1 = dbeta(seq(0.00001, 0.99999,length.out = 1000), shape1 = out_new$a, shape2 = out_new$b)
  }
  return(d1)
}


#=============================================================================#
# pooling samples of p_Zc and p_Zs
#=============================================================================#

p_overall = list()
for(tt in 1:nrow(zik_c)){
  p_overall[[tt]] = list()
}


f = '../output/4f_p_distribution_overall_confirmed-igm.RData'
if(!file.exists(f)){
  for(tt in 1:nrow(zik_c)){
    if(is.na(sum(unlist(p_Zs[[tt]])))){
      p_Zs_tmp = unlist(p_Zs[[tt]])[-which(is.na(unlist(p_Zs[[tt]])))]
    }else{
      p_Zs_tmp = p_Zs[[tt]]
    }
    if(is.na(sum(unlist(p_Zc[[tt]])))){
      p_Zc_tmp = unlist(p_Zc[[tt]])[-which(is.na(unlist(p_Zc[[tt]])))]
    }else{
      p_Zc_tmp = p_Zc[[tt]]
    }
    
    if(length(p_Zc_tmp) > 1 | length(p_Zs_tmp) > 1){
      if(length(p_Zc_tmp) > 1 & length(p_Zs_tmp) > 1){
        d_Zc = calc_dist(p_in = p_Zc_tmp)
        d_Zs = calc_dist(p_in = p_Zs_tmp)
        d_avg = d_Zc * d_Zs
        p_overall[[tt]] = sample(seq(0.00001, 0.99999,length.out = 1000), prob = d_avg / sum(d_avg), replace = T)
      }else if(length(p_Zc_tmp) <= 1){ # should this be p_Zc or p_Zc_tmp
        d_Zs = calc_dist(p_in = p_Zs_tmp)
        p_overall[[tt]] = sample(seq(0.00001, 0.99999,length.out = 1000), prob = d_Zs / sum(d_Zs), replace = T)
      }else{
        d_Zc = calc_dist(p_in = p_Zc_tmp)
        p_overall[[tt]] = sample(seq(0.00001, 0.99999,length.out = 1000), prob = d_Zc / sum(d_Zc), replace = T)
      }
      if(var(p_overall[[tt]]) == 0){
        p_overall[[tt]] = jitter(p_overall[[tt]])
      }
    }else{
      p_overall[[tt]] = NA
    }
  }
  save(p_overall, file = f)
}else{
  load(f)
}


#=============================================================================#
# plotting distributions of p through time 
#=============================================================================#

props_through_time_c = 
  sapply(1:nrow(zik_c), function(ff) sum(zik_c[ff,3:ncol(zik_c)], na.rm = T) / 
           (sum(zik_c[ff,3], na.rm = T) 
            + sum(other_c[ff,3], na.rm = T)))

props_through_time_s = 
  sapply(1:nrow(zik_s), function(ff) sum(zik_s[ff,3:ncol(zik_c)], na.rm = T) / 
           (sum(zik_s[ff,3:ncol(zik_c)], na.rm = T) 
            + sum(other_s[ff,3:ncol(zik_c)], na.rm = T)))


# pdf('~/Dropbox/misdiagnosis/output/results_final/supp_4f_revised_p_distributions_spatially_monthly_aggregated.pdf',
#     height = 6, width = 6)
# par(mar = c(2,3,3,0), oma = c(2,2,2,4))
# for(tt in 1:nrow(zik_c)){
#   if(is.na(sum(unlist(p_Zs[[tt]])))){
#     p_Zs_tmp = unlist(p_Zs[[tt]])[-which(is.na(unlist(p_Zs[[tt]])))]
#   }else{
#     p_Zs_tmp = p_Zs[[tt]]
#   }
#   if(is.na(sum(unlist(p_Zc[[tt]])))){
#     p_Zc_tmp = unlist(p_Zc[[tt]])[-which(is.na(unlist(p_Zc[[tt]])))]
#   }else{
#     p_Zc_tmp = p_Zc[[tt]]
#   }
# 
#   if(!is.na(p_overall[[tt]])){
#     layout(matrix(1:3, nrow = 3))
#     if(length(p_Zc_tmp) > 1 & length(p_Zs_tmp) > 1){
#       plot(density(p_Zc_tmp), xlim = c(0,1), main = '',
#            col = adjustcolor('#6495ED', alpha.f = 0.8), las = 1)
#       mtext(side = 3,  text = expression('Posterior distribution of p'['Z,c']))
#       polygon(density(p_Zc_tmp), col = adjustcolor('#6495ED', alpha.f = 0.5), border = NA)
#       abline(v = props_through_time_c[tt], col = adjustcolor('#6495ED', alpha.f = 0.8))
#       mtext(side = 3, outer = T, text = paste0('Month ', tt))
# 
#       plot(density(p_Zs_tmp), xlim = c(0,1), main = '',
#            col = adjustcolor('#F08080', alpha.f = 0.8), las = 1)
#       mtext(side = 3,  text = expression('Posterior distribution of p'['Z,s']))
#       polygon(density(p_Zs_tmp), col = adjustcolor('#F08080', alpha.f = 0.5), border = NA)
#       abline(v = props_through_time_s[tt], col = adjustcolor('#F08080', alpha.f = 0.8))
# 
#     }else if(length(p_Zc_tmp) <= 1){ # should this be p_Zc or p_Zc_tmp
# 
#       plot(density(p_Zs_tmp), xlim = c(0,1), main = '',
#            col = adjustcolor('#F08080', alpha.f = 0.8), las = 1)
#       mtext(side = 3,  text = expression('Posterior distribution of p'['Z,s']))
#       polygon(density(p_Zs_tmp), col = adjustcolor('#F08080', alpha.f = 0.5), border = NA)
#       abline(v = props_through_time_s[tt], col = adjustcolor('#F08080', alpha.f = 0.8))
#       mtext(side = 3, outer = T, text = paste0('Month ', tt))
#       
#     }else{
#       plot(density(p_Zc_tmp), xlim = c(0,1), main = '',
#            col = adjustcolor('#6495ED', alpha.f = 0.8), las = 1)
#       mtext(side = 3,  text = expression('Posterior distribution of p'['Z,c']))
#       polygon(density(p_Zc_tmp), col = adjustcolor('#6495ED', alpha.f = 0.5), border = NA)
#       abline(v = props_through_time_c[tt], col = adjustcolor('#6495ED', alpha.f = 0.8))
#       mtext(side = 3, outer = T, text = paste0('Month ', tt))
#       
#     }
# 
#     plot(density(p_overall[[tt]]), xlim = c(0, 1), main = '',
#          col = adjustcolor('#9370DB', alpha.f = 0.8), las = 1)
#     mtext(side = 3,  text = expression('Posterior distribution of p'['Z']))
#     polygon(density(p_overall[[tt]]), col = adjustcolor('#9370DB', alpha.f = 0.5), border = NA)
#   }
#   mtext(side = 2, outer = T, text = 'Density')
#   mtext(side = 1, outer = T, 'Proportion of Zika in total Zika, chikungunya, and dengue', line = 0.5)
# }
# dev.off()


#=============================================================================#
# plotting distributions of p relative to empirical value of p
# using violin plots
#=============================================================================#

# pdf('../output/supp/supp_4f_violin_spatially_monthly_aggregated_IgM.pdf', 
#     width = 8, height = 6)
layout(matrix(1:3, nrow = 3))

# layout(1)
par(mar = c(3.5,5,1,1), oma = c(2,0,1,0))
plot(-100, -100, xlim = c(1,nrow(zik_c)), 
     ylim = c(0, 1), xaxt = 'n',
     las = 1, xlab = '', ylab = '')
mtext(side = 2,  text = expression('Posterior distribution of p'['Z,c']), line = 2.5, cex = 0.75)

l = 1
for(tt in 1:nrow(zik_c)){
  
  if(is.na(sum(unlist(p_Zc[[tt]])))){
    p_Zc_tmp = unlist(p_Zc[[tt]][-which(is.na(p_Zc[[tt]]))])
  }else{
    p_Zc_tmp = p_Zc[[tt]]
  }
  
  if(length(p_Zc_tmp)!=0){
    vioplot(p_Zc_tmp, add = T, at= l, 
            col = adjustcolor('#6495ED', alpha.f = 0.5), 
            border = adjustcolor('#6495ED', alpha.f = 0.6), 
            rectCol = 'navy', 
            colMed = 'darkolivegreen1')
  }
  segments(y0 =  props_through_time_c[l], y1 = props_through_time_c[l], 
           x0 = l-0.25, 
           x1 = l+0.25, lwd = 2, col = 'navy')
  l = l + 1
}
axis(side = 1, at = seq(1, 24, length.out = 8), labels = c('Q4', 'Q1', 'Q2', 'Q3', 'Q4',
                                                           'Q1', 'Q2', 'Q3'))
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, 24, length.out = 8)[1], 
                                                  seq(1, 24, length.out = 8)[2])), text = 2015)
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, 24, length.out = 8)[3], 
                                                  seq(1, 24, length.out = 8)[4])), text = 2016)
mtext(side = 1, line = 2, cex = 0.75, at = seq(1, 24, length.out = 8)[7], text = 2017)

axis(side = 1, tick = T, lwd.ticks = 0, labels = NA)




# p_Z,s
plot(-100, -100, xlim = c(1,nrow(zik_c)), 
     ylim = c(0, 1), xaxt = 'n',
     las = 1, xlab = '', ylab = '')
mtext(side = 2,  text = expression('Posterior distribution of p'['Z,s']), line = 2.5, cex = 0.75)

l = 1
for(tt in 1:nrow(zik_c)){
  
  if(is.na(sum(unlist(p_Zs[[tt]])))){
    p_Zs_tmp = unlist(p_Zs[[tt]][-which(is.na(p_Zs[[tt]]))])
  }else{
    p_Zs_tmp = p_Zs[[tt]]
  }
  
  if(length(p_Zs_tmp)!=0){
    vioplot(p_Zs_tmp, add = T, at= l, 
            col = adjustcolor('#F08080', alpha.f = 0.5), 
            border = adjustcolor('#F08080', alpha.f = 0.6), 
            rectCol = 'navy', 
            colMed = '#FFE800')
  }
  segments(y0 =  props_through_time_s[l], y1 = props_through_time_s[l], 
           x0 = l-0.25, 
           x1 = l+0.25, lwd = 2, col = 'deeppink4')
  l = l + 1
}
axis(side = 1, at = seq(1, 24, length.out = 8), labels = c('Q4', 'Q1', 'Q2', 'Q3', 'Q4',
                                                           'Q1', 'Q2', 'Q3'))
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, 24, length.out = 8)[1], 
                                                  seq(1, 24, length.out = 8)[2])), text = 2015)
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, 24, length.out = 8)[3], 
                                                  seq(1, 24, length.out = 8)[4])), text = 2016)
mtext(side = 1, line = 2, cex = 0.75, at = seq(1, 24, length.out = 8)[7], text = 2017)

axis(side = 1, tick = T, lwd.ticks = 0, labels = NA)




# p_Z
plot(-100, -100, xlim = c(1,nrow(zik_c)), 
     ylim = c(0, 1), xaxt = 'n',
     las = 1, xlab = '', ylab = '')
mtext(side = 2,  text = expression('Posterior distribution of p'['Z']), line = 2.5, cex = 0.75)

l = 1
for(tt in 1:nrow(zik_c)){
  vioplot(p_overall[[tt]], add = T, at= l, 
          col = adjustcolor('#C39BD3', alpha.f = 0.8), 
          border = adjustcolor('#C39BD3', alpha.f = 0.9), 
          rectCol = 'navy', 
          colMed = '#D4E6F1')
  segments(y0 =  props_through_time_c[l], y1 = props_through_time_c[l], x0 = l-0.25, 
           x1 = l+0.25, lwd = 2, col = 'navy')
  segments(y0 =  props_through_time_s[l], y1 = props_through_time_s[l], 
           x0 = l-0.25, x1 = l+0.25, lwd = 2, col = 'deeppink4')
  l = l + 1
}

mtext(side = 1, text = 'Time', line = 3.5)
axis(side = 1, at = seq(1, 24, length.out = 8), labels = c('Q4', 'Q1', 'Q2', 'Q3', 'Q4',
                                                           'Q1', 'Q2', 'Q3'))
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, 24, length.out = 8)[1], 
                                                  seq(1, 24, length.out = 8)[2])), text = 2015)
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, 24, length.out = 8)[3], 
                                                  seq(1, 24, length.out = 8)[4])), text = 2016)
mtext(side = 1, line = 2, cex = 0.75, at = seq(1, 24, length.out = 8)[7], text = 2017)

axis(side = 1, tick = T, lwd.ticks = 0, labels = NA)

# dev.off()




#=============================================================================#
# monte carlo samples for the number of missed infections from 2014-2017
#=============================================================================#

# for one country
# take 100 samples per day 
# for each sample on each day, calculate number of missed infections

num_samples = 1000
mc_samples = matrix(NA, nrow = num_samples, ncol = nrow(zik_c))

for(tt in 1:nrow(zik_c)){
  samps = sample(p_overall[[tt]], num_samples, replace = T)
  mc_samples[,tt] = calc_missed_zika(c_Z_in = sum(zik_c[tt,2:ncol(zik_c)]) + sum(zik_c[tt,2:ncol(zik_s)]),
                                     c_O_in = sum(other_c[tt,2:ncol(other_c)]) + sum(other_s[tt,2:ncol(other_s)]),
                                     p_z_in = samps, browse = F) 
}
mc_samples[which(is.na(mc_samples))] = 0


mc_samples_sum = sapply(1:num_samples, function(ff) sum(mc_samples[ff,], na.rm = T))

summary(mc_samples_sum)
CIs = quantile(mc_samples_sum, probs = c(0.0125, 0.5, 0.975))
print(CIs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -145464  103394  200907  211009  305771  799990 

# 1.25%       50%     97.5% 
# -75511.68 200907.12 497986.41 

#=============================================================================#
# region-wide estimates of missed zika cases over time
#=============================================================================#


time_missed_zika_region = matrix(NA, nrow = dim(mc_samples)[1], 
                                 ncol = dim(mc_samples)[2])
for(tt in 1:dim(mc_samples)[2]){
  time_missed_zika_region[,tt] = sapply(1:dim(mc_samples)[1], function(ff)
    sum(mc_samples[ff, tt]))
}


# aggregated Zika
zik_region = rep(NA, dim(mc_samples)[2])
l = 1
for(tt in 1:nrow(zik_c)){
  zik_region[l] = sum(zik_c[tt, 3:ncol(zik_c)], zik_s[tt, 3:ncol(zik_c)])
  l = l + 1
}

# aggregated chik
other_region = rep(NA, dim(mc_samples)[2])
l = 1
for(tt in 1:nrow(zik_c)){
  other_region[l] = sum(other_c[tt, 3:ncol(other_c)], other_s[tt, 3:ncol(other_c)])
  l = l + 1
}


time_missed_zika_region = matrix(NA, nrow = dim(mc_samples)[1], 
                                 ncol = dim(mc_samples)[2])
for(tt in 1:dim(mc_samples)[2]){
  time_missed_zika_region[,tt] = sapply(1:dim(mc_samples)[1], function(ff)
    sum(mc_samples[ff, tt]))
}



# revised zika 
revised_zika = time_missed_zika_region
for(ii in 1:nrow(time_missed_zika_region)){
  revised_zika[ii,] = revised_zika[ii,] + zik_region
}
revised_zika_CI = write_CI(revised_zika)


not_keep = c()
for(tt in 1:ncol(mc_samples)){
  if(var(time_missed_zika_region[,tt]) == 0){
    not_keep = c(not_keep, tt)
  }
}

# pdf('../output/supp/supp_f_misdiagnosed_zika_region_spatial_monthly_aggregation_confirmed-igm.pdf', 
#     height = 6, width = 8.5)
layout(matrix(c(1,1,2,2,2), nrow = 5))
par(mar = c(1,5.75,3,1), oma = c(3.5,0,1,0))

plot(-100, -100, xlim = c(1, dim(mc_samples)[2]), ylim = c(min(time_missed_zika_region), max(time_missed_zika_region)),
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
rect(xleft = -10, xright = 200, ybottom = -1000000, ytop = 0, border = rgb(0,0,0,0.1), col = rgb(0,0,0,0.1))
axis(side = 1, at = seq(1, 24, length.out = 8), labels = c('Q4', 'Q1', 'Q2', 'Q3', 'Q4',
                                                           'Q1', 'Q2', 'Q3'))
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, 24, length.out = 8)[1], 
                                                  seq(1, 24, length.out = 8)[2])), text = 2015)
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, 24, length.out = 8)[3], 
                                                  seq(1, 24, length.out = 8)[4])), text = 2016)
mtext(side = 1, line = 2, cex = 0.75, at = seq(1, 24, length.out = 8)[7], text = 2017)


axis(side = 2, at = seq(min(time_missed_zika_region), max(time_missed_zika_region), length.out = 8), las = 1,
     labels = round_any(seq(min(time_missed_zika_region), max(time_missed_zika_region), length.out = 8) / 1000, 5))
vioplot(time_missed_zika_region[,1], time_missed_zika_region[,2], time_missed_zika_region[,3], 
        time_missed_zika_region[,4],
        time_missed_zika_region[,5], 
        time_missed_zika_region[,6],
        time_missed_zika_region[,7], time_missed_zika_region[,8],
        time_missed_zika_region[,9], 
        time_missed_zika_region[,10], time_missed_zika_region[,11], time_missed_zika_region[,12],
        time_missed_zika_region[,13], time_missed_zika_region[,14], time_missed_zika_region[,15],
        time_missed_zika_region[,16], time_missed_zika_region[,17], time_missed_zika_region[,18],
        time_missed_zika_region[,19], time_missed_zika_region[,20], time_missed_zika_region[,21],
        time_missed_zika_region[,22], time_missed_zika_region[,23], time_missed_zika_region[,24],
        col = adjustcolor('navy', alpha.f = 0.5), colMed = 'darkolivegreen1', border = 'navy', rectCol = 'navy',
        add = T, at = which(1:nrow(zik_c) %!in% not_keep))
mtext(side = 2, text = 'Zika cases misdiagnosed as dengue', line = 4.25)
mtext(side = 2, text = 'or chikungunya cases (thousands)', line = 2.75)


plot(-100, -100, xlim = c(1, dim(revised_zika_CI)[2]), ylim = c(0, max(revised_zika_CI, other_region)),
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
lines(zik_region, col = '#1e84d4', lwd = 2)
lines(other_region, col = '#adc7db', lwd = 2)
polygon(x = c(1:ncol(revised_zika_CI), rev(1:ncol(revised_zika_CI))),
        y = c(revised_zika_CI[1,], rev(revised_zika_CI[2,])), 
        col = adjustcolor('#8069d4', alpha.f = 0.5), border = NA)
lines(revised_zika_CI[3,], lwd = 3, col = 'darkolivegreen1')
axis(side = 2, at = seq(min(revised_zika_CI), 
                        max(revised_zika_CI, other_region), length.out = 8), las = 1,
     labels = round_any(seq(min(revised_zika_CI), 
                            max(revised_zika_CI, other_region), length.out = 8) / 1000, 10))
mtext(side = 2, line = 3.25, text = 'Incidence (thousands)')

legend(x =15, y =600000, legend = c('Reported Zika cases', 
                                    'Reported chikungunya + dengue cases', 
                                    'Revised Zika cases'),
       lty = c(1,1,1), lwd = c(2,2,10), col = c('#1e84d4', '#adc7db', adjustcolor('#8069d4', alpha.f=0.5)), bty = 'n')
legend(x = 15, y =600000, legend = rep('', 3), lty = c(1,1,1), lwd = c(0,0,3), 
       col = c('#1e84d4', '#adc7db', 'darkolivegreen1'), bty = 'n')
axis(side = 1, at = seq(1, 24, length.out = 8), labels = c('Q4', 'Q1', 'Q2', 'Q3', 'Q4',
                                                           'Q1', 'Q2', 'Q3'))
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, 24, length.out = 8)[1], 
                                                  seq(1, 24, length.out = 8)[2])), text = 2015)
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, 24, length.out = 8)[3], 
                                                  seq(1, 24, length.out = 8)[4])), text = 2016)
mtext(side = 1, line = 2, cex = 0.75, at = seq(1, 24, length.out = 8)[7], text = 2017)

mtext(side = 1, outer = T, text = 'Time', line = 2.25)

# dev.off()



#=============================================================================#
# applying p_Zs and p_Zc independently and then aggregating
#=============================================================================#

num_samples = 1000
mc_samples_p_Zc = matrix(NA, nrow = num_samples, ncol = nrow(zik_c))
mc_samples_p_Zs = matrix(NA, nrow = num_samples, ncol = nrow(zik_c))

for(tt in which(other_s$Week == 39 & other_s$Year == 2015):nrow(zik_c)){
  
  # if(is.na(sum(unlist(p_Zc[[tt]])))){
  #   samps = sample(unlist(p_Zc[[tt]])[-which(is.na(unlist(p_Zc[[tt]])))], num_samples, replace = T) 
  # }
  if(tt %in% c(142, 155, 156, 186)){
    samps = sample(unlist(p_Zc[[tt]])[-which(is.na(unlist(p_Zc[[tt]])))], num_samples, replace = T) 
  }else{
    print('love')
    samps = sample(p_Zc[[tt]], num_samples, replace = T)
  }
  if(sum(is.na(unlist(samps))) != 0){
    print('bye')
    if(sum(unlist(samps)[-which(is.na(unlist(samps)))]) != 0){
      mc_samples_p_Zc[,tt] = calc_missed_zika(c_Z_in = sum(zik_c[tt,2:ncol(zik_c)]),
                                              c_O_in = sum(other_c[tt,2:ncol(other_c)]),
                                              p_z_in = samps, browse = F) 
    }
  }else{
    print('hi')
    mc_samples_p_Zc[,tt] = calc_missed_zika(c_Z_in = sum(zik_c[tt,2:ncol(zik_c)]),
                                            c_O_in = sum(other_c[tt,2:ncol(other_c)]),
                                            p_z_in = samps, browse = F) 
  }
  
  
  samps = sample(p_Zs[[tt]], num_samples, replace = T)
  mc_samples_p_Zs[,tt] = calc_missed_zika(c_Z_in = sum(zik_s[tt,2:ncol(zik_c)]),
                                          c_O_in = sum(other_s[tt,2:ncol(other_c)]),
                                          p_z_in = samps, browse = F) 
  
}

mc_samples_p_Zc[which(is.na(mc_samples_p_Zc))] = 0
mc_samples_p_Zs[which(is.na(mc_samples_p_Zs))] = 0

mc_samples_sum_p_Zc = sapply(1:num_samples, function(ff) sum(mc_samples_p_Zc[ff,], na.rm = T))
mc_samples_sum_p_Zs = sapply(1:num_samples, function(ff) sum(mc_samples_p_Zs[ff,], na.rm = T))

summary(mc_samples_sum_p_Zc)
summary(mc_samples_sum_p_Zs)

mc_samples_sum_total = mc_samples_sum_p_Zc + mc_samples_sum_p_Zs
quantile(mc_samples_sum_total, probs = c(0.0125, 0.5, 0.975))

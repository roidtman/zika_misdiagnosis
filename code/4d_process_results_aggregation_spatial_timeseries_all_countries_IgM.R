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

load('../output/3d_revised_p_temporal_spatial_aggregated.RData')
load('../data/processed/processed_arbo_americas.RData')
load('../output/diagnostic_distributions.RData')

'%!in%' <- function(x,y)!('%in%'(x,y))

other_c = deng_c
other_c[, 3:ncol(deng_c)] = deng_c[, 3:ncol(deng_c)] + chik_c[, 3:ncol(chik_c)]
other_s = deng_s
other_s[, 3:ncol(deng_s)] = deng_s[, 3:ncol(deng_s)] + chik_s[, 3:ncol(chik_s)]


#=============================================================================#
# start time series in week 39 of 2015
#=============================================================================#

other_s = other_s[which(other_s$Week == 39 & other_s$Year == 2015): nrow(other_s),]
other_c = other_c[which(other_c$Week == 39 & other_c$Year == 2015): nrow(other_c),]
zik_c = zik_c[which(zik_c$Week == 39 & zik_c$Year == 2015): nrow(zik_c),]
zik_s = zik_s[which(zik_s$Week == 39 & zik_s$Year == 2015): nrow(zik_s),]


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

# p_overall = 

f = '../output/4d_p_distribution_overall.RData'
if(!file.exists(f)){
  if(is.na(sum(p_Zs))){
    p_Zs_tmp = p_Zs[-which(is.na(p_Zs))]
  }else{
    p_Zs_tmp = p_Zs
  }
  if(is.na(sum(p_Zc))){
    p_Zc_tmp = p_Zc[-which(is.na(p_Zc))]
  }else{
    p_Zc_tmp = p_Zc
  }
  if(length(p_Zc_tmp) > 1 & length(p_Zs_tmp) > 1){
    d_Zc = calc_dist(p_in = p_Zc_tmp)
    d_Zs = calc_dist(p_in = p_Zs_tmp)
    d_avg = d_Zc * d_Zs
    p_overall = sample(seq(0.00001, 0.99999,length.out = 1000), prob = d_avg / sum(d_avg), replace = T)
  }else if(length(p_Zc_tmp) <= 1){ # should this be p_Zc or p_Zc_tmp
    d_Zs = calc_dist(p_in = p_Zs_tmp)
    p_overall = sample(seq(0.00001, 0.99999,length.out = 1000), prob = d_Zs / sum(d_Zs), replace = T)
  }else{
    d_Zc = calc_dist(p_in = p_Zc_tmp)
    p_overall = sample(seq(0.00001, 0.99999,length.out = 1000), prob = d_Zc / sum(d_Zc), replace = T)
  }
  if(var(p_overall) == 0){
    p_overall = jitter(p_overall)
  }
}else{
  load(f)
}



# pdf('~/Dropbox/misdiagnosis/output/results_20190508/revised_p_distributions_cumulative_cases.pdf')
# par(mar = c(3,3,3,1), oma = c(0,0,3,1))
# for(cc in countries_keep_ind){
#   if(is.na(sum(p_Zs[[cc]]))){
#     p_Zs_tmp = p_Zs[[cc]][-which(is.na(p_Zs[[cc]]))]
#   }else{
#     p_Zs_tmp = p_Zs[[cc]]
#   }
#   if(is.na(sum(p_Zc[[cc]]))){
#     p_Zc_tmp = p_Zc[[cc]][-which(is.na(p_Zc[[cc]]))]
#   }else{
#     p_Zc_tmp = p_Zc[[cc]]
#   }
#   
#   layout(matrix(1:3, nrow = 3))
#   
#   if(length(p_Zc_tmp) > 1 & length(p_Zs_tmp) > 1){
#     plot(density(p_Zc_tmp), xlim = c(0,1), main = 'Distribution of p_Zc')
#     mtext(side = 3, outer = T, text = countries[cc])
#     plot(density(p_Zs_tmp), xlim = c(0,1), main = 'Distribution of p_Zs')
#   }else if(length(p_Zc_tmp) <= 1){ # should this be p_Zc or p_Zc_tmp
#     plot(density(p_Zs_tmp), xlim = c(0,1), main = 'Distribution of p_Zs')
#     mtext(side = 3, outer = T, text = countries[cc])
#   }else{
#     plot(density(p_Zc_tmp), xlim = c(0,1), main = 'Distribution of p_Zc')
#     mtext(side = 3, outer = T, text = countries[cc])
#   }
#   
#   plot(density(p_overall[[cc]]), xlim = c(0, 1), main = 'Overall distribution of p_Z')
# }
# 
# dev.off()
# 
# # calculate 95% CI 
# p_overall_CI = array(NA, dim = c(3, nrow(zik_c), length(countries)))
# for(cc in 1:length(countries)){
#   for(tt in 1:length(p_Zc[[cc]])){
#     p_overall_CI[,tt,cc] = write_CI(param_posterior = as.matrix(p_overall[[cc]][[tt]]))
#   }
# }



#=============================================================================#
# monte carlo samples for the number of missed infections from 2014-2017
#=============================================================================#

# for one country
# take 100 samples per day 
# for each sample on each day, calculate number of missed infections

num_samples = 1000
mc_samples = rep(NA, num_samples)

  
samps = sample(p_overall, num_samples, replace = T)
mc_samples = calc_missed_zika(c_Z_in = sum(zik_c[,-(1:2)]) + sum(zik_s[,-c(1:2)]),
                              c_O_in = sum(other_c[,-c(1:2)]) + sum(other_s[,-c(1:2)]),
                              p_z_in = samps, browse = F) 


summary(mc_samples)
quantile(mc_samples, probs = c(0.0125, 0.5, 0.975))


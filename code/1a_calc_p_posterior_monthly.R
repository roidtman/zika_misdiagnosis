# Author: Rachel Oidtman

# This code is to estimate the posterior distributions of p'_Z for the misdiagnosis
# project using the beta-binomial conjugate prior relationship. We will do this for
# every country at every time point in the Americas for which we have data.


#=============================================================================#
# load data and libraries
#=============================================================================#
setwd('~/Dropbox/misdiagnosis/code/')

rm(list = ls())
load('../data/processed/processed_arbo_americas.RData')

other_c = deng_c
other_c[, 3:ncol(deng_c)] = deng_c[, 3:ncol(deng_c)] + chik_c[, 3:ncol(chik_c)]
other_s = deng_s
other_s[, 3:ncol(deng_s)] = deng_s[, 3:ncol(deng_s)] + chik_s[, 3:ncol(chik_s)]

#=============================================================================#
# aggregating monthly data
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
  
  other_s_month[mm,3:ncol(other_s_month)] = colSums(
    other_s[which(other_s$Year == other_s_month[mm,1] & months_dat == other_s_month[mm,2]),3:ncol(other_s)])
  
  other_c_month[mm,3:ncol(other_s_month)] = colSums(
    other_c[which(other_s$Year == other_s_month[mm,1] & months_dat == other_s_month[mm,2]),3:ncol(other_s)])
  
  zik_s_month[mm,3:ncol(other_s_month)] = colSums(
    zik_s[which(other_s$Year == other_s_month[mm,1] & months_dat == other_s_month[mm,2]),3:ncol(other_s)])
  
  zik_c_month[mm,3:ncol(other_s_month)] = colSums(
    zik_c[which(other_s$Year == other_s_month[mm,1] & months_dat == other_s_month[mm,2]),3:ncol(other_s)])
}


zik_c = zik_c_month
zik_s = zik_s_month
other_c = other_c_month
other_s = other_s_month

#=============================================================================#
# aggregating monthly data
#=============================================================================#

countries = colnames(deng_c[3:ncol(deng_c)])

alpha_prior = 1
beta_prior = 1


observed_p_prime_Zc = matrix(NA, nrow = nrow(other_c), ncol = length(countries))
observed_p_prime_Zs = matrix(NA, nrow = nrow(other_s), ncol = length(countries))
l = 1
for(cc in unique(countries)){
  observed_p_prime_Zc[,l] = (zik_c[,which(colnames(zik_c) == cc)]) / 
    ((zik_c[, which(colnames(zik_c) == cc)]) + other_c[,which(colnames(other_c) == cc)])
  
  observed_p_prime_Zs[,l] = (zik_s[,which(colnames(zik_s) == cc)]) / 
    ((zik_s[, which(colnames(zik_s) == cc)]) + other_s[,which(colnames(other_s) == cc)])
  l = l + 1
}

#=============================================================================#
# functions
#=============================================================================#

# alpha_in = alpha_prior from beta prior; beta_in = beta_prior from beta prior 
# num_suc = # of successes in binomial dist; num_tri = # trials in binomial dist
# returns alpha_posterior and beta_posterior for beta distribution 
bb_conj = function(alpha_in, beta_in, num_suc, num_tri){
  alpha_out = alpha_in + num_suc
  beta_out = beta_in + num_tri - num_suc
  return(c(beta_param1 = alpha_out, beta_param2 = beta_out))
}


#=============================================================================#
# calculating p'_Zc and p'_Zs for every country at every time step
#=============================================================================#

# confirmed cases
p_prime_Zc = list()
length(p_prime_Zc) = length(unique(countries))

l = 1
for(cc in unique(countries)){
  p_prime_Zc[[l]] = data.frame(beta_param1 = rep(NA, nrow(zik_c)),
                               beta_param2 = rep(NA, nrow(zik_c)))
  for(tt in 1:nrow(zik_c)){
    tmp_out = bb_conj(alpha_in = alpha_prior, beta_in = beta_prior, 
                      num_suc = zik_c[tt, which(colnames(zik_c) == cc)], 
                      num_tri = other_c[tt, which(colnames(other_c) == cc)] + zik_c[tt, which(colnames(zik_c) == cc)])
    p_prime_Zc[[l]][tt, ] = tmp_out
    
  }
  l = l + 1
}


# suspected cases
p_prime_Zs = list()
length(p_prime_Zs) = length(unique(countries))

l = 1
for(cc in unique(countries)){
  p_prime_Zs[[l]] = data.frame(beta_param1 = rep(NA, nrow(zik_s)),
                               beta_param2 = rep(NA, nrow(zik_s)))
  for(tt in 1:nrow(zik_s)){
    tmp_out = bb_conj(alpha_in = alpha_prior, beta_in = beta_prior, 
                      num_suc = zik_s[tt, which(colnames(zik_s) == cc)], 
                      num_tri = other_s[tt, which(colnames(other_s) == cc)] + zik_s[tt, which(colnames(zik_s) == cc)])
    p_prime_Zs[[l]][tt, ] = tmp_out
    
  }
  l = l + 1
}


#=============================================================================#
# save posterior distributions
#=============================================================================#

save(p_prime_Zc, p_prime_Zs, file = '../output/p_prime_time_country_posteriors_monthly_agg.RData')



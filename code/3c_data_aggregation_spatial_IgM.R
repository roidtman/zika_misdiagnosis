# Author: Rachel Oidtman

# This script is to take our estimates of p'_Z,c and p'_Z,s and determine
# allowable values of se and sp at each time point. 


#=============================================================================#
# load data and libraries
#=============================================================================#

rm(list = ls())

load('../data/processed/processed_arbo_americas.RData')
load('../output/diagnostic_distributions.RData')


zik_c = cbind(zik_c[,1:2], rowSums(zik_c[,3:ncol(zik_c)]))
zik_s = cbind(zik_s[,1:2], rowSums(zik_s[,3:ncol(zik_s)]))

other_c = deng_c
other_c[, 3:ncol(deng_c)] = deng_c[, 3:ncol(deng_c)] + chik_c[, 3:ncol(chik_c)]
other_c = cbind(other_c[,1:2], rowSums(other_c[,3:ncol(other_c)]))
other_s = deng_s
other_s[, 3:ncol(deng_s)] = deng_s[, 3:ncol(deng_s)] + chik_s[, 3:ncol(chik_s)]
other_s = cbind(other_s[,1:2], rowSums(other_s[,3:ncol(other_s)]))

countries = colnames(deng_c[3:ncol(deng_c)])


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
# calculate p' for region-wide time series
#=============================================================================#

# observed p' for region
observed_p_prime_Zc = rep(NA, length = nrow(other_c))
observed_p_prime_Zs = rep(NA, length = nrow(other_s))


for(tt in 1:nrow(other_c)){
  observed_p_prime_Zc[tt] = zik_c[tt,3] / (sum(zik_c[tt,3], other_c[tt,3]))
  observed_p_prime_Zs[tt] = zik_s[tt,3] / (sum(zik_s[tt,3], other_s[tt,3]))
}


# beta-binomial conjugate relationship to bring uncertainty into estimates of p'
alpha_prior = 1
beta_prior = 1

bb_conj = function(alpha_in, beta_in, num_suc, num_tri){
  alpha_out = alpha_in + num_suc
  beta_out = beta_in + num_tri - num_suc
  return(c(beta_param1 = alpha_out, beta_param2 = beta_out))
}


# confirmed and suspected cases p'
p_prime_Zc = matrix(NA, nrow = nrow(other_s), ncol = 2)
p_prime_Zs = matrix(NA, nrow = nrow(other_s), ncol = 2)

for(tt in 1:nrow(zik_c)){
  p_prime_Zc[tt,] = bb_conj(alpha_in = alpha_prior, beta_in = beta_prior,
                            num_suc = zik_c[tt,3],
                            num_tri = sum(other_c[tt,3], zik_c[tt,3]))
  
  p_prime_Zs[tt,] = bb_conj(alpha_in = alpha_prior, beta_in = beta_prior,
                            num_suc = zik_s[tt,3],
                            num_tri = sum(other_s[tt,3], zik_s[tt,3]))
}


#=============================================================================#
# data structure
#=============================================================================#

# with more than 1000 samples, subsequent analyses take very long period of time. start with 10.
# number of samples for posterior draws
num_samples = 10

# draw num_sample samples from posterior of p'_Zc at each time step for each country
p_prime_Zc_samples = matrix(NA, nrow = num_samples, ncol = nrow(zik_c))
p_prime_Zs_samples = matrix(NA, nrow = num_samples, ncol = nrow(zik_c))

for(tt in 1:nrow(zik_c)){
  p_prime_Zc_samples[, tt] = rbeta(num_samples, shape1 = p_prime_Zc[tt, 1], shape2 = p_prime_Zc[tt, 2])
  p_prime_Zs_samples[, tt] = rbeta(num_samples, shape1 = p_prime_Zs[tt, 1], shape2 = p_prime_Zs[tt, 2])
}



# data structure to hold variable length estimates of p_Zc and p_Zs for each time point 
p_Zc = list()
length(p_Zc) = nrow(zik_c)

p_Zs = list()
length(p_Zs) = nrow(zik_c)


#=============================================================================#
# functions
#=============================================================================#


# se and sp constraints 1 & 2 --> 1 = se > p'; 2 = se + sp cannot equal 1
constrain_diag_1_2 = function(p_prime_in, diag_mat_in, browse = F){
  if(browse) browser()
  return(which((diag_mat_in[,1] > p_prime_in) & diag_mat_in[,1] + diag_mat_in[,2] != 1))
}


# se and sp constraint 3 --> numerator and denomiator of p must either both be 
# positive or negative to result in a postive p value
constrain_diag_3 = function(p_prime_in, diag_mat_in, browse = F){
  if(browse) browser()
  
  if(length(diag_mat_in) == 0){
    diag_mat_in = matrix(NA, ncol= 2)
  }
  
  if(length(nrow(diag_mat_in)) == 0){
    diag_mat_in = matrix(diag_mat_in, ncol = 2)
  }
  
  a = sapply(1:nrow(diag_mat_in), function(ff) p_prime_in - 1 + diag_mat_in[ff,2])
  b = sapply(1:nrow(diag_mat_in), function(ff) diag_mat_in[ff,1] - 1 + diag_mat_in[ff,2])
  
  return(which(a/b <= 1 & a/b >= 0))
}


# calculate pZ function 
calc_pZ = function(se_in, sp_in, p_prime_in, browse = F){
  if(browse) browser()
  return(
    sapply(1:length(se_in), function(ff) 
      ((p_prime_in - 1 + sp_in[ff]) / (sp_in[ff] - 1 + se_in[ff])))
  )  
}


#=============================================================================#
# calculating allowable se and sp and estimating p_Z
#=============================================================================#

# l = 1
for(tt in 1:nrow(zik_c)){
  
  # determine indices of allowable sensitivity and specificity values for confirmed and suspected cases for constraints 1 and 2
  diag_c_time_country_1_2 = lapply(1:num_samples, function(ff) 
    diag_dist_c[lapply(p_prime_Zc_samples[,tt], constrain_diag_1_2, diag_dist_c)[[ff]],])
  
  diag_s_time_country_1_2 = lapply(1:num_samples, function(ff) 
    diag_dist_s[lapply(p_prime_Zs_samples[,tt], constrain_diag_1_2, diag_dist_s)[[ff]],])
  
  
  # input allowable se and sp values for constraints 1 and 2 into constraint 3 to determine final allowable se and sp values
  # for confirmed and suspected cases
  diag_c_time_country_3 = 
    lapply(1:num_samples, function(kk) 
      matrix(diag_c_time_country_1_2[[kk]],ncol=2)[lapply(1:num_samples, 
                                                          function(ff) constrain_diag_3(p_prime_Zc_samples[ff,tt], diag_c_time_country_1_2[[ff]]))[[kk]],])
  
  diag_s_time_country_3 = 
    lapply(1:num_samples, function(kk) 
      matrix(diag_s_time_country_1_2[[kk]],ncol=2)[lapply(1:num_samples, 
                                                          function(ff) constrain_diag_3(p_prime_Zs_samples[ff,tt], diag_s_time_country_1_2[[ff]]))[[kk]],])
  
  
  # calculate estimates p_Zc
  vec_tmp = c()
  for(ii in 1:num_samples){
    out_tmp = calc_pZ(se_in = matrix(diag_c_time_country_3[[ii]], ncol = 2)[,1], 
                      sp_in = matrix(diag_c_time_country_3[[ii]], ncol = 2)[,2], 
                      p_prime_in = p_prime_Zc_samples[ii, tt])
    vec_tmp = c(vec_tmp, out_tmp) 
  }
  p_Zc[[tt]] = vec_tmp
  
  # calculate estimates p_Zs
  vec_tmp = c()
  for(ii in 1:num_samples){
    out_tmp = calc_pZ(se_in = matrix(diag_s_time_country_3[[ii]], ncol = 2)[,1], 
                      sp_in = matrix(diag_s_time_country_3[[ii]], ncol = 2)[,2], 
                      p_prime_in = p_prime_Zs_samples[ii, tt], browse = F)
    vec_tmp = c(vec_tmp, out_tmp)
  }
  p_Zs[[tt]] = vec_tmp
}


#=============================================================================#
# save output
#=============================================================================#

save(p_Zs, p_Zc, file = '../output/3f_revised_p_spatial_monthly_aggregated.RData')


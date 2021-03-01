# Author: Rachel Oidtman

# This script is to demonstrate the differences in revised estimates of p'_Z,c 
# and p'_Z,s when you calculate them at each time point versuse the aggregated 
# time series. 


#=============================================================================#
# load data and libraries
#=============================================================================#

rm(list = ls())

load('../data/processed/processed_arbo_americas.RData')
load('../output/diagnostic_distributions.RData')


other_c = deng_c
other_c[, 3:ncol(deng_c)] = deng_c[, 3:ncol(deng_c)] + chik_c[, 3:ncol(chik_c)]
other_s = deng_s
other_s[, 3:ncol(deng_s)] = deng_s[, 3:ncol(deng_s)] + chik_s[, 3:ncol(chik_s)]


countries = colnames(deng_c[3:ncol(deng_c)])


#=============================================================================#
# start time series in week 39 of 2015
#=============================================================================#

other_s = other_s[which(other_s$Week == 39 & other_s$Year == 2015): nrow(other_s),]
other_c = other_c[which(other_c$Week == 39 & other_c$Year == 2015): nrow(other_c),]
zik_c = zik_c[which(zik_c$Week == 39 & zik_c$Year == 2015): nrow(zik_c),]
zik_s = zik_s[which(zik_s$Week == 39 & zik_s$Year == 2015): nrow(zik_s),]

#=============================================================================#
# calculate p' for cumulative time series 
#=============================================================================#

# observed cumulative p' for every country
observed_p_prime_Zc = rep(NA, length = length(countries))
observed_p_prime_Zs = rep(NA, length = length(countries))

l = 1
for(cc in unique(countries)){
  observed_p_prime_Zc[l] = (sum(zik_c[,which(colnames(zik_c) == cc)])) / 
    (sum((zik_c[, which(colnames(zik_c) == cc)])) + sum(other_c[,which(colnames(other_c) == cc)]))
  
  observed_p_prime_Zs[l] = (sum(zik_s[,which(colnames(zik_s) == cc)])) / 
    (sum((zik_s[, which(colnames(zik_s) == cc)])) + sum(other_s[,which(colnames(other_s) == cc)]))
  l = l + 1
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
p_prime_Zc = matrix(NA, nrow = 2, ncol = length(countries))
p_prime_Zs = matrix(NA, nrow = 2, ncol = length(countries))

l = 1
for(cc in unique(countries)){
  p_prime_Zc[,l] = bb_conj(alpha_in = alpha_prior, beta_in = beta_prior, 
                           num_suc = sum(zik_c[, which(colnames(zik_c) == cc)]), 
                           num_tri = sum(other_c[, which(colnames(other_c) == cc)] + zik_c[, which(colnames(zik_c) == cc)]))
  
  p_prime_Zs[,l] = bb_conj(alpha_in = alpha_prior, beta_in = beta_prior, 
                           num_suc = sum(zik_s[, which(colnames(zik_c) == cc)]), 
                           num_tri = sum(other_s[, which(colnames(other_c) == cc)] + zik_s[, which(colnames(zik_c) == cc)]))
  
  l = l + 1
}


#=============================================================================#
# data structure
#=============================================================================#

# with more than 1000 samples, subsequent analyses take very long period of time. start with 10.
# number of samples for posterior draws
num_samples = 10

# draw num_sample samples from posterior of p'_Zc at each time step for each country
p_prime_Zc_samples = matrix(NA, nrow = num_samples, ncol = length(countries))
p_prime_Zs_samples = matrix(NA, nrow = num_samples, ncol = length(countries))

for(cc in 1:length(countries)){
    p_prime_Zc_samples[, cc] = rbeta(num_samples, shape1 = p_prime_Zc[1,cc], shape2 = p_prime_Zc[2,cc])
    p_prime_Zs_samples[, cc] = rbeta(num_samples, shape1 = p_prime_Zs[1,cc], shape2 = p_prime_Zs[2,cc])
}




# data structure to hold variable length estimates of p_Zc and p_Zs for each country and each time point 
p_Zc = list()
length(p_Zc) = length(countries)

p_Zs = list()
length(p_Zs) = length(countries)



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
for(cc in 1:length(countries)){
    
    # determine indices of allowable sensitivity and specificity values for confirmed and suspected cases for constraints 1 and 2
    diag_c_time_country_1_2 = lapply(1:num_samples, function(ff) 
      diag_dist_c[lapply(p_prime_Zc_samples[,cc], constrain_diag_1_2, diag_dist_c)[[ff]],])
    
    diag_s_time_country_1_2 = lapply(1:num_samples, function(ff) 
      diag_dist_s[lapply(p_prime_Zs_samples[,cc], constrain_diag_1_2, diag_dist_s)[[ff]],])
    
    
    # input allowable se and sp values for constraints 1 and 2 into constraint 3 to determine final allowable se and sp values
    # for confirmed and suspected cases
    diag_c_time_country_3 = 
      lapply(1:num_samples, function(kk) 
        matrix(diag_c_time_country_1_2[[kk]],ncol=2)[lapply(1:num_samples, 
                                                            function(ff) constrain_diag_3(p_prime_Zc_samples[ff,cc], diag_c_time_country_1_2[[ff]]))[[kk]],])
    
    diag_s_time_country_3 = 
      lapply(1:num_samples, function(kk) 
        matrix(diag_s_time_country_1_2[[kk]],ncol=2)[lapply(1:num_samples, 
                                                            function(ff) constrain_diag_3(p_prime_Zs_samples[ff,cc], diag_s_time_country_1_2[[ff]]))[[kk]],])
    
    
    # calculate estimates p_Zc
    vec_tmp = c()
    for(ii in 1:num_samples){
      out_tmp = calc_pZ(se_in = matrix(diag_c_time_country_3[[ii]], ncol = 2)[,1], 
                        sp_in = matrix(diag_c_time_country_3[[ii]], ncol = 2)[,2], 
                        p_prime_in = p_prime_Zc_samples[ii, cc])
      vec_tmp = c(vec_tmp, out_tmp) 
    }
    p_Zc[[cc]] = unlist(vec_tmp)
    
    # calculate estimates p_Zs
    vec_tmp = c()
    for(ii in 1:num_samples){
      out_tmp = calc_pZ(se_in = matrix(diag_s_time_country_3[[ii]], ncol = 2)[,1], 
                        sp_in = matrix(diag_s_time_country_3[[ii]], ncol = 2)[,2], 
                        p_prime_in = p_prime_Zs_samples[ii, cc], browse = F)
      vec_tmp = c(vec_tmp, out_tmp)
    }
    p_Zs[[cc]] = unlist(vec_tmp)
}


#=============================================================================#
# save output
#=============================================================================#

save(p_Zs, p_Zc, file = '../output/3b_revised_p_temporal_aggregated.RData')

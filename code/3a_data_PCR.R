# Author: Rachel Oidtman

# This script is to take our estimates of p'_Z,c and p'_Z,s and determine
# allowable values of se and sp at each time point. 


#=============================================================================#
# load data and libraries
#=============================================================================#

rm(list = ls())

load('../data/processed/processed_arbo_americas.RData')
load('../output/diagnostic_distributions.RData')
load('../output/diagnostic_igg_distribution.RData')
load('../output/p_prime_time_country_posteriors_monthly_agg.RData')

rm(diag_dist_c)

other_c = deng_c
other_c[, 3:ncol(deng_c)] = deng_c[, 3:ncol(deng_c)] + chik_c[, 3:ncol(chik_c)]
other_s = deng_s
other_s[, 3:ncol(deng_s)] = deng_s[, 3:ncol(deng_s)] + chik_s[, 3:ncol(chik_s)]


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
# data structure
#=============================================================================#

# with more than 1000 samples, subsequent analyses take very long period of time. start with 10.
# number of samples for posterior draws
num_samples = 10

# draw num_sample samples from posterior of p'_Zc at each time step for each country
p_prime_Zc_samples = array(NA, dim = c(num_samples, nrow(p_prime_Zc[[1]]), length(countries)))

for(cc in 1:length(countries)){
  for(tt in 1:nrow(p_prime_Zc[[1]])){
    p_prime_Zc_samples[, tt, cc] = 
      rbeta(num_samples, shape1 = p_prime_Zc[[cc]]$beta_param1[tt], shape2 = p_prime_Zc[[cc]]$beta_param2[tt])
  }
}


# draw num_sample samples from posterior of p'_Zs at each time step for each country
p_prime_Zs_samples = array(NA, dim = c(num_samples, nrow(p_prime_Zs[[1]]), length(countries)))

for(cc in 1:length(countries)){
  for(tt in 1:nrow(p_prime_Zs[[1]])){
    p_prime_Zs_samples[, tt, cc] = 
      rbeta(num_samples, shape1 = p_prime_Zs[[cc]]$beta_param1[tt], shape2 = p_prime_Zs[[cc]]$beta_param2[tt])
  }
}


# data structure to hold variable length estimates of p_Zc and p_Zs for each country and each time point 
p_Zc = list()
length(p_Zc) = length(countries)

p_Zs = list()
length(p_Zs) = length(countries)


# data structure to hold variable length revised estimates of C_Z, C_O, S_Z, and S_O
C_Z = list()
length(C_Z) = length(countries)
C_O = list()
length(C_O) = length(countries)

S_Z = list()
length(S_Z) = length(countries)
S_O = list()
length(S_O) = length(countries)


#=============================================================================#
# functions
#=============================================================================#

# # se / sp constraint function
# constrain_diag = function(p_prime_in, se_vec_in, browse = F){
#   if(browse) browser()
#   tmp = list()
#   length(tmp) = length(p_prime_in)
#     
#   for(ff in 1:length(p_prime_in)){
#     tmp[[ff]] = which(se_vec_in >= p_prime_in[ff])
#   }
#   
#   return(tmp)
# }


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
# getting allowable se / sp for different points in time for Brazil plot
#=============================================================================#

# barplot(other_s[,which(colnames(other_s) == 'COL')])
# barplot(zik_s[,which(colnames(other_s) == 'COL')], add = T)
# abline(v = times_to_plot)
# 
# cc = which(countries == 'COL')
# times_to_plot = c(1, 5, 10, 15)
# allowable_se_sp = list()
# length(allowable_se_sp) = length(times_to_plot)
# 
# l = 1
# for(tt in times_to_plot){
# 
#   diag_s_time_country_1_2 = lapply(1:num_samples, function(ff)
#     diag_dist_s[lapply(p_prime_Zs_samples[,tt,cc], constrain_diag_1_2, diag_dist_s)[[ff]],])
# 
#   diag_s_time_country_3 =
#     lapply(1:num_samples, function(kk)
#       matrix(diag_s_time_country_1_2[[kk]],ncol=2)[lapply(1:num_samples,
#                                                           function(ff) constrain_diag_3(p_prime_Zs_samples[ff,tt,cc],
#                                                                                         diag_s_time_country_1_2[[ff]]))[[kk]],])
# 
#   allowable_se_sp[[l]] = diag_s_time_country_3
# 
#   l = l + 1
# }
# 
# save(allowable_se_sp, file = '~/Dropbox/misdiagnosis/output/allowable_se_sp_plotting_colombia.RData')


#=============================================================================#
# calculating allowable se and sp and estimating p_Z
#=============================================================================#

# l = 1
for(cc in 1:length(countries)){
  for(tt in 1:nrow(p_prime_Zc[[1]])){
    
    # determine indices of allowable sensitivity and specificity values for confirmed and suspected cases for constraints 1 and 2
    diag_c_time_country_1_2 = lapply(1:num_samples, function(ff) 
      diag_dist_igg[lapply(p_prime_Zc_samples[,tt,cc], constrain_diag_1_2, diag_dist_igg)[[ff]],])
    
    diag_s_time_country_1_2 = lapply(1:num_samples, function(ff) 
      diag_dist_s[lapply(p_prime_Zs_samples[,tt,cc], constrain_diag_1_2, diag_dist_s)[[ff]],])
    
    
    # input allowable se and sp values for constraints 1 and 2 into constraint 3 to determine final allowable se and sp values
    # for confirmed and suspected cases
    diag_c_time_country_3 = 
      lapply(1:num_samples, function(kk) 
        matrix(diag_c_time_country_1_2[[kk]],ncol=2)[lapply(1:num_samples, 
                                             function(ff) constrain_diag_3(p_prime_Zc_samples[ff,tt,cc], 
                                                                           diag_c_time_country_1_2[[ff]]))[[kk]],])
    
    diag_s_time_country_3 = 
      lapply(1:num_samples, function(kk) 
        matrix(diag_s_time_country_1_2[[kk]],ncol=2)[lapply(1:num_samples, 
                                             function(ff) constrain_diag_3(p_prime_Zs_samples[ff,tt,cc], 
                                                                           diag_s_time_country_1_2[[ff]]))[[kk]],])
    
    
    # calculate estimates p_Zc
    vec_tmp = c()
    for(ii in 1:num_samples){
      out_tmp = calc_pZ(se_in = matrix(diag_c_time_country_3[[ii]], ncol = 2)[,1], 
                        sp_in = matrix(diag_c_time_country_3[[ii]], ncol = 2)[,2], 
                        p_prime_in = p_prime_Zc_samples[ii, tt, cc])
      vec_tmp = c(vec_tmp, out_tmp) 
    }
    p_Zc[[cc]][[tt]] = unlist(vec_tmp)
    
    # calculate estimates p_Zs
    vec_tmp = c()
    for(ii in 1:num_samples){
      out_tmp = calc_pZ(se_in = matrix(diag_s_time_country_3[[ii]], ncol = 2)[,1], 
                        sp_in = matrix(diag_s_time_country_3[[ii]], ncol = 2)[,2], 
                        p_prime_in = p_prime_Zs_samples[ii, tt, cc], browse = F)
      vec_tmp = c(vec_tmp, out_tmp)
    }
    p_Zs[[cc]][[tt]] = unlist(vec_tmp)
  }
}


#=============================================================================#
# revising estimates of Zika cases at each time point in each country
#=============================================================================#

# first calculate total cases at each time point
# cases_total_c = chik_c
# cases_total_c[3:ncol(chik_c)] = other_c[3:ncol(other_c)] + zik_c[3:ncol(zik_c)]
# 
# cases_total_s = chik_s
# cases_total_s[3:ncol(chik_s)] = other_s[3:ncol(other_s)] + zik_s[3:ncol(zik_s)]
# 
# l = 1
# for(cc in 3:ncol(other_c)){
#   C_Z[[l]] = list()
#   C_O[[l]] = list()
#   
#   S_Z[[l]] = list()
#   S_O[[l]] = list()
#   
#   for(tt in 1:nrow(p_prime_Zc[[1]])){
#     if(cases_total_c[tt,cc] == 0){
#       C_Z[[l]][[tt]] = 0
#       C_Z[[l]][[tt]] = 0
#     }else{
#       C_Z[[l]][[tt]] = unlist(p_Zc[[l]][[tt]]) * cases_total_c[tt,cc]
#       C_O[[l]][[tt]] = (1-unlist(p_Zc[[l]][[tt]])) * cases_total_c[tt,cc]
#     }
#     
#     if(cases_total_s[tt,cc] == 0){
#       S_Z[[l]][[tt]] = 0
#       S_O[[l]][[tt]] = 0
#     }else{
#       S_Z[[l]][[tt]] = unlist(p_Zs[[l]][[tt]]) * cases_total_s[tt,cc]
#       S_O[[l]][[tt]] = (1-unlist(p_Zs[[l]][[tt]])) * cases_total_s[tt,cc]
#     }
#   }
#   l = l + 1
# }


#=============================================================================#
# save output
#=============================================================================#

save(p_Zs, p_Zc, file = '../output/3e_revised_case_estimates_monthly_agg_confirmed-igm.RData')


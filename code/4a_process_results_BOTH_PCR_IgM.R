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


load('../output/3e_revised_case_estimates_monthly_agg.RData')
load('../data/processed/processed_arbo_americas.RData')
load('../output/diagnostic_distributions.RData')
load('../output/p_prime_time_country_posteriors_monthly_agg.RData')

'%!in%' <- function(x,y)!('%in%'(x,y))

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
# functions
#=============================================================================#

# function to calculate the number of missed Zika infections
# if output is positive --> there were Zika cases missed
# if output is negative --> there were dengue/chik cases missed
calc_missed_zika = function(c_Z_in, c_O_in, p_z_in){
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

# p_overall = p_Zc
p_overall = list()
length(p_overall) = length(countries)

f = '../output/4e_p_distribution_overall.RData'
if(!file.exists(f)){
  for(cc in 1:length(countries)){
    l = 1
    p_overall[[cc]] = list()
    
    for(tt in 1:nrow(zik_c)){
      if(is.na(sum(p_Zs[[cc]][[tt]]))){
        p_Zs_tmp = p_Zs[[cc]][[tt]][-which(is.na(p_Zs[[cc]][[tt]]))]
      }else{
        p_Zs_tmp = p_Zs[[cc]][[tt]]
      }
      if(is.na(sum(p_Zc[[cc]][[tt]]))){
        p_Zc_tmp = p_Zc[[cc]][[tt]][-which(is.na(p_Zc[[cc]][[tt]]))]
      }else{
        p_Zc_tmp = p_Zc[[cc]][[tt]]
      }
      
      p_overall[[cc]][[l]] = list()
      
      if(length(p_Zc_tmp) > 1 | length(p_Zs_tmp) > 1){
        if(length(p_Zc_tmp) > 1 & length(p_Zs_tmp) > 1){
          d_Zc = calc_dist(p_in = p_Zc_tmp)
          d_Zs = calc_dist(p_in = p_Zs_tmp)
          d_avg = d_Zc * d_Zs
          if(!sum(d_avg)==0){
            p_overall[[cc]][[l]] = sample(seq(0.00001, 0.99999,length.out = 1000), prob = d_avg / sum(d_avg), replace = T)
          }else{
            p_overall[[cc]][[l]] = rep(0, 1000)
          }
        }else if(length(p_Zc_tmp) <= 1){ # should this be p_Zc or p_Zc_tmp
          d_Zs = calc_dist(p_in = p_Zs_tmp)
          p_overall[[cc]][[l]] = sample(seq(0.00001, 0.99999,length.out = 1000), prob = d_Zs / sum(d_Zs), replace = T)
        }else{
          d_Zc = calc_dist(p_in = p_Zc_tmp)
          p_overall[[cc]][[l]] = sample(seq(0.00001, 0.99999,length.out = 1000), prob = d_Zc / sum(d_Zc), replace = T)
        }
        if(var(p_overall[[cc]][[l]]) == 0){
          p_overall[[cc]][[l]] = jitter(p_overall[[cc]][[l]])
        }
      }else{
        p_overall[[cc]][[l]] = NA
      }
      l = l + 1
    }
  }
  save(p_overall, file = f)
}else{
  load(f)
}


#=============================================================================#
# plotting samples of p_Zc and p_Zs
#=============================================================================#
props_through_time_c = as.matrix(zik_c[,3:ncol(zik_c)] / (zik_c[,3:ncol(zik_c)] + other_c[,3:ncol(zik_c)]))
props_through_time_c[which(is.na(props_through_time_c))] = 0

props_through_time_s = as.matrix(zik_s[,3:ncol(zik_s)] / (zik_s[,3:ncol(zik_s)] + other_s[,3:ncol(zik_s)]))
props_through_time_s[which(is.na(props_through_time_s))] = 0


library(viridis)
cols = plasma(nrow(zik_c))


# legend for countries
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  # dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  # axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}


#=============================================================================#
# monte carlo samples for the number of missed infections from 2014-2017
#=============================================================================#

# for one country
# take 100 samples per day 
# for each sample on each day, calculate number of missed infections
num_samples = 1000
mc_samples = array(NA, dim = c(num_samples, nrow(zik_c), length(countries)))

for(cc in 1:length(countries)){
  l = 1
  for(tt in 1:nrow(zik_c)){
    samps = sample(p_overall[[cc]][[l]], num_samples, replace = T)
    mc_samples[,tt,cc] = calc_missed_zika(c_Z_in = zik_c[tt,cc+2] + zik_s[tt,cc+2],
                                          c_O_in = other_c[tt,cc+2] + other_s[tt,cc+2],
                                          p_z_in = samps) 
    l = l + 1
  }
}
mc_samples[which(is.na(mc_samples))] = 0


#=============================================================================#
# only include output from the 32nd week of 2015 forward
#=============================================================================#


mc_samples_sum = sapply(1:num_samples, function(ff) sum(mc_samples[ff,,], na.rm = T))
length(which(mc_samples_sum > 0))


# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 223491  303094  326513  327770  351089  430640 
summary(mc_samples_sum)
CIs = quantile(mc_samples_sum, probs = c(0.0125, 0.5, 0.975))
print(CIs)

# 1.25%      50%    97.5% 
# 207596.7 266536.3 317107.7 

# sum of zik incidence prior to 2016
mc_samples_prior_2016 = sapply(1:num_samples, function(ff)
  sum(mc_samples[ff, 1 : length(which(zik_c[,1] == 2015)[1]: which(zik_c[,1] == 2016)[1]),]))
quantile(mc_samples_prior_2016, probs = c(0.0125,0.5, 0.975))

# 1.25%       50%     97.5% 
# 74595.25 113459.56 148575.71 

#=============================================================================#
# actual prop_zika, estimated prop_zika country level
#=============================================================================#

# final zika epidemic size in every country for Zika
total_zika_ctry = sapply(3:ncol(zik_c), function(ff) sum(zik_c[,ff], 
                                                         zik_s[,ff]))


# revised final size of zika epidemic across region
quantile(sum(total_zika_ctry) + mc_samples_sum, probs = c(0.0125, 0.5, 0.975))


# log total 
log_total_zika_ctry = log10(total_zika_ctry)
log_total_zika_ctry[which(is.infinite(log_total_zika_ctry))] = 0

# final epidemic size in every country for dengue/chik
total_other_ctry = sapply(3:ncol(other_c), function(ff) sum(other_c[,ff], 
                                                            other_s[,ff]))

# prop_zika at end of epidemic
end_prop = total_zika_ctry / (total_zika_ctry + total_other_ctry)

# prop_zika from MC estimates
total_missed_zika_ctry = matrix(NA, nrow = dim(mc_samples)[1], ncol = length(countries))
for(cc in 1:length(countries)){
  total_missed_zika_ctry[,cc] = sapply(1:dim(mc_samples)[1], function(ff) sum(mc_samples[ff,,cc]))
}
end_prop_samples = matrix(NA, nrow = dim(mc_samples)[1], ncol = length(countries))
for(cc in 1:length(countries)){
  end_prop_samples[,cc] = sapply(1:dim(mc_samples)[1], function(ff)
    (total_zika_ctry[cc] + total_missed_zika_ctry[ff,cc]) / 
      (total_other_ctry[cc] + total_zika_ctry[cc] + total_missed_zika_ctry[ff,cc]))
}


end_prop_cum = sum(total_zika_ctry) / sum((total_zika_ctry + total_other_ctry))
end_prop_sample_cum = sapply(1:dim(mc_samples)[1], 
                             function(ff) (sum(total_zika_ctry + total_missed_zika_ctry[ff,])) 
                             / sum((total_other_ctry + total_zika_ctry + total_missed_zika_ctry[ff,])))


sorted_prop = sort(end_prop, index.return = T)

load('../output/figures_out/fig5_IGM.RData')

# pdf('../output/revised_prop_with_barplot_PCR-IgM.pdf', height = 5, width = 8)
layout(matrix(c(1,1,2,2,2), nrow = 5))
par(mar = c(3,6,1,0.5), oma = c(2,0,0,0))
x = sapply(1:(length(countries)),
           function(ff) -0.5 + ff)
plot(-100, -100, ylim = c(0, 6), xlim = c(1, length(countries)), las =1, xaxt = 'n', 
     ylab = '', xlab = '', cex.axis = 1, yaxt = 'n')
abline(h = 1:6, col = rgb(0,0,0,0.2))
barplot(c(log_total_zika_ctry[sorted_prop$ix], log10(sum(total_zika_ctry))), ylim = c(0,6), 
        space = c(rep(0, length(log_total_zika_ctry)), 0.5),
        xaxs = 'i', add = T, xaxt = 'n', yaxt = 'n', bty = 'n', col = '#8069d4', border = 'navy')
abline(v = 7, lwd = 2, col = rgb(0,0,0,0.75), lty = 2)
# abline(v = 43, lwd = 2, col = 'coral')
axis(side = 1, at = x, 
     labels = countries[sorted_prop$ix], las = 2)
axis(side =1, at = 44, las = 2, labels = 'Region')
axis(side = 2, at = c(0:6), las = 1, labels = c(0, 10, 100, 1000, 10000, 100000, 1e6))
mtext(side = 2, 'Total Zika incidence', line = 3.5)

x = seq(0, (length(countries) * 4), by = 4)[-1]
plot(-100, -100, ylim = c(0, 1), xlim = c(5, (length(countries) * 4) + 4), las =1, 
     xaxt = 'n', ylab = '', xlab = '', cex.axis = 1)

# PCR
vioplot(
  end_prop_samples[,7], end_prop_samples[,11], end_prop_samples[,14],
  end_prop_samples[,20],
  end_prop_samples[,24], 
  end_prop_samples[,35], end_prop_samples[,39],
  end_prop_samples[,34], end_prop_samples[,3],
  end_prop_samples[,27], end_prop_samples[,22], 
  end_prop_samples[,30], end_prop_samples[,32],
  end_prop_samples[,8], end_prop_samples[,9], 
  end_prop_samples[,17], end_prop_samples[,13], 
  end_prop_samples[,6], end_prop_samples[,1], 
  end_prop_samples[,37], end_prop_samples[,16],
  end_prop_samples[,31], end_prop_samples[,12], 
  end_prop_samples[,10], end_prop_samples[,23],
  end_prop_samples[,5], end_prop_samples[,2], 
  end_prop_samples[,28], end_prop_samples[,41],
  end_prop_samples[,42], end_prop_samples[,15], 
  end_prop_samples[,19], end_prop_samples[,40], 
  end_prop_samples[,21], end_prop_samples[,38],
  end_prop_samples[,26], end_prop_samples[,25], 
  end_prop_samples[,4], end_prop_samples[,18],
  end_prop_samples[,33], end_prop_samples[,29], 
  end_prop_samples[,36], end_prop_samples[,43], 
  at = x + 0.75,
  add = T, col = '#8069d4', border = 'navy', rectCol = 'navy', 
  colMed = 'darkolivegreen1', horizontal = F)
vioplot(end_prop_sample_cum, at = x[length(x)] + 6.5, col = '#8069d4', border = 'navy', rectCol = 'navy', 
        colMed = 'darkolivegreen1', horizontal = F, add = T)


# IgM
vioplot(
  end_prop_samples_IgM[,7], end_prop_samples_IgM[,11], end_prop_samples_IgM[,14],
  end_prop_samples_IgM[,20],
  end_prop_samples_IgM[,24], 
  end_prop_samples_IgM[,35], end_prop_samples_IgM[,39],
  end_prop_samples_IgM[,34], end_prop_samples_IgM[,3],
  end_prop_samples_IgM[,27], end_prop_samples_IgM[,22], 
  end_prop_samples_IgM[,30], end_prop_samples_IgM[,32],
  end_prop_samples_IgM[,8], end_prop_samples_IgM[,9], 
  end_prop_samples_IgM[,17], end_prop_samples_IgM[,13], 
  end_prop_samples_IgM[,6], end_prop_samples_IgM[,1], 
  end_prop_samples_IgM[,37], end_prop_samples_IgM[,16],
  end_prop_samples_IgM[,31], end_prop_samples_IgM[,12], 
  end_prop_samples_IgM[,10], end_prop_samples_IgM[,23],
  end_prop_samples_IgM[,5], end_prop_samples_IgM[,2], 
  end_prop_samples_IgM[,28], end_prop_samples_IgM[,41],
  end_prop_samples_IgM[,42], end_prop_samples_IgM[,15], 
  end_prop_samples_IgM[,19], end_prop_samples_IgM[,40], 
  end_prop_samples_IgM[,21], end_prop_samples_IgM[,38],
  end_prop_samples_IgM[,26], end_prop_samples_IgM[,25], 
  end_prop_samples_IgM[,4], end_prop_samples_IgM[,18],
  end_prop_samples_IgM[,33], end_prop_samples_IgM[,29], 
  end_prop_samples_IgM[,36], end_prop_samples_IgM[,43], 
  at = x-0.75,
  add = T, col = 'mistyrose3', border = 'mistyrose4', rectCol = 'mistyrose4', 
  colMed = 'lavender', horizontal = F)
vioplot(end_prop_sample_cum_IgM, at = x[length(x)] + 5.5, col = 'mistyrose3', border = 'mistyrose4', rectCol = 'mistyrose4', 
        colMed = 'lavender', horizontal = F, add = T)

abline(v = (7 * 4 + 2), lwd = 2, col = rgb(0,0,0,0.75), lty = 2)

legend(x =1, y =1, legend = c(expression(p[Z] ~ (PCR)),
                              expression(p[Z] ~ (IgM)),
                              expression(hat(p[Z])~'')),
       pch = c(16, 16, NA), lty = c(0, 0, 1), lwd = c(0,0,2), pt.cex = 2.5, col = c('#8069d4', 'mistyrose3', 'coral'), bg = 'white')
legend(x =1.7, y =0.9999, legend = c('','',''),
       pch = c(16, 16, NA), pt.cex = 1, col = c('darkolivegreen1', 'lavender', 'white'), bty = 'n')



l = 1
for(cc in sorted_prop$ix){
  segments(x0 = x[l]-0.75, x1 = x[l]+0.75, y0 = end_prop[cc], y1 = end_prop[cc], lwd = 3, 
           col = adjustcolor('coral', alpha.f = 0.8))
  l = l + 1
}
segments(x0 = x[length(x)] + 5.25, x[length(x)] + 6.75, y0= end_prop_cum, y1 = end_prop_cum, lwd = 3, col = adjustcolor('coral', alpha.f = 0.8))
axis(side = 1, at = x, 
     labels = countries[sorted_prop$ix], las = 2)
axis(side =1, at = x[length(x)] + 6, las = 2, labels = 'Region')
mtext(side = 2, text = 'Cumulative proportion of Zika in', line = 4)
mtext(side = 2, text = 'total Zika, chik, and dengue ', line = 2.5)
mtext(side = 1, text = 'Country', line = 3.5)

# dev.off()


#=============================================================================#
# region-wide estimates of missed zika cases over time
#=============================================================================#


time_missed_zika_region = matrix(NA, nrow = dim(mc_samples)[1], 
                                 ncol = dim(mc_samples)[2])
for(tt in 1:dim(mc_samples)[2]){
  time_missed_zika_region[,tt] = sapply(1:dim(mc_samples)[1], function(ff)
    sum(mc_samples[ff, tt, ]))
}


# aggregated Zika
zik_region = rep(NA, nrow(zik_c))
l = 1
for(tt in 1: nrow(zik_c)){
  zik_region[l] = sum(zik_c[tt, 3:ncol(zik_c)], zik_s[tt, 3:ncol(zik_c)])
  l = l + 1
}

# aggregated dengue / chik
other_region = rep(NA, nrow(zik_c))
l = 1
for(tt in 1: nrow(zik_c)){
  other_region[l] = sum(other_c[tt, 3:ncol(zik_c)], other_s[tt, 3:ncol(zik_c)])
  l = l + 1
}

# revised zika 
revised_zika = time_missed_zika_region
for(ii in 1:nrow(time_missed_zika_region)){
  revised_zika[ii,] = revised_zika[ii,] + zik_region
}
revised_zika_CI = write_CI(revised_zika)

load('../output/figures_out/fig4_IGM.RData')

# pdf('../output/revised_estimates_timeseries_PCR-IgM.pdf',
#     width = 8, height = 6)
layout(matrix(c(1,1,2,2,2), nrow = 5))
par(mar = c(1,5.75,3,1), oma = c(3.5,0,1,0))

x = seq(0, ncol(time_missed_zika_region) * 2, by = 2)[-1]

plot(-100, -100, xlim = c(1, x[length(x)]), 
     ylim = c(min(time_missed_zika_region, time_missed_zika_region_IgM), max(time_missed_zika_region, time_missed_zika_region_IgM)),
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
rect(xleft = -10, xright = 200, ybottom = -100000000, ytop = 0, border = rgb(0,0,0,0.1), col = rgb(0,0,0,0.1))
axis(side = 2, at = seq(min(time_missed_zika_region, time_missed_zika_region_IgM), 
                        max(time_missed_zika_region, time_missed_zika_region_IgM), length.out = 8), las = 1,
     labels = round_any(seq(min(time_missed_zika_region, time_missed_zika_region_IgM), 
                            max(time_missed_zika_region, time_missed_zika_region_IgM), length.out = 8) / 1000, 5))
axis(side = 2, at =0, las = 1)

# PCR
vioplot(time_missed_zika_region[,1], time_missed_zika_region[,2], time_missed_zika_region[,3], 
        time_missed_zika_region[,4], time_missed_zika_region[,5], time_missed_zika_region[,6], 
        time_missed_zika_region[,7], time_missed_zika_region[,8], time_missed_zika_region[,9], 
        time_missed_zika_region[,10], time_missed_zika_region[,11], time_missed_zika_region[,12], 
        time_missed_zika_region[,13], time_missed_zika_region[,14], time_missed_zika_region[,15], 
        time_missed_zika_region[,16], time_missed_zika_region[,17], time_missed_zika_region[,18], 
        time_missed_zika_region[,19], time_missed_zika_region[,20], time_missed_zika_region[,21], 
        time_missed_zika_region[,22], time_missed_zika_region[,23], time_missed_zika_region[,24], 
        at = x - 0.35,
        col = '#8069d4', colMed = 'darkolivegreen1', border = 'navy', rectCol = 'navy',
        add = T, 
        pchMed = 20)

# IgM
vioplot(time_missed_zika_region_IgM[,1], time_missed_zika_region_IgM[,2], time_missed_zika_region_IgM[,3], 
        time_missed_zika_region_IgM[,4], time_missed_zika_region_IgM[,5], time_missed_zika_region_IgM[,6], 
        time_missed_zika_region_IgM[,7], time_missed_zika_region_IgM[,8], time_missed_zika_region_IgM[,9], 
        time_missed_zika_region_IgM[,10], time_missed_zika_region_IgM[,11], time_missed_zika_region_IgM[,12], 
        time_missed_zika_region_IgM[,13], time_missed_zika_region_IgM[,14], time_missed_zika_region_IgM[,15], 
        time_missed_zika_region_IgM[,16], time_missed_zika_region_IgM[,17], time_missed_zika_region_IgM[,18], 
        time_missed_zika_region_IgM[,19], time_missed_zika_region_IgM[,20], time_missed_zika_region_IgM[,21], 
        time_missed_zika_region_IgM[,22], time_missed_zika_region_IgM[,23], time_missed_zika_region_IgM[,24], 
        at = x + 0.35,
        col = 'mistyrose3', border = 'mistyrose4', rectCol = 'mistyrose4', 
        colMed = 'lavender',
        add = T, 
        pchMed = 20)

legend(x =30, y =220000, legend = c('Missed Zika cases (PCR)',
                                    'Missed Zika cases (IgM)'),
       pch = 16, pt.cex = 2.5, col = c('#8069d4', 'mistyrose3'), bty = 'n')
legend(x =30, y =220000, legend = c('',
                                    ''),
       pch = 16, pt.cex = 1, col = c('darkolivegreen1', 'lavender'), bty = 'n')

mtext(side = 2, text = 'Zika cases misdiagnosed as dengue', line = 4.25)
mtext(side = 2, text = 'or chikungunya cases (thousands)', line = 2.75)
axis(side = 1, at = seq(1, x[length(x)], length.out = 8), labels = c('Q4', 'Q1', 'Q2', 'Q3', 'Q4',
                                                           'Q1', 'Q2', 'Q3'))
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, x[length(x)], length.out = 8)[1], 
                                                  seq(1, x[length(x)], length.out = 8)[2])), text = 2015)
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, x[length(x)], length.out = 8)[3], 
                                                  seq(1, x[length(x)], length.out = 8)[4])), text = 2016)
mtext(side = 1, line = 2, cex = 0.75, at = seq(1, x[length(x)], length.out = 8)[7], text = 2017)


plot(-100, -100, xlim = c(1, dim(revised_zika_CI)[2]), ylim = c(0, max(revised_zika_CI, other_region)),
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
lines(zik_region, col = '#1e84d4', lwd = 2)
lines(other_region, col = '#adc7db', lwd = 2)

# IgM
polygon(x = c(1:ncol(revised_zika_CI), rev(1:ncol(revised_zika_CI))),
        y = c(revised_zika_CI_IgM[1,], rev(revised_zika_CI_IgM[2,])), 
        col = adjustcolor('mistyrose3', alpha.f = 0.8), border = 'mistyrose3')
lines(revised_zika_CI_IgM[3,], lwd = 3, col = 'lavender')


# PCR
polygon(x = c(1:ncol(revised_zika_CI), rev(1:ncol(revised_zika_CI))),
        y = c(revised_zika_CI[1,], rev(revised_zika_CI[2,])), 
        col = adjustcolor('#8069d4', alpha.f = 0.5), border = '#8069d4')
lines(revised_zika_CI[3,], lwd = 3, col = 'darkolivegreen1')

lines(zik_region, col = '#1e84d4', lwd = 2)
lines(other_region, col = '#adc7db', lwd = 2)

axis(side = 2, at = seq(min(revised_zika_CI), 
                        max(revised_zika_CI, other_region), length.out = 8), las = 1,
     labels = round_any(seq(min(revised_zika_CI), 
                            max(revised_zika_CI, other_region), length.out = 8) / 1000, 10))
mtext(side = 2, line = 3.25, text = 'Incidence (thousands)')

legend(x =15, y =600000, legend = c('Reported Zika cases', 
                                    'Reported chikungunya + dengue cases', 
                                    'Revised Zika cases (PCR)',
                                    'Revised Zika cases (IgM)'),
       lty = c(1,1,1), lwd = c(2,2,10, 10), col = c('#1e84d4', '#adc7db', adjustcolor('#8069d4', alpha.f=0.5), 
                                                adjustcolor('mistyrose3', alpha.f = 0.8)), bty = 'n')
legend(x = 15, y =600000, legend = rep('', 4), lty = c(1,1,1,1), lwd = c(0,0,3,3), 
       col = c('#1e84d4', '#adc7db', 'darkolivegreen1', 'lavender'), bty = 'n')
axis(side = 1, at = seq(1, 24, length.out = 8), labels = c('Q4', 'Q1', 'Q2', 'Q3', 'Q4',
                                                           'Q1', 'Q2', 'Q3'))
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, 24, length.out = 8)[1], 
                                                  seq(1, 24, length.out = 8)[2])), text = 2015)
mtext(side = 1, line = 2, cex = 0.75, at = mean(c(seq(1, 24, length.out = 8)[3], 
                                                  seq(1, 24, length.out = 8)[4])), text = 2016)
mtext(side = 1, line = 2, cex = 0.75, at = seq(1, 24, length.out = 8)[7], text = 2017)

mtext(side = 1, outer = T, text = 'Time', line = 2.25)
# axis(side = 1, at = c(15, 67), labels = c(2016, 2017))
# dev.off()


#=============================================================================#
# country-specific estimates of missed zika through time
#=============================================================================#

# pdf('../output/supplement/supp_revised_estimates_timeseries_by_ctry.pdf',
#     width = 8, height = 6)

for(cc in 1:length(countries)){
  # aggregated country-specific Zika
  zik_tmp = rep(NA, nrow(zik_c[which(zik_c$Year == 2015 & zik_c$Week == 39): nrow(zik_c),]))
  l = 1
  for(tt in which(zik_c$Year==2015 & zik_c$Week == 39): nrow(zik_c)){
    zik_tmp[l] = sum(zik_c[tt, which(colnames(zik_c) == countries[cc])], 
                     zik_s[tt, which(colnames(zik_c) == countries[cc])])
    l = l + 1
  }
  
  # aggregated country-specific dengue / chik
  other_tmp = rep(NA, nrow(zik_c[which(zik_c$Year == 2015 & zik_c$Week == 39): nrow(zik_c),]))
  l = 1
  for(tt in which(zik_c$Year==2015 & zik_c$Week == 39): nrow(zik_c)){
    other_tmp[l] = sum(other_c[tt, which(colnames(zik_c) == countries[cc])], 
                       other_s[tt, which(colnames(zik_c) == countries[cc])])
    l = l + 1
  }
  
  # revised zika 
  revised_zika_tmp = mc_samples[,,cc]
  for(ii in 1:nrow(revised_zika_tmp)){
    revised_zika_tmp[ii,] = revised_zika_tmp[ii,] + zik_tmp
  }
  revised_zika_CI = write_CI(revised_zika_tmp)
  
  layout(matrix(c(1,1,2,2,2), nrow = 5))
  par(mar = c(1,6.5,3,1), oma = c(2,0,0.5,0))
  plot(-100, -100, xlim = c(1, dim(mc_samples)[2]), ylim = c(min(mc_samples[,,cc], na.rm = T),
                                                             max(mc_samples[,,cc], na.rm = T)),
       xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
  mtext(side = 3, text = countries[cc])
  axis(side = 2, at = seq(min(mc_samples[,,cc], na.rm = T),
                          max(mc_samples[,,cc], na.rm = T), 
                          length.out = 5), las = 1,
       labels = round(seq(min(mc_samples[,,cc], na.rm = T),
                          max(mc_samples[,,cc], na.rm = T), 
                          length.out = 5)))
  
  vioplot(mc_samples[1:100,1,cc], mc_samples[1:100,2,cc], mc_samples[1:100,3,cc], 
          mc_samples[1:100,4,cc], mc_samples[1:100,5,cc], mc_samples[1:100,6,cc], 
          mc_samples[1:100,7,cc], mc_samples[1:100,8,cc], mc_samples[1:100,9,cc], 
          mc_samples[1:100,10,cc], mc_samples[1:100,11,cc], mc_samples[1:100,12,cc], 
          mc_samples[1:100,13,cc], mc_samples[1:100,14,cc], mc_samples[1:100,15,cc], 
          mc_samples[1:100,16,cc], mc_samples[1:100,17,cc], mc_samples[1:100,18,cc], 
          mc_samples[1:100,19,cc], mc_samples[1:100,20,cc], mc_samples[1:100,21,cc], 
          mc_samples[1:100,22,cc], mc_samples[1:100,23,cc], mc_samples[1:100,24,cc], 
          mc_samples[1:100,25,cc], mc_samples[1:100,26,cc], mc_samples[1:100,27,cc], 
          mc_samples[1:100,28,cc], mc_samples[1:100,29,cc], mc_samples[1:100,30,cc], 
          mc_samples[1:100,31,cc], 
          mc_samples[1:100,32,cc],
          mc_samples[1:100,33,cc], 
          mc_samples[1:100,34,cc], mc_samples[1:100,35,cc], mc_samples[1:100,36,cc], 
          mc_samples[1:100,37,cc], mc_samples[1:100,38,cc], mc_samples[1:100,39,cc], 
          mc_samples[1:100,40,cc], mc_samples[1:100,41,cc], mc_samples[1:100,42,cc], 
          mc_samples[1:100,43,cc], mc_samples[1:100,44,cc], mc_samples[1:100,45,cc], 
          mc_samples[1:100,46,cc], mc_samples[1:100,47,cc], mc_samples[1:100,48,cc], 
          mc_samples[1:100,49,cc], mc_samples[1:100,50,cc], mc_samples[1:100,51,cc], 
          mc_samples[1:100,52,cc], mc_samples[1:100,53,cc], mc_samples[1:100,54,cc], 
          mc_samples[1:100,55,cc], mc_samples[1:100,56,cc], mc_samples[1:100,57,cc], 
          mc_samples[1:100,58,cc], mc_samples[1:100,59,cc], mc_samples[1:100,60,cc], 
          mc_samples[1:100,61,cc], mc_samples[1:100,62,cc], mc_samples[1:100,63,cc], 
          mc_samples[1:100,64,cc], mc_samples[1:100,65,cc], mc_samples[1:100,66,cc], 
          mc_samples[1:100,67,cc], mc_samples[1:100,68,cc], mc_samples[1:100,69,cc], 
          mc_samples[1:100,70,cc], mc_samples[1:100,71,cc], mc_samples[1:100,72,cc], 
          mc_samples[1:100,73,cc], mc_samples[1:100,74,cc], mc_samples[1:100,75,cc], 
          mc_samples[1:100,76,cc],
          mc_samples[1:100,77,cc], mc_samples[1:100,78,cc], 
          mc_samples[1:100,79,cc], mc_samples[1:100,80,cc], mc_samples[1:100,81,cc], 
          mc_samples[1:100,82,cc], mc_samples[1:100,83,cc], mc_samples[1:100,84,cc], 
          mc_samples[1:100,85,cc], mc_samples[1:100,86,cc], mc_samples[1:100,87,cc], 
          mc_samples[1:100,88,cc], mc_samples[1:100,89,cc], mc_samples[1:100,90,cc], 
          mc_samples[1:100,91,cc], mc_samples[1:100,92,cc], mc_samples[1:100,93,cc], 
          mc_samples[1:100,94,cc], mc_samples[1:100,95,cc], mc_samples[1:100,96,cc], 
          mc_samples[1:100,97,cc], mc_samples[1:100,98,cc], 
          col = '#8069d4', colMed = 'darkolivegreen1', border = 'navy', rectCol = 'navy',
          add = T, 
          pchMed = 20)
  mtext(side = 2, text = 'ZIKV cases misdiagnosed as', line = 5)
  mtext(side = 2, text = 'DENV/CHIKV cases', line = 3.5)
  axis(side = 1, at = c(15, 67), labels = c(2016, 2017))
  
  
  plot(-100, -100, xlim = c(1, dim(revised_zika_CI)[2]), ylim = c(0, max(revised_zika_CI, other_tmp)),
       xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
  lines(zik_tmp, col = '#1e84d4', lwd = 2)
  lines(other_tmp, col = '#adc7db', lwd = 2)
  polygon(x = c(1:ncol(revised_zika_CI), rev(1:ncol(revised_zika_CI))),
          y = c(revised_zika_CI[1,], rev(revised_zika_CI[2,])), 
          col = adjustcolor('#8069d4', alpha.f = 0.5), border = NA)
  lines(revised_zika_CI[3,], lwd = 3, col = 'darkolivegreen1')
  axis(side = 2, at = seq(min(revised_zika_CI), 
                          max(revised_zika_CI, other_tmp), length.out = 5), las = 1,
       labels = round(seq(min(revised_zika_CI), 
                          max(revised_zika_CI, other_tmp), length.out = 5)))
  mtext(side = 2, line = 4, text = 'Incidence')
  
  legend(x = 70, y = max(revised_zika_CI, other_tmp) * 0.95, legend = c('reported ZIKV cases', 'reported CHIKV + DENV cases', 'revised ZIKV cases'),
         lty = c(1,1,1), lwd = c(2,2,10), col = c('#1e84d4', '#adc7db', adjustcolor('#8069d4', alpha.f=0.5)), bty = 'n')
  legend(x = 70, y =max(revised_zika_CI, other_tmp) * 0.95, 
         legend = rep('', 3), lty = c(1,1,1), lwd = c(0,0,3), col = c(rgb(0,0,0,0), rgb(0,0,0,0),
                                                                      'darkolivegreen1'), bty = 'n')
  axis(side = 1, at = c(15, 67), labels = c(2016, 2017))
}

# dev.off()


#=============================================================================#
# applying p_Zs and p_Zc independently and then aggregating
#=============================================================================#


num_samples = 1000
mc_samples_p_Zc = array(NA, dim = c(num_samples, nrow(zik_c), length(countries)))
mc_samples_p_Zs = array(NA, dim = c(num_samples, nrow(zik_c), length(countries)))

for(cc in 1:length(countries)){
  l = 1
  for(tt in which(zik_c$Year == 2015 & zik_c$Week == 39):nrow(zik_c)){
    samps = sample(p_Zc[[cc]][[l]], num_samples, replace = T)
    mc_samples_p_Zc[,tt,cc] = calc_missed_zika(c_Z_in = zik_c[tt,cc+2],
                                               c_O_in = other_c[tt,cc+2],
                                               p_z_in = samps)
    
    samps = sample(p_Zs[[cc]][[l]], num_samples, replace = T)
    mc_samples_p_Zs[,tt,cc] = calc_missed_zika(c_Z_in = zik_s[tt,cc+2],
                                               c_O_in = other_s[tt,cc+2],
                                               p_z_in = samps) 
    l = l + 1
  }
}
mc_samples_p_Zc[which(is.na(mc_samples_p_Zc))] = 0
mc_samples_p_Zs[which(is.na(mc_samples_p_Zs))] = 0

mc_samples_p_Zc = mc_samples_p_Zc[,which(zik_c$Year == 2015 & zik_c$Week == 39):nrow(zik_c),]
mc_samples_p_Zs = mc_samples_p_Zs[,which(zik_c$Year == 2015 & zik_c$Week == 39):nrow(zik_c),]

mc_samples_sum_p_Zc = sapply(1:num_samples, function(ff) sum(mc_samples_p_Zc[ff,,], na.rm = T))
mc_samples_sum_p_Zs = sapply(1:num_samples, function(ff) sum(mc_samples_p_Zs[ff,,], na.rm = T))

summary(mc_samples_sum_p_Zc)
summary(mc_samples_sum_p_Zs)

mc_samples_sum_total = mc_samples_sum_p_Zc + mc_samples_sum_p_Zs
quantile(mc_samples_sum_total, probs = c(0.0125, 0.5, 0.975))

# sum of zik incidence prior to 2016
mc_samples_prior_2016_p_Zc = sapply(1:num_samples, function(ff)
  sum(mc_samples_p_Zc[ff, 1 : length(which(zik_c$Year == 2015 & zik_c$Week == 39): which(zik_c$Year == 2016 & zik_c$Week == 1)),]))
quantile(mc_samples_prior_2016_p_Zc, probs = c(0.0125,0.5, 0.975))

mc_samples_prior_2016_p_Zs = sapply(1:num_samples, function(ff)
  sum(mc_samples_p_Zs[ff, 1 : length(which(zik_c$Year == 2015 & zik_c$Week == 39): which(zik_c$Year == 2016 & zik_c$Week == 1)),]))
quantile(mc_samples_prior_2016_p_Zs, probs = c(0.0125,0.5, 0.975))



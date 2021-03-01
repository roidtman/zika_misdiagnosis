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


load('../output/3b_revised_p_temporal_aggregated.RData')
load('../data/processed/processed_arbo_americas.RData')
load('../output/diagnostic_distributions.RData')
# load('../output/p_prime_time_country_posteriors.RData')

'%!in%' <- function(x,y)!('%in%'(x,y))

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
for(cc in 1:length(countries)){
  p_overall[[cc]] = list()
}

# issues with SLV and GTM
f = '../output/4b_p_distribution_overall.RData'
if(!file.exists(f)){
  for(cc in 1:length(countries)){
    if(is.na(sum(p_Zs[[cc]]))){
      p_Zs_tmp = p_Zs[[cc]][-which(is.na(p_Zs[[cc]]))]
    }else{
      p_Zs_tmp = p_Zs[[cc]]
    }
    if(is.na(sum(p_Zc[[cc]]))){
      p_Zc_tmp = p_Zc[[cc]][-which(is.na(p_Zc[[cc]]))]
    }else{
      p_Zc_tmp = p_Zc[[cc]]
    }
    
    if(length(p_Zc_tmp) <= 1 & length(p_Zs_tmp) <= 1){
      p_overall[[cc]] = NA
    }else{
      if(length(p_Zc_tmp) > 1 & length(p_Zs_tmp) > 1){
        d_Zc = calc_dist(p_in = p_Zc_tmp)
        d_Zs = calc_dist(p_in = p_Zs_tmp)
        d_avg = d_Zc * d_Zs
        p_overall[[cc]] = sample(seq(0.00001, 0.99999,length.out = 1000), prob = d_avg / sum(d_avg), replace = T)
      }else if(length(p_Zc_tmp) <= 1){ 
        d_Zs = calc_dist(p_in = p_Zs_tmp)
        p_overall[[cc]] = sample(seq(0.00001, 0.99999,length.out = 1000), prob = d_Zs / sum(d_Zs), replace = T)
      }else{
        d_Zc = calc_dist(p_in = p_Zc_tmp)
        p_overall[[cc]] = sample(seq(0.00001, 0.99999,length.out = 1000), prob = d_Zc / sum(d_Zc), replace = T)
      }
      if(var(p_overall[[cc]]) == 0){
        p_overall[[cc]] = jitter(p_overall[[cc]])
      }
    }
  }
  save(p_overall, file = f)
}else{
  load(f)
}


#=============================================================================#
# plotting distributions of p relative to empirical value of p
#=============================================================================#

# final zika epidemic size in every country for Zika
zika_ctry_c = sapply(3:ncol(zik_c), function(ff) sum(zik_c[,ff]))
zika_ctry_s = sapply(3:ncol(zik_c), function(ff) sum(zik_s[,ff]))

# final epidemic size in every country for dengue/chik
other_ctry_c = sapply(3:ncol(other_c), function(ff) sum(other_c[,ff]))
other_ctry_s = sapply(3:ncol(other_c), function(ff) sum(other_s[,ff]))

# empirical prop_zika at end of epidemic
end_prop_c = zika_ctry_c / (zika_ctry_c + other_ctry_c)
end_prop_s = zika_ctry_s / (zika_ctry_s + other_ctry_s)


# pdf('../output/supp/supp_4b_revised_p_distributions_cumulative_cases.pdf',
#     height = 6, width = 6)
countries_keep = 1:length(countries)
countries_keep = countries_keep[-c(20, 35)]
par(mar = c(2,3,3,0), oma = c(2,2,2,4))
for(cc in countries_keep){
  if(is.na(sum(p_Zs[[cc]]))){
    p_Zs_tmp = p_Zs[[cc]][-which(is.na(p_Zs[[cc]]))]
  }else{
    p_Zs_tmp = p_Zs[[cc]]
  }
  if(is.na(sum(p_Zc[[cc]]))){
    p_Zc_tmp = p_Zc[[cc]][-which(is.na(p_Zc[[cc]]))]
  }else{
    p_Zc_tmp = p_Zc[[cc]]
  }

  layout(matrix(1:3, nrow = 3))

  if(length(p_Zc_tmp) > 1 & length(p_Zs_tmp) > 1){
    plot(density(p_Zc_tmp), xlim = c(0,1), main = '',
         col = adjustcolor('#6495ED', alpha.f = 0.8), las = 1)
    mtext(side = 3,  text = expression('Posterior distribution of p'['Z,c']))
    polygon(density(p_Zc_tmp), col = adjustcolor('#6495ED', alpha.f = 0.5), border = NA)
    mtext(side = 3, line = 2, text = countries[cc])
    abline(v = end_prop_c[cc], col = adjustcolor('#6495ED', alpha.f = 0.8))

    plot(density(p_Zs_tmp), xlim = c(0,1), main = '',
         col = adjustcolor('#F08080', alpha.f = 0.8), las = 1)
    mtext(side = 3,  text = expression('Posterior distribution of p'['Z,s']))
    polygon(density(p_Zs_tmp), col = adjustcolor('#F08080', alpha.f = 0.5), border = NA)
    abline(v = end_prop_s[cc], col = adjustcolor('#F08080', alpha.f = 0.8))

  }else if(length(p_Zc_tmp) <= 1){ # should this be p_Zc or p_Zc_tmp
    plot(density(p_Zs_tmp), xlim = c(0,1), main = '',
         col = adjustcolor('#F08080', alpha.f = 0.8), las = 1)
    mtext(side = 3,  text = expression('Posterior distribution of p'['Z,s']))
    polygon(density(p_Zs_tmp), col = adjustcolor('#F08080', alpha.f = 0.5), border = NA)
    abline(v = end_prop_s[cc], col = adjustcolor('#F08080', alpha.f = 0.8))
    mtext(side = 3, line = 2, text = countries[cc])

  }else{
    plot(density(p_Zc_tmp), xlim = c(0,1), main = '',
         col = adjustcolor('#6495ED', alpha.f = 0.8), las = 1)
    mtext(side = 3,  text = expression('Posterior distribution of p'['Z,c']))
    polygon(density(p_Zc_tmp), col = adjustcolor('#6495ED', alpha.f = 0.5), border = NA)
    abline(v = end_prop_c[cc], col = adjustcolor('#6495ED', alpha.f = 0.8))
    mtext(side = 3, line = 2, text = countries[cc])

  }

  plot(density(p_overall[[cc]]), xlim = c(0, 1), main = '',
       col = adjustcolor('#9370DB', alpha.f = 0.8), las = 1)
  mtext(side = 2, outer = T, text = 'Density')
  mtext(side = 1, outer = T, 'Proportion of Zika in total Zika, chikungunya, and dengue', line = 0.5)
  mtext(side = 3,  text = expression('Posterior distribution of p'['Z']))
  polygon(density(p_overall[[cc]]), col = adjustcolor('#9370DB', alpha.f = 0.5), border = NA)
}
# dev.off()



#=============================================================================#
# plotting distributions of p relative to empirical value of p
# using violin plots
#=============================================================================#

# pdf('../output/supp/supp_4b_violin_temporally_aggregated.pdf', width = 8, height = 6)
layout(matrix(1:3, nrow = 3))

# layout(1)
par(mar = c(3.5,5,1,1), oma = c(2,0,1,0))
plot(-100, -100, xlim = c(2, length(countries)-1), ylim = c(0, 1), xaxt = 'n',
     las = 1, xlab = '', ylab = '')
mtext(side = 2,  text = expression('Posterior distribution of p'['Z,c']), line = 2.5, cex = 0.75)
axis(side = 1, at = seq(1, length(countries), length.out = 43),
     labels = countries, las = 2)
abline(v = 1:length(countries) + 0.5, col = rgb(0,0,0,0.4))
for(cc in countries_keep){
  
  if(is.na(sum(p_Zc[[cc]]))){
    p_Zc_tmp = p_Zc[[cc]][-which(is.na(p_Zc[[cc]]))]
  }else{
    p_Zc_tmp = p_Zc[[cc]]
  }
  
  if(length(p_Zc_tmp)!=0){
    vioplot(p_Zc_tmp, add = T, at= cc, 
            col = adjustcolor('#6495ED', alpha.f = 0.5), 
            border = adjustcolor('#6495ED', alpha.f = 0.6), 
            rectCol = 'navy', 
            colMed = 'darkolivegreen1')
  }
  segments(y0 =  end_prop_c[cc], y1 = end_prop_c[cc], x0 = cc-0.25, 
           x1 = cc+0.25, lwd = 3, col = 'navy')
}


# p_Z,s
plot(-100, -100, xlim = c(2, length(countries)-1), ylim = c(0, 1), xaxt = 'n',
     las = 1, xlab = '', ylab = '')
mtext(side = 2,  text = expression('Posterior distribution of p'['Z,s']), line = 2.5, cex = 0.75)
axis(side = 1, at = seq(1, length(countries), length.out = 43),
     labels = countries, las = 2)

abline(v = 1:length(countries) + 0.5, col = rgb(0,0,0,0.4))
for(cc in countries_keep){
  
  if(is.na(sum(p_Zs[[cc]]))){
    p_Zs_tmp = p_Zs[[cc]][-which(is.na(p_Zs[[cc]]))]
  }else{
    p_Zs_tmp = p_Zs[[cc]]
  }
  
  if(length(p_Zs_tmp)!=0){
    vioplot(p_Zs_tmp, add = T, at= cc, 
            col = adjustcolor('#F08080', alpha.f = 0.5), 
            border = adjustcolor('#F08080', alpha.f = 0.6), 
            rectCol = 'navy', 
            colMed = '#FFE800')
  }
  segments(y0 =  end_prop_s[cc], y1 = end_prop_s[cc], 
           x0 = cc-0.25, x1 = cc+0.25, lwd = 3, col = 'deeppink4')
}


# p_Z
plot(-100, -100, xlim = c(2, length(countries)-1), ylim = c(0, 1), xaxt = 'n',
     las = 1, xlab = '', ylab = '')
mtext(side = 2,  text = expression('Posterior distribution of p'['Z']), line = 2.5, cex = 0.75)
axis(side = 1, at = seq(1, length(countries), length.out = 43),
     labels = countries, las = 2)

abline(v = 1:length(countries) + 0.5, col = rgb(0,0,0,0.4))
for(cc in countries_keep){
  vioplot(p_overall[[cc]], add = T, at= cc, 
          col = adjustcolor('#C39BD3', alpha.f = 0.8), 
          border = adjustcolor('#C39BD3', alpha.f = 0.9), 
          rectCol = 'navy', 
          colMed = '#D4E6F1')
  segments(y0 =  end_prop_c[cc], y1 = end_prop_c[cc], x0 = cc-0.25, 
           x1 = cc+0.25, lwd = 3, col = 'navy')
  segments(y0 =  end_prop_s[cc], y1 = end_prop_s[cc], 
           x0 = cc-0.25, x1 = cc+0.25, lwd = 3, col = 'deeppink4')
}

mtext(side = 1, text = 'Countries', line = 3.5)

# dev.off()





#=============================================================================#
# monte carlo samples for the number of missed infections from 2014-2017
#=============================================================================#

# for one country
# take 100 samples per day 
# for each sample on each day, calculate number of missed infections

num_samples = 1000
mc_samples = matrix(NA, nrow = num_samples, ncol = length(countries))

for(cc in 1:length(countries)){
 
  samps = sample(p_overall[[cc]], num_samples, replace = T)
  mc_samples[,cc] = calc_missed_zika(c_Z_in = sum(zik_c[,cc+2]) + sum(zik_s[,cc+2]),
                                    c_O_in = sum(other_c[,cc+2]) + sum(other_s[,cc+2]),
                                    p_z_in = samps, browse = F) 
}


# for(cc in 1:length(countries)){
#   if(!is.na(sum(mc_samples[,cc]))){
#     plot(density(mc_samples[,cc]), main = countries[cc])
#   }
# }


mc_samples_sum = sapply(1:num_samples, function(ff) sum(mc_samples[ff,], na.rm = T))

summary(mc_samples_sum)
quantile(mc_samples_sum, probs = c(0.0125, 0.5, 0.975))
# 1.25%       50%     97.5% 
# -532293.0 -439910.1 -287126.7 old se/sp 
# 1.25%       50%     97.5% 
# 45007.97 165209.24 277880.91 

#=============================================================================#
# actual prop_zika, estimated prop_zika country level
#=============================================================================#

# final zika epidemic size in every country for Zika
total_zika_ctry = sapply(3:ncol(zik_c), function(ff) sum(zik_c[,ff], zik_s[,ff]))

# log total 
log_total_zika_ctry = log10(total_zika_ctry)
log_total_zika_ctry[which(is.infinite(log_total_zika_ctry))] = 0


# final epidemic size in every country for dengue/chik
total_other_ctry = sapply(3:ncol(other_c), function(ff) sum(other_c[,ff], other_s[,ff]))

# empirical prop_zika at end of epidemic
end_prop = total_zika_ctry / (total_zika_ctry + total_other_ctry)

# samples of prop zika at end of epidemic
end_prop_samples = matrix(NA, nrow = dim(mc_samples)[1], ncol = length(countries))
for(cc in 1:length(countries)){
  end_prop_samples[,cc] = sapply(1:dim(mc_samples)[1], function(ff)
    (total_zika_ctry[cc] + mc_samples[ff,cc]) / (total_other_ctry[cc] + total_zika_ctry[cc]))
  if(length(which(is.na(end_prop_samples[,cc]))) == dim(end_prop_samples)[1]){
    end_prop_samples[,cc] = rep(0, dim(end_prop_samples)[1])
  }
}

# cumulative proportion over all countries
end_prop_cum = sum(total_zika_ctry) / sum((total_zika_ctry + total_other_ctry))
end_prop_sample_cum = sapply(1:dim(mc_samples)[1], 
                             function(ff) (sum(total_zika_ctry + mc_samples[ff,], na.rm = T)) 
                             / sum(total_other_ctry + total_zika_ctry + mc_samples[ff,], na.rm = T))



# ordered by end_prop
# pdf('../output/supp/supp_b_revised_prop_z_cumulative_cases.pdf', height = 5, width = 8)
sorted_prop = sort(end_prop, index.return = T)
x = sapply(1:(length(countries)), 
           function(ff) -0.25 + ff)


layout(matrix(c(1,1,2,2,2), nrow = 5))
par(mar = c(3,6,1,0.5), oma = c(2,0,0,0))
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


plot(-100, -100, ylim = c(0, 1), xlim = c(1, length(countries)), las =1, 
     xaxt = 'n', ylab = '', xlab = '', cex.axis = 1)
vioplot(
  end_prop_samples[,7], end_prop_samples[,11], end_prop_samples[,14],
  end_prop_samples[1:100,20],
  end_prop_samples[,24], 
  end_prop_samples[1:100,35], end_prop_samples[,39],
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
  at = x,
  add = T, col = '#8069d4', border = 'navy', rectCol = 'navy', 
  colMed = 'darkolivegreen1', horizontal = F)

vioplot(end_prop_sample_cum, at = 44, col = '#8069d4', border = 'navy', rectCol = 'navy', 
        colMed = 'darkolivegreen1', horizontal = F, add = T)

l = 1
for(cc in sorted_prop$ix[1:length(sorted_prop$ix)]){
  # segments(x0 = 0.5 +cc, x1 = 1+cc, y0 = end_prop[cc], y1 = end_prop[cc])
  # segments(x0 = x[cc]-0.25, x1 = x[cc]+0.25, y0 = end_prop[cc], y1 = end_prop[cc], lwd = 3, col = 'magenta')
  # segments(x0 = x[l]-0.25, x1 = x[l]+0.25, y0 = end_prop[cc], y1 = end_prop[cc], lwd = 3, col = 'coral')
  segments(x0 = x[l]-0.25, x1 = x[l]+0.25, y0 = end_prop[cc], y1 = end_prop[cc], lwd = 3, col = 'coral')
  l = l + 1
}
segments(x0 = 44-0.25, 44+0.24, y0= end_prop_cum, y1 = end_prop_cum, lwd = 3, col = adjustcolor('coral', alpha.f = 0.8))

axis(side = 1, at = x, 
     labels = countries[sorted_prop$ix], las = 2)
axis(side =1, at = 44, las = 2, labels = 'Region')
mtext(side = 2, text = 'Cumulative proportion of Zika in', line = 3.5)
mtext(side = 2, text = 'total Zika, chik, and dengue', line = 2.4)
# mtext(side = 3, text = 'Cumulative cases only', outer = T)
# dev.off()




#=============================================================================#
# applying p_Zs and p_Zc independently and then aggregating
#=============================================================================#

num_samples = 1000
mc_samples_p_Zc = matrix(NA, nrow = num_samples, ncol = length(countries))
mc_samples_p_Zs = matrix(NA, nrow = num_samples, ncol = length(countries))

for(cc in 1:length(countries)){
  
  samps = sample(p_Zc[[cc]], num_samples, replace = T)
  mc_samples_p_Zc[,cc] = calc_missed_zika(c_Z_in = sum(zik_c[,cc+2]),
                                     c_O_in = sum(other_c[,cc+2]),
                                     p_z_in = samps, browse = F) 
  
  samps = sample(p_Zs[[cc]], num_samples, replace = T)
  mc_samples_p_Zs[,cc] = calc_missed_zika(c_Z_in = sum(zik_s[,cc+2]),
                                          c_O_in = sum(other_s[,cc+2]),
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

# Processing PAHO DENV data
#=============================================================================#

#=============================================================================#
# setting up workspace
#=============================================================================#
rm(list = ls())
setwd('.')
require(tidyverse)
require(stringr)

#=============================================================================#
# loading in data
#=============================================================================#

dat_inc_c = read.csv('../data/raw/zika_report_timeseries_confirmed.csv')
dat_inc_s = read.csv('../data/raw/zika_report_timeseries_suspected.csv')


for(cc in 3:ncol(dat_inc_c)){
  if(sum(dat_inc_s[,cc] != 0)){
    dat_inc_s[,cc] = dat_inc_s[,cc] - dat_inc_c[,cc]
  }
}


layout(matrix(1:55, nrow = 5))
for(ii in 3:ncol(dat_inc_c)){
  plot(dat_inc_c[,ii], type = 'l', col = 'red', ylim = c(0, max(dat_inc_s[,ii], dat_inc_c[,ii])),
       xaxt = 'n', yaxt = 'n', main = colnames(dat_inc_c)[ii])
  lines(dat_inc_s[,ii])
}


write.csv(dat_inc_c, '../data/processed/zika_confirmed_timeseries.csv', 
          row.names = FALSE)
write.csv(dat_inc_s, '../data/processed/zika_suspected_timeseries.csv', 
          row.names = FALSE)
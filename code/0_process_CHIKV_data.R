# Processing CHIKV data

#=============================================================================#
# setting up workspace
#=============================================================================#
rm(list = ls())
setwd('.')
require(tidyverse)


#=============================================================================#
# loading in data
#=============================================================================#

dat_cum_c = read.csv(file = '../data/raw/Chikungunya_confirmed_cases.csv')
dat_cum_s = read.csv(file = '../data/raw/Chikungunya_suspected_cases.csv')

# get rid of 2013 data
dat_cum_c = data.frame(dat_cum_c)
dat_cum_c = dat_cum_c[-which(dat_cum_c$Year==2013),]
dat_cum_s = data.frame(dat_cum_s)
dat_cum_s = dat_cum_s[-which(dat_cum_s$Year==2013),]

# extract ISO country names
countries = colnames(dat_cum_c)[-c(1:2)]


#=============================================================================#
# converting the data from cumulative incidence to weekly incidence for 
# suspected and confirmed cases
#=============================================================================#

## confirmed cases
dat_inc_c = dat_cum_c[-nrow(dat_cum_c), ]
dat_inc_c[3:ncol(dat_inc_c)] = NA

for(rr in 1:(nrow(dat_inc_c) - 1)){
  for(cc in 3:ncol(dat_inc_c)){
    dat_inc_c[rr, cc] = dat_cum_c[rr + 1, cc] - dat_cum_c[rr, cc]
    if(dat_inc_c[rr, cc] < 0){
      dat_inc_c[rr, cc] = 0
    }
  }
}

write.csv(dat_inc_c, '../data/processed/chikungunya_confirmed_timeseries.csv', 
          row.names = FALSE)

layout(matrix(1:54,9,6))
par(mar = c(1,1,0.5,0.5), oma = c(2,2,1,1))

for(ii in 3:ncol(dat_cum_c)){
  plot(dat_inc_c[,ii], type = 'l', xaxt = 'n')
}


## suspected cases
dat_inc_s = dat_cum_s[-nrow(dat_cum_s), ]
dat_inc_s[3:ncol(dat_inc_s)] = NA

for(rr in 1:(nrow(dat_inc_s) - 1)){
  for(cc in 3:ncol(dat_inc_s)){
    dat_inc_s[rr, cc] = dat_cum_s[rr + 1, cc] - dat_cum_s[rr, cc]
    if(dat_inc_s[rr, cc] < 0){
      dat_inc_s[rr, cc] = 0
    }
  }
}

write.csv(dat_inc_s, '../data/processed/chikungunya_suspected_timeseries.csv', 
          row.names = FALSE)



layout(matrix(1:54,9,6))
par(mar = c(1,1,0.5,0.5), oma = c(2,2,1,1))

for(ii in 3:ncol(dat_cum_s)){
  plot(dat_inc_s[,ii], type = 'l', xaxt = 'n')
}

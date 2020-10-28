# Processing CHIKV, ZIKV, and DENV data together
#=============================================================================#

#=============================================================================#
# setting up workspace
#=============================================================================#
rm(list = ls())
setwd('.')
require(tidyverse)


#=============================================================================#
# loading in data
#=============================================================================#

chik_c = read_csv(file = '../data/processed/chikungunya_confirmed_timeseries.csv')
chik_s = read_csv(file = '../data/processed/chikungunya_suspected_timeseries.csv')
zik_c = read_csv(file = '../data/processed/zika_confirmed_timeseries.csv')
zik_s = read_csv(file = '../data/processed/zika_suspected_timeseries.csv')
deng_c = read_csv(file = '../data/processed/dengue_confirmed_timeseries.csv')
deng_s = read_csv(file = '../data/processed/dengue_suspected_timeseries.csv')


#=============================================================================#
# processing data
#=============================================================================#

# only keep places that are available for all 3 data sets
places_zik = data_frame(places = colnames(zik_s)[3:ncol(zik_c)])
places_chik = data_frame(places = colnames(chik_s)[3:ncol(chik_s)])
places_deng = data_frame(places = colnames(deng_s)[3:ncol(deng_s)])

places_common = inner_join(places_zik, places_chik) %>%
  inner_join(places_deng)
places_common = places_common$places


# subset datasets to only include those that have all the countries for all arboviruses
chik_c = cbind(chik_c[,1:2], chik_c[,which(colnames(chik_c) %in% places_common)])
chik_s = cbind(chik_s[,1:2], chik_s[,which(colnames(chik_s) %in% places_common)])
zik_c = cbind(zik_c[,1:2], zik_c[,which(colnames(zik_c) %in% places_common)])
zik_s = cbind(zik_s[,1:2], zik_s[,which(colnames(zik_s) %in% places_common)])
deng_c = cbind(deng_c[,1:2], deng_c[,which(colnames(deng_c) %in% places_common)])
deng_s = cbind(deng_s[,1:2], deng_s[,which(colnames(deng_s) %in% places_common)])


# need to add week 1 2014 - week 31 2015 to dat_zik... all 0's
zik_2014 = chik_s
zik_2014[, 3:ncol(zik_2014)] = 0
zik_2014 = zik_2014[1:(52 + 31), ]
zik_c = rbind(zik_2014, zik_c)
zik_s = rbind(zik_2014, zik_s)


# need to get rid of end of 2017 so dengue, chik, and zika all have same time frames
deng_c = deng_c[-which(deng_c$EW > 32 & deng_c$year == 2017),]
deng_s = deng_s[-which(deng_s$EW > 32 & deng_s$year == 2017), ]
chik_c = chik_c[-which(chik_c$Week > 32 & chik_c$Year == 2017), ]
chik_s = chik_s[-which(chik_s$Week > 32 & chik_s$Year == 2017), ]


# ordering dengue data the same way as the zika and chik data
deng_c_reordered = chik_c
deng_s_reordered = chik_s
for(dd in 3:ncol(chik_s)){
  deng_c_reordered[, dd] = deng_c[,which(colnames(deng_c) == colnames(chik_s)[dd])]
  deng_s_reordered[, dd] = deng_s[,which(colnames(deng_s) == colnames(chik_s)[dd])]
}
deng_c = deng_c_reordered
deng_s = deng_s_reordered


save(deng_c, deng_s, zik_c, zik_s, chik_s, chik_c, file = '../data/processed/processed_arbo_americas.RData')

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
countries = read_csv('../data/processed/country_codes.csv')

for(yy in 2014:2017){
  for(ww in 1:53){
    f = paste0('../data/raw/DENV_', yy, '_csv/W_By_Last_Available_EpiWeek_data-', ww, '.csv') 
    dat = read_csv(file = f) %>%
      select(c(`Country or Subregion`, EW, `Total of Dengue Cases (b)`, Confirmed)) %>%
      rename(place = `Country or Subregion`) %>%
      rename(total = `Total of Dengue Cases (b)`) %>%
      rename(confirmed = Confirmed)
    dat$place[which(str_detect(string = dat$place, pattern = "Cura\x8dao"))] = 'Curaçao'
    dat$year = rep(yy, nrow(dat))
    
    #=============================================================================#
    # joining ISO codes with dat
    #=============================================================================#
    
    dat_full = left_join(dat, countries, by = 'place')
    dat_full = dat_full[-which(is.na(dat_full$ISO)), ]
    
    if(yy == 2014 & ww == 1){
      dat_cum = dat_full 
    }else{
      dat_cum = bind_rows(dat_cum, dat_full)
    }
  }
}

#=============================================================================#
# processing data
#=============================================================================#

# if a country never reports data, drop that country
places_discard = c()
for(ii in 1:nrow(countries)){
  if(sum(dat_cum[which(dat_cum$ISO == countries$ISO[ii]), ]$total, 
         na.rm = T) == 0){
    places_discard = c(places_discard, countries$ISO[ii]) 
  }
}


dat_cum = dat_cum[-which(dat_cum$ISO == places_discard[1]),]
dat_cum = dat_cum[-which(dat_cum$ISO == places_discard[2]),]
countries = countries[-which(countries$ISO == places_discard[1]),]
countries = countries[-which(countries$ISO == places_discard[2]),]
countries = countries[-which(countries$ISO == places_discard[3]),]

# get rid of countries that report incomplete data
num_reports = data_frame(ISO = countries$ISO,
                         num_weeks = rep(NA, nrow(countries)))
for(ii in 1:nrow(num_reports)){
  num_reports$num_weeks[ii] = nrow(dat_cum[which(dat_cum$ISO==num_reports$ISO[ii]),])
}
places_discard = num_reports$ISO[which(num_reports$num_weeks != 53*4)]
dat_cum = dat_cum[-which(dat_cum$ISO==places_discard[1]),]
dat_cum = dat_cum[-which(dat_cum$ISO==places_discard[2]),]
dat_cum = dat_cum[-which(dat_cum$ISO==places_discard[3]),]
countries = countries[-which(countries$ISO == places_discard[1]),]
countries = countries[-which(countries$ISO == places_discard[2]),]
countries = countries[-which(countries$ISO == places_discard[3]),]

#=============================================================================#
# converting the data from cumulative incidence to weekly incidence
#=============================================================================#

# total cases 
dat_inc_c = select(dat_cum, -c(total, confirmed))
dat_inc_c$total_cases = rep(NA, nrow(dat_inc_c))
dat_inc_c$conf_cases = rep(NA, nrow(dat_inc_c))

dat_inc_s = select(dat_cum, -c(total, confirmed))
dat_inc_s$susp_cases = rep(NA, nrow(dat_inc_s))


for(cc in 1:nrow(countries)){
  for(yy in 2014:2017){
    for(ww in 1:52){
      ind1 = which(dat_inc_c$EW == ww 
                   & dat_inc_c$year == yy
                   & dat_inc_c$ISO == countries$ISO[cc])
      ind2 = which(dat_inc_c$EW == ww + 1 
                   & dat_inc_c$year == yy
                   & dat_inc_c$ISO == countries$ISO[cc])
      
      # weekly total cases
      dat_inc_c$total_cases[ind1] = 
        dat_cum$total[ind2] - dat_cum$total[ind1]
      
      if(is.na(dat_inc_c$total_cases[ind1])){
        dat_inc_c$total_cases[ind1] = 0
      }else if(dat_inc_c$total_cases[ind1] < 0){
        dat_inc_c$total_cases[ind1] = 0
      }
      
      # weekly confirmed cases
      dat_inc_c$conf_cases[ind1] = 
        dat_cum$confirmed[ind2] - dat_cum$confirmed[ind1]
      
      if(is.na(dat_inc_c$conf_cases[ind1])){
        dat_inc_c$conf_cases[ind1] = 0
      }else if(dat_inc_c$conf_cases[ind1] < 0){
        dat_inc_c$conf_cases[ind1] = 0
      }
      
      # weekly suspected cases
      dat_inc_s$susp_cases[ind1] = dat_inc_c$total_cases[ind1] - dat_inc_c$conf_cases[ind1]
      
      if(is.na(dat_inc_s$susp_cases[ind1])){
        dat_inc_s$susp_cases[ind1] = 0
      }else if(dat_inc_s$susp_cases[ind1] < 0){
        dat_inc_s$susp_cases[ind1] = 0
      }
    }
  }
}

dat_inc_c = dat_inc_c[-which(is.na(dat_inc_c$conf_cases)),]
dat_inc_s = dat_inc_s[-which(is.na(dat_inc_s$susp_cases)),]
# write_csv(dat_inc, '../data/processed/dengue_report_timeseries_long.csv')

# wide, timeseries data
dat_wide_c = matrix(NA, ncol = nrow(countries) + 2, nrow = (52*4))
dat_wide_c[, 1] = rep(2014:2017, each = 52)
dat_wide_c[, 2] = rep((1:52), 4)
colnames(dat_wide_c) = c('year', 'EW', countries$ISO)

dat_wide_s = dat_wide_c

for(ii in 1:nrow(countries)){
  # confirmed caes
  dat_wide_c[,which(colnames(dat_wide_c) == countries$ISO[ii])] = 
    dat_inc_c$conf_cases[which(dat_inc_c$ISO == countries$ISO[ii])]
  
  # suspected cases
  dat_wide_s[,which(colnames(dat_wide_s) == countries$ISO[ii])] = 
    dat_inc_s$susp_cases[which(dat_inc_s$ISO == countries$ISO[ii])]
}

write.csv(dat_wide_c, '../data/processed/dengue_confirmed_timeseries.csv', row.names = FALSE)
write.csv(dat_wide_s, '../data/processed/dengue_suspected_timeseries.csv', row.names = FALSE)

#=============================================================================#
# investigating the data
#=============================================================================#
# 
# # looking at data
# layout(matrix(1:48, 8,6))
# par(mar = c(1,1,0.5,0.5), oma = c(2,2,1,1))
# 
# for(ii in 1:nrow(countries)){
#   plot(dat_inc[which(dat_inc$ISO==countries$ISO[ii]),]$cases, type = 'l',
#        xaxt = 'n') 
#   text(10,1, countries$ISO[ii])
# }



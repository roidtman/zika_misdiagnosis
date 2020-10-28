# Processing PAHO DENV data
#=============================================================================#

#=============================================================================#
# setting up workspace
#=============================================================================#
rm(list = ls())
setwd('.')
require(tidyverse)
require(ISOcodes)
require(stringr)

#=============================================================================#
# loading in data
#=============================================================================#
countries_names = read_csv(file = '../data/raw/DENV_2014_csv/W_By_Last_Available_EpiWeek_data-1.csv') %>%
  select(`Country or Subregion`) %>%
  rename(place = `Country or Subregion`)

# find the ISO for each place
inds = rep(NA, nrow(countries_names))
for(ii in 1:nrow(countries_names)){
  if(! sum(str_detect(countries_names[ii,], ISO_3166_1$Name)) == 0){
    inds[ii] = which(str_detect(countries_names[ii,], ISO_3166_1$Name))
    # print(c(which(str_detect(countries_names[ii,], ISO_3166_1$Name)), ii))
  }
}

# create new dataset with ISO and place name
countries_combined = data_frame(place = c(ISO_3166_1$Name[inds]),
                                ISO = c(ISO_3166_1$Alpha_3[inds]))

# fill in place name where NA values located
for(ii in 1:length(which(is.na(inds)))){
  countries_combined$place[which(is.na(inds))[ii]] = 
    countries_names$place[which(is.na(inds))[ii]]
}

# manually adding ISO values
countries_combined$ISO[10] = 'BOL'
countries_combined$ISO[11] = 'BES'
countries_combined$ISO[20] = 'CUW'
countries_combined$place[22] = 'Dominican Republic'
countries_combined$ISO[22] = 'DOM'
countries_combined$place[20] = 'Curaçao'
countries_combined$ISO[47] = 'MAF'
countries_combined$ISO[44] = 'BLM'
countries_combined$ISO[56] = 'VEN'
countries_combined$ISO[57] = 'VGB'
countries_combined$ISO[58] = 'VIR'
countries_combined = countries_combined[-15,]

# get rid of places that don't have an ISO value
countries_combined = drop_na(countries_combined)

# produced csv file
write_csv(x = countries_combined, path = '../data/processed/country_codes.csv')


# Create ILO conditional forecasts from provided files
library(dplyr)

# this file contains a cubic panel dynamic dataset
all_cv    <- read.csv("inst/varia/ilo_cubic_panel.csv")
colnames(all_cv) = c("year", "iso3code", "country", "UR", "EPR", "LFPR", "dgdp")

# all variables all countries
data_cv <- all_cv %>% 
  filter(year >= 2024) %>% 
  select(year, iso3code, UR, EPR, LFPR, dgdp)

# Create a list with the country data
countries = unique(data_cv$iso3code)
countries = countries[order(countries)]
ilo_conditional_forecasts = list()
for (i in 1:length(countries)) {
  ilo_conditional_forecasts[[i]] <- data_cv %>% 
    filter(iso3code == countries[i]) %>% 
    select(UR, EPR, LFPR, dgdp) %>% 
    ts(start = 2024, frequency = 1)
  names(ilo_conditional_forecasts)[i] <- countries[i]
}

save(ilo_conditional_forecasts, file = "data/ilo_conditional_forecasts.rda")


# Create ILO dataset from provided files
library(dplyr)

# this file contains a dynamic panel dynamic dataset
all_cv    <- read.csv("inst/varia/ilo_dynamic_panel.csv")
colnames(all_cv) = c("year", "iso3code", "country", "GDP", "UR", "EPR", "LFPR")

# all variables all countries
data_cv <- all_cv %>% 
  filter(year <= 2023) %>% 
  mutate(gdp = log(GDP)) %>% 
  select(year, iso3code, gdp, UR, EPR, LFPR)

# Create a list with the country data
countries = unique(data_cv$iso3code)
countries = countries[order(countries)]
ilo_dynamic_panel = list()
for (i in 1:length(countries)) {
  ilo_dynamic_panel[[i]] <- data_cv %>% 
    filter(iso3code == countries[i]) %>% 
    select(gdp, UR, EPR, LFPR) %>% 
    ts(start = 1991, frequency = 1)
  names(ilo_dynamic_panel)[i] <- countries[i]
}

save(ilo_dynamic_panel, file = "data/ilo_dynamic_panel.rda")

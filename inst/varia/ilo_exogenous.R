
devtools::load_all()

# dummy panel variables
################################################
data(ilo_cubic_panel)
countries     = names(ilo_cubic_panel)

dummies       = matrix(0, nrow(ilo_cubic_panel$POL), 3)
colnames(dummies) = c("2008", "2020", "2021")
dummies       = ts(dummies, start = 1991, frequency = 1)
dummies[18,1] = dummies[30,2] = dummies[31,3] = 1

ilo_exogenous_variables             = list()
for (c in 1:length(ilo_cubic_panel)) {
  ilo_exogenous_variables[[c]]      = dummies
  names(ilo_exogenous_variables)[c] = countries[c]
}



# dummy panel forecasts
################################################
data(ilo_conditional_forecasts)

dummy_forecasts     = matrix(0, nrow(ilo_conditional_forecasts$POL), 3)
colnames(dummy_forecasts) = c("2024", "2025", "2026")
dummy_forecasts     = ts(dummy_forecasts, start = 2024, frequency = 1)

ilo_exogenous_forecasts             = list()
for (c in 1:length(ilo_conditional_forecasts)) {
  ilo_exogenous_forecasts[[c]]      = dummy_forecasts
  names(ilo_exogenous_forecasts)[c] = countries[c]
}


# save
################################################
save(ilo_exogenous_variables, file = "data/ilo_exogenous_variables.rda")
save(ilo_exogenous_forecasts, file = "data/ilo_exogenous_forecasts.rda")


data(ilo_dynamic_panel)

# for bsvar
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bvarPANEL$new(ilo_dynamic_panel)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2)

set.seed(1)
suppressMessages(
  ff2              <- ilo_dynamic_panel |>
    specify_bvarPANEL$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2)
)


expect_identical(
  ff[[1]]$forecasts[1,1,1], ff2[[1]]$forecasts[1,1,1],
  info = "forecast: forecast identical for normal and pipe workflow."
)

expect_true(
  is.numeric(ff[[1]]$forecasts) & is.array(ff[[1]]$forecasts),
  info = "forecast: returns numeric array."
)


expect_error(
  specify_bvarPANEL$new(ilo_dynamic_panel) |> forecast(horizon = 3),
  info = "forecast: wrong input provided."
)

expect_error(
  forecast(run_no1, horizon = 1.5),
  info = "forecast: specify horizon as integer."
)


# conditional forecasting
data(ilo_conditional_forecasts)

set.seed(1)
suppressMessages(
  specification_no1 <- specify_bvarPANEL$new(ilo_dynamic_panel)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, 6, conditional_forecast = ilo_conditional_forecasts)

set.seed(1)
suppressMessages(
  ff2              <- ilo_dynamic_panel |>
    specify_bvarPANEL$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 6, conditional_forecast = ilo_conditional_forecasts)
)


expect_identical(
  ff[[1]]$forecasts[1,1,1], ff2[[1]]$forecasts[1,1,1],
  info = "conditional forecast: forecast identical for normal and pipe workflow."
)

expect_equivalent(
  ilo_conditional_forecasts[[1]][1,1], ff2[[1]]$forecasts[1,1,1],
  info = "conditional forecast: forecasts and provided conditional forecasts identical."
)

expect_true(
  is.numeric(ff[[1]]$forecasts) & is.array(ff[[1]]$forecasts),
  info = "conditional forecast: returns numeric array."
)

expect_error(
  forecast(run_no1, horizon = 4, conditional_forecast = ilo_conditional_forecasts),
  pattern = "horizon",
  info = "conditional forecast: provided forecasts different from horizon."
)

expect_error(
  forecast(run_no1, horizon = 6, conditional_forecast = ilo_conditional_forecasts[-1]),
  info = "conditional forecast: uneven number of countries in forecasts and data."
)


# exogenous variables
data("ilo_exogenous_variables")
data("ilo_exogenous_forecasts")

set.seed(1)
suppressMessages(
  specification_no1 <- specify_bvarPANEL$new(ilo_dynamic_panel, exogenous = ilo_exogenous_variables)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, 6, exogenous_forecast = ilo_exogenous_forecasts)

set.seed(1)
suppressMessages(
  ff2              <- ilo_dynamic_panel |>
    specify_bvarPANEL$new(exogenous = ilo_exogenous_variables) |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 6, exogenous_forecast = ilo_exogenous_forecasts)
)


expect_identical(
  ff[[1]]$forecasts[1,1,1], ff2[[1]]$forecasts[1,1,1],
  info = "exogenous forecast: forecast identical for normal and pipe workflow."
)

expect_true(
  is.numeric(ff[[1]]$forecasts) & is.array(ff[[1]]$forecasts),
  info = "exogenous forecast: returns numeric array."
)


expect_error(
  forecast(run_no1, horizon = 4, exogenous_forecast = ilo_exogenous_forecasts),
  pattern = "horizon",
  info = "exogenous forecast: provided forecasts different from horizon."
)

expect_error(
  forecast(run_no1, horizon = 4),
  pattern = "missing",
  info = "exogenous forecast: missing exogenous forecasts for model with exogenous variables."
)

ilo_exogenous_forecasts[[1]][1,1] = NA
expect_error(
  forecast(run_no1, horizon = 6, exogenous_forecast = ilo_exogenous_forecasts),
  pattern = "values",
  info = "exogenous forecast: provided exogenous forecasts contain missing values."
)



# truncated forecasts

expect_identical(
  class(specify_bvarPANEL$new(ilo_dynamic_panel, type = c("real",rep("rate",3))))[1], 
  "BVARPANEL",
  info = "truncated forecast: good specification of argument type."
)

expect_error(
  specify_bvarPANEL$new(ilo_dynamic_panel, type = rep("rate",3)),
  pattern = "length",
  info = "truncated forecast: wrong specification of argument type."
)

set.seed(1)
suppressMessages(
  specification_no1 <- specify_bvarPANEL$new(ilo_dynamic_panel, type = c("real",rep("rate",3)))
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
suppressMessages(
  ff                  <- forecast(run_no1, 6)
)

expect_true(
  all(ff$POL$forecasts[2,,] >= 0),
  info = "truncated forecast: unemployment rates forecasts are non-negative."
)
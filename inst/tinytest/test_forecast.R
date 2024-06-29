
data(ilo_cubic_panel)

# for bsvar
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bvarPANEL$new(ilo_cubic_panel)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2)

set.seed(1)
suppressMessages(
  ff2              <- ilo_cubic_panel |>
    specify_bvarPANEL$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2)
)


expect_identical(
  ff$forecasts[1,1,1,1], ff2$forecasts[1,1,1,1],
  info = "forecast: forecast identical for normal and pipe workflow."
)

expect_true(
  is.numeric(ff$forecasts) & is.array(ff$forecasts),
  info = "forecast: returns numeric array."
)


expect_error(
  specify_bvarPANEL$new(ilo_cubic_panel) |> forecast(horizon = 3),
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
  specification_no1 <- specify_bvarPANEL$new(ilo_cubic_panel)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, 6, conditional_forecast = ilo_conditional_forecasts)

set.seed(1)
suppressMessages(
  ff2              <- ilo_cubic_panel |>
    specify_bvarPANEL$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 6, conditional_forecast = ilo_conditional_forecasts)
)


expect_identical(
  ff$forecasts[1,1,1,1], ff2$forecasts[1,1,1,1],
  info = "conditional forecast: forecast identical for normal and pipe workflow."
)

expect_true(
  is.numeric(ff$forecasts) & is.array(ff$forecasts),
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



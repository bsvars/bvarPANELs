data(ilo_dynamic_panel)

set.seed(1)
specification_no1   <- specify_bvarPANEL$new(ilo_dynamic_panel)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)

set.seed(1)
specification_no2   <- specify_bvarPANEL$new(ilo_dynamic_panel)
run_no2             <- estimate(specification_no2, 3, 1, show_progress = FALSE)

set.seed(1)
run_no3             <- ilo_dynamic_panel |>
  specify_bvarPANEL$new() |>
  estimate(S = 3, thin = 1, show_progress = FALSE)

expect_identical(
  run_no1$last_draw$starting_values$A[1,1],
  run_no2$last_draw$starting_values$A[1,1],
  info = "estimate_bvarPANEL: the last_draw(s) of two runs to be identical."
)

expect_identical(
  run_no1$posterior$A[1,1,1],
  run_no2$posterior$A[1,1,1],
  info = "estimate_bvarPANEL: the first draws of two runs to be identical."
)

expect_identical(
  run_no1$last_draw$starting_values$A[1,1],
  run_no3$last_draw$starting_values$A[1,1],
  info = "estimate_bvarPANEL: the last_draw(s) of a normal and pipe run to be identical."
)

# a test of a good setting of S and thin
expect_error(
  estimate(specification_no1, 3, 2, show_progress = FALSE),
  info = "Argument S is not a positive integer multiplication of argument thin."
)

expect_error(
  estimate(specification_no1, 2, 3, show_progress = FALSE),
  info = "Argument S is not a positive integer multiplication of argument thin."
)


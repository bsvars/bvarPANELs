
#' R6 Class Representing PriorBVARs
#'
#' @description
#' The class PriorBVARs presents a prior specification for the Bayesian
#' VAR model for each country.
#' 
#' @examples 
#' prior = specify_prior_bvars$new(C = 2, N = 3, p = 1)
#' prior$M
#' 
#' @export
specify_prior_bvars = R6::R6Class(
  "PriorBVARs",
  
  public = list(
    
    #' @field M an \code{KxN} matrix, the mean of the MNIW prior
    #' distribution for the autoregressive matrices \eqn{\mathbf{A}_c} 
    M           = matrix(),
    
    #' @field W a \code{KxK} column-specific covariance matrix of the
    #' MNIW prior distribution for the autoregressive matrices \eqn{\mathbf{A}_c}
    W           = matrix(),
    
    #' @field S_inv an \code{NxN} row-specific precision matrix of the
    #' MNIW prior distribution for the covariance matrices \eqn{\mathbf{\Sigma}_c}
    S_inv       = matrix(),
    
    #' @field lambda a positive shape of the exponential prior distribution 
    #' for the shape parameter \eqn{\nu}.
    lambda  = NA,
    
    #' @field mu_m a scalar mean of the normal prior distribution
    #' for the average persistence parameter \eqn{m}.
    mu_m   = NA,
    
    #' @field sigma2_m a positive scalar variance of the normal prior 
    #' distribution for the average persistence parameter \eqn{m}.
    sigma2_m  = NA,
    
    #' @field s_w a positive scalar scale of the inverse-gamma 2 prior 
    #' distribution for parameter \eqn{w}.
    s_w  = NA,
    
    #' @field nu_w a positive scalar shape of the inverse-gamma 2  prior 
    #' distribution for parameter \eqn{w}.
    nu_w  = NA,
    
    #' @field s_s a positive scalar scale parameter of the 
    #' gamma prior distribution for parameter \eqn{s}.
    s_s  = NA,
    
    #' @field a_s a positive scalar shape parameter of the
    #' gamma prior distribution for parameter \eqn{s}.
    a_s  = NA,
    
    #' @description
    #' Create a new prior specification PriorBVARs.
    #' @param C a positive integer - the number of countries in the data.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @param stationary an \code{N} logical vector - its element set to 
    #' \code{FALSE} sets the prior mean for the autoregressive parameters of the 
    #' \code{N}th equation to the random walk process, otherwise to white noise.
    #' @return A new prior specification PriorBVARs.
    #' @examples 
    #' # a prior for 2-country, 3-variable example with one lag and stationary data
    #' prior = specify_prior_bvars$new(C = 2, N = 3, p = 1)
    #' prior$M
    #' 
    initialize = function(C, N, p, d = 0, stationary = rep(FALSE, N)) {
      stopifnot("Argument C must be a positive integer number." = C > 0 & C %% 1 == 0)
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)
      stopifnot("Argument stationary must be a logical vector of length N." = length(stationary) == N & is.logical(stationary))
    
      K                 = N * p + 1 + d
      self$M            = t(cbind(diag(as.numeric(!stationary)), matrix(0, N, K - N)))
      self$W            = diag(c(kronecker((1:p)^2, rep(1, N) ), rep(10, 1 + d)))
      self$S_inv        = diag(N)
      self$lambda       = 30
      self$mu_m         = 1
      self$sigma2_m     = 0.1
      self$s_w          = 1
      self$nu_w         = 3
      self$s_s          = 1
      self$a_s          = 1
    }, # END initialize
    
    #' @description
    #' Returns the elements of the prior specification PriorBVARs as a \code{list}.
    #' 
    #' @examples 
    #' # a prior for 2-country, 3-variable example with four lags
    #' prior = specify_prior_bvars$new(C = 2, N = 3, p = 4)
    #' prior$get_prior() # show the prior as list
    #' 
    get_prior = function(){
      list(
        M        = self$M,
        W        = self$W,
        S_inv    = self$S_inv,
        lambda   = self$lambda,
        mu_m     = self$mu_m,
        sigma2_m = self$sigma2_m,
        s_w      = self$s_w,
        nu_w     = self$nu_w,
        s_s      = self$s_s,
        a_s      = self$a_s
      )
    } # END get_prior
  ) # END public
) # END specify_prior_bvars







#' R6 Class Representing StartingValuesBVARs
#'
#' @description
#' The class StartingValuesBVARs presents starting values for the Bayesian
#' hierarchical panel VAR model.
#' 
#' @examples 
#' # starting values for a Bayesian Panel VAR
#' sv = specify_starting_values_bvars$new(C = 2, N = 3, p = 1)
#' 
#' @export
specify_starting_values_bvars = R6::R6Class(
  "StartingValuesBVARs",
  
  public = list(
    
    #' @field A_c an \code{KxNxC} array of starting values for the local parameter 
    #' \eqn{\mathbf{A}_c}. 
    A_c           = array(),
    
    #' @field Sigma_c an \code{NxNxC} array of starting values for the local
    #' parameter \eqn{\mathbf{\Sigma}_c}. 
    Sigma_c       = array(),
    
    #' @field nu a \code{C}-vector of positive starting values for the parameter
    #' \eqn{\nu}.
    nu            = numeric(),
    
    #' @field m a \code{C}-vector of starting values for the parameter
    #' \eqn{m}.
    m             = numeric(),
    
    #' @field w a \code{C}-vector of positive starting values for the parameter
    #' \eqn{w}.
    w             = numeric(),
    
    #' @field s  a \code{C}-vector of positive starting values for the parameter
    #' \eqn{s}.
    s             = numeric(),
    
    #' @description
    #' Create new starting values StartingValuesBVARs
    #' @param C a positive integer - the number of countries in the data.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @return Starting values StartingValuesBVARs
    #' @examples 
    #' # starting values for Bayesian VARs 2-country model with 4 lags for a 3-variable system.
    #' sv = specify_starting_values_bvars$new(C = 2, N = 3, p = 4)
    #' 
    initialize = function(C, N, p, d = 0){
      stopifnot("Argument C must be a positive integer number." = C > 0 & C %% 1 == 0)
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)
      
      K               = N * p + 1 + d
      self$A_c        = array(stats::rnorm(C * K * N, sd = 0.001), c(K, N, C))
      self$Sigma_c    = stats::rWishart(C, N + 1, diag(N))
      self$nu         = rep(N + 1 + 0.1, C)
      self$m          = stats::rnorm(C, sd = 0.001)
      self$w          = stats::rgamma(C, 1)
      self$s          = stats::rgamma(C, 1)
    }, # END initialize
    
    #' @description
    #' Returns the elements of the starting values StartingValuesBVARs as 
    #' a \code{list}.
    #' 
    #' @examples 
    #' # starting values for bvars with 1 lag for a 3-variable system
    #' sv = specify_starting_values_bvars$new(C = 2, N = 3, p = 1)
    #' sv$get_starting_values()   # show starting values as list
    #' 
    get_starting_values   = function(){
      list(
        A_c           = self$A_c,
        Sigma_c       = self$Sigma_c,
        nu            = self$nu,
        m             = self$m,
        w             = self$w,
        s             = self$s
      )
    }, # END get_starting_values
    
    #' @description
    #' Returns the elements of the starting values StartingValuesBVARs as a \code{list}.
    #' @param last_draw a list containing the same elements as object StartingValuesBVARs.
    #' @return An object of class StartingValuesBVARs including the last draw 
    #' of the current MCMC as the starting value to be passed to the continuation 
    #' of the MCMC estimation.
    #' 
    #' @examples
    #' sv = specify_starting_values_bvars$new(C = 2, N = 3, p = 1)
    #' 
    #' # Modify the starting values by:
    #' sv_list = sv$get_starting_values()   # getting them as list
    #' sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
    #' sv$set_starting_values(sv_list)      # providing to the class object
    #' 
    set_starting_values   = function(last_draw) {
      self$A_c            = last_draw$A_c
      self$Sigma_c        = last_draw$Sigma_c
      self$nu             = last_draw$nu
      self$m              = last_draw$m
      self$w              = last_draw$w
      self$s              = last_draw$s
    } # END set_starting_values
  ) # END public
) # END specify_starting_values_bvars



#' R6 Class representing the specification of the BVARs model
#'
#' @description
#' The class BVARs presents complete specification for the Bayesian
#' Vector Autoregressions for cubic data.
#' 
#' @examples 
#' spec = specify_bvars$new(
#'    data = ilo_dynamic_panel
#' )
#' 
#' @export
specify_bvars = R6::R6Class(
  "BVARs",
  
  private = list(
    type  = "wozniak" # a value from set, \code{wozniak}, \code{jarocinski}, indicating model specification
  ),
  
  public = list(
    
    #' @field p a non-negative integer specifying the autoregressive lag order of the model. 
    p                      = numeric(),
    
    #' @field prior an object PriorBSVAR with the prior specification. 
    prior                  = list(),
    
    #' @field data_matrices an object DataMatricesBVARPANEL with the data matrices.
    data_matrices          = list(),
    
    #' @field starting_values an object StartingValuesBVARPANEL with the starting values.
    starting_values        = list(),
    
    #' @field adaptiveMH a vector of four values setting the adaptive MH sampler 
    #' for nu: adaptive rate, target acceptance rate, the iteration at which to 
    #' start adapting, the initial scaling rate
    adaptiveMH             = numeric(),
    
    #' @description
    #' Create a new specification of the Bayesian Panel VAR model BVARPANEL.
    #' @param data a list with \code{C} elements of \code{(T_c+p)xN} matrices 
    #' with time series data.
    #' @param p a positive integer providing model's autoregressive lag order.
    #' @param exogenous a \code{(T+p)xd} matrix of exogenous variables. 
    #' @param stationary an \code{N} logical vector - its element set to 
    #' \code{FALSE} sets the prior mean for the autoregressive parameters of the 
    #' \code{N}th equation to the random walk process, otherwise to white noise.
    #' @param type an \code{N} character vector with elements set to "rate" or "real"
    #' determining the truncation of the predictive density to \code{[0, 100]} and
    #' \code{(-Inf, Inf)} (no truncation) for each of the variables.
    #' @return A new complete specification for the Bayesian Panel VAR model BVARPANEL.
    initialize = function(
    data,
    p = 1L,
    exogenous = NULL,
    stationary = rep(FALSE, ncol(data[[1]])),
    type = rep("real", ncol(data[[1]]))
    ) {
      stopifnot("Argument data has to contain matrices with the same number of columns." 
                = length(unique(simplify2array(lapply(data, ncol)))) == 1)
      stopifnot("Argument p has to be a positive integer." 
                = ((p %% 1) == 0 & p > 0))
      
      self$p    = p
      C         = length(data)
      N         = unique(simplify2array(lapply(data, ncol)))
      d             = 0
      if (!is.null(exogenous)) {
        d           = ncol(exogenous[[1]])
      }
    
      self$data_matrices   = specify_panel_data_matrices$new(data, self$p, exogenous, type)
      self$prior           = specify_prior_bvars$new(C, N, self$p, d, stationary)
      self$starting_values = specify_starting_values_bvars$new(C, N, self$p, d)
      self$adaptiveMH      = c(0.44, 0.6)
      
    }, # END initialize
    
    #' @description
    #' Returns the data matrices as the DataMatricesBVARPANEL object.
    #' 
    #' @examples
    #' spec = specify_bvarPANEL$new(
    #'    data = ilo_dynamic_panel
    #' )
    #' spec$get_data_matrices()
    #' 
    get_data_matrices = function() {
      self$data_matrices$clone()
    }, # END get_data_matrices
    
    #' @description
    #' Returns the prior specification as the PriorBVARPANEL object.
    #' 
    #' @examples 
    #' spec = specify_bvarPANEL$new(
    #'    data = ilo_dynamic_panel
    #' )
    #' spec$get_prior()
    #' 
    get_prior = function() {
      self$prior$clone()
    }, # END get_prior
    
    #' @description
    #' Returns the starting values as the StartingValuesBVARPANEL object.
    #' 
    #' @examples 
    #' spec = specify_bvarPANEL$new(
    #'    data = ilo_dynamic_panel
    #' )
    #' spec$get_starting_values()
    #' 
    get_starting_values = function() {
      self$starting_values$clone()
    }, # END get_starting_values
    
    #' @description
    #' Returns the type of the model.
    #' 
    #' @examples 
    #' spec = specify_bvarPANEL$new(
    #'    data = ilo_dynamic_panel
    #' )
    #' spec$get_type()
    #' 
    get_type = function() {
      private$type
    }, # END get_starting_values
    
    
    #' @description
    #' Sets the VAR model priors to objective prior by Zellner (1972).
    #' 
    #' @references
    #' Zellner (1971). \emph{An Introduction to Bayesian Inference in Econometrics}. 
    #' John Wiley & Sons.
    #' 
    #' @examples
    #' spec = specify_bvars$new(
    #'    data = ilo_dynamic_panel
    #' )
    #' spec$set_prior2objective()
    #' 
    set_prior2objective = function() {
      
      stopifnot("This model cannot be adjusted. Only the default specification can." 
                = private$type == "wozniak")
      
      message("Setting the model priors to objective prior by Zellner (1972).")
      private$type                = "zellner"
      
      C = length(self$data_matrices$Y)
      N = ncol(self$data_matrices$Y[[1]])
      
      # values used in computations
      self$starting_values$nu   = rep(-(N + 1), C)
      self$starting_values$m    = rep(0, C)
      self$starting_values$w    = rep(0, C)
      self$starting_values$s    = rep(0, C)
      
    }, # END set_prior2objective
    
    
    #' @description
    #' Sets the prior mean of the global autoregressive parameters to the OLS 
    #' pooled panel estimator following Zellner, Hong (1989).
    #' 
    #' @param x a vector of four values setting the adaptive MH sampler for nu:
    #' adaptive rate, target acceptance rate, the iteration at which to 
    #' start adapting, the initial scaling rate
    #' 
    #' @references
    #' Zellner, Hong (1989). Forecasting international growth rates using 
    #' Bayesian shrinkage and other procedures. \emph{Journal of Econometrics}, 
    #' \bold{40}(1), 183–202, \doi{10.1016/0304-4076(89)90036-5}.
    #' 
    #' @examples 
    #' spec = specify_bvarPANEL$new(
    #'    data = ilo_dynamic_panel
    #' )
    #' spec$set_global2pooled()
    #' 
    set_global2pooled = function(x) {
      stopifnot("This model cannot be adjusted. Only the default specification can." 
                = private$type == "wozniak")
      
      C = length(self$data_matrices$Y)
      N = ncol(self$data_matrices$Y[[1]])
      d = ncol(self$data_matrices$exogenous[[1]])
      p = self$p
      K = N * p + d
      
      XX = matrix(0, K, K)
      XY = matrix(0, K, N)
      for (c in 1:C) {
        XcYc_tmp = .Call(`_bpvars_Y_c_and_X_c`, 
                         self$data_matrices$Y[[c]], 
                         self$data_matrices$exogenous[[c]],
                         p)
        
        Yc = XcYc_tmp[1,][[1]]
        Xc = XcYc_tmp[2,][[1]]
        
        not_missing = !apply(self$data_matrices$missing[[c]], 1, \(x)(any(x == 1)))
        if (sum(not_missing) != 0) {
          Yc = Yc[not_missing,]
          Xc = Xc[not_missing,]
          XX = XX + crossprod(Xc)
          XY = XY + crossprod(Xc, Yc)
        }
      }
      self$prior$M = solve(XX, XY)
    }, # END set_adaptiveMH
    
    #' @description
    #' Sets the parameters of adaptive Metropolis-Hastings sampler for the parameter nu.
    #' 
    #' @param x a vector of four values setting the adaptive MH sampler for nu:
    #' adaptive rate, target acceptance rate, the iteration at which to 
    #' start adapting, the initial scaling rate
    #' 
    #' @examples 
    #' spec = specify_bvarPANEL$new(
    #'    data = ilo_dynamic_panel
    #' )
    #' spec$set_adaptiveMH(c(0.6, 0.4, 10, 0.1))
    #' 
    set_adaptiveMH = function(x) {
      stopifnot("Argument x has to be a numeric vector of length 4." = length(x) == 4 & is.numeric(x))
      stopifnot("Argument x must contain positive values." = all(x > 0))
      stopifnot("The second element of argument x must be less than 1." = x[2] < 1)
      stopifnot("The third element of argument x must greater or equal to 1." = x[3] >= 1)
      x[3]            = floor(x[3])
      self$adaptiveMH = x
    } # END set_adaptiveMH
  ) # END public
) # END specify_bvarPANEL




#' R6 Class Representing PosteriorBVARs
#'
#' @description
#' The class PosteriorBVARs contains posterior output and the specification 
#' including the last MCMC draw for the Bayesian Panel VAR model. 
#' Note that due to the thinning of the MCMC output the starting value in element 
#' \code{last_draw} might not be equal to the last draw provided in 
#' element \code{posterior}.
#' 
#' @seealso \code{\link{specify_bvars}}
#' 
#' @examples 
#' specification = specify_bvars$new(
#'    data = ilo_dynamic_panel[1:5]
#' )
#' posterior       = estimate(specification, 5)
#' class(posterior)
#' 
#' @export
specify_posterior_bvars = R6::R6Class(
  "PosteriorBVARs",
  
  private = list(
    normalised = FALSE
  ), # END private
  
  public = list(
    
    #' @field last_draw an object of class BVARs with the last draw of the 
    #' current MCMC run as the starting value to be passed to the continuation 
    #' of the MCMC estimation using \code{estimate()}. 
    last_draw = list(),
    
    #' @field posterior a list containing Bayesian estimation output.
    posterior = list(),
    
    #' @description
    #' Create a new posterior output PosteriorBVARs.
    #' @param specification_bvarPANEL an object of class BVARs with the last 
    #' draw of the current MCMC run as the starting value.
    #' @param posterior_bvarPANEL a list containing Bayesian estimation output.
    #' @return A posterior output PosteriorBVARs.
    initialize = function(specification_bvarPANEL, posterior_bvarPANEL) {
      
      stopifnot("Argument specification_bvarPANEL must be of class BVARs." = any(class(specification_bvarPANEL) == "BVARs"))
      stopifnot("Argument posterior_bvarPANEL must must contain MCMC output." = is.list(posterior_bvarPANEL))
      
      self$last_draw    = specification_bvarPANEL
      self$posterior    = posterior_bvarPANEL
      
      N = dim(specification_bvarPANEL$starting_values$A_c)[2]
      K = dim(specification_bvarPANEL$starting_values$A_c)[1]
      C = dim(specification_bvarPANEL$starting_values$A_c)[3]
      S = dim(posterior_bvarPANEL$nu)[2]
      
      Sigma_c           = array(NA, c(N, N, C, S))
      A_c               = array(NA, c(K, N, C, S))
      for (s in 1:S) {
        A_c[,,,s]       = posterior_bvarPANEL$A_c_cpp[s,1][[1]]
        Sigma_c[,,,s]   = posterior_bvarPANEL$Sigma_c_cpp[s,1][[1]]
      }
      self$posterior$Sigma_c   = Sigma_c
      self$posterior$A_c       = A_c
      
      
    }, # END initialize
    
    #' @description
    #' Returns a list containing Bayesian estimation output.
    #' 
    #' @examples 
    #' specification = specify_bvars$new(
    #'    data = ilo_dynamic_panel[1:5]
    #' )
    #' posterior       = estimate(specification, 5)
    #' posterior$get_posterior()
    #' 
    get_posterior       = function(){
      self$posterior
    }, # END get_posterior
    
    #' @description
    #' Returns an object of class BVARs with the last draw of the current 
    #' MCMC run as the starting value to be passed to the continuation of the 
    #' MCMC estimation using \code{estimate()}.
    #' 
    #' @examples
    #' specification = specify_bvars$new(
    #'    data = ilo_dynamic_panel[1:5]
    #' )
    #' burn_in        = estimate(specification, 5)
    #' posterior      = estimate(burn_in, 5)
    #' 
    get_last_draw      = function(){
      self$last_draw$clone()
    } # END get_last_draw
    
  ) # END public
) # END specify_posterior_bvars



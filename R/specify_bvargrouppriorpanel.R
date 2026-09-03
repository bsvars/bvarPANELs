

#' R6 Class Representing StartingValuesBVARGROUPPRIORPANEL
#'
#' @description
#' The class StartingValuesBVARGROUPPRIORPANEL presents starting values for the Bayesian
#' hierarchical panel VAR model with country grouping
#' 
#' @examples 
#' # starting values for a Bayesian Panel VAR
#' sv = specify_starting_values_bvarGroupPriorPANEL$new(rep(1,2), C = 2, G = 1, N = 3, p = 1)
#' 
#' @export
specify_starting_values_bvarGroupPriorPANEL = R6::R6Class(
  "StartingValuesBVARGROUPPRIORPANEL",
  
  public = list(
    
    #' @field group_allocation a numeric vector with integer numbers denoting group allocations
    group_allocation = numeric(),
    
    #' @field A_c an \code{KxNxC} array of starting values for the local parameter 
    #' \eqn{\mathbf{A}_c}. 
    A_c           = array(),
    
    #' @field Sigma_c an \code{NxNxC} array of starting values for the local
    #' parameter \eqn{\mathbf{\Sigma}_c}. 
    Sigma_c       = array(),
    
    #' @field A_g an \code{KxNxG} array of starting values for the group parameter 
    #' \eqn{\mathbf{A}_g}. 
    A_g           = array(),
    
    #' @field Sigma_g an \code{NxNxG} array of starting values for the group
    #' parameter \eqn{\mathbf{\Sigma}_g}. 
    Sigma_g       = array(),
    
    #' @field V an \code{KxK} matrix of starting values for the global parameter 
    #' \eqn{\mathbf{V}}. 
    V             = matrix(),
    
    #' @field nu a positive scalar with starting values for the global parameter
    #' \eqn{\nu}.
    nu            = NA,
    
    #' @field m a positive scalar with starting values for the global hyper-parameter
    #' \eqn{m}.
    m             = NA,
    
    #' @field w a positive scalar with starting values for the global hyper-parameter
    #' \eqn{w}.
    w             = NA,
    
    #' @field s a positive scalar with starting values for the global hyper-parameter
    #' \eqn{s}.
    s             = NA,
    
    
    #' @description
    #' Create new starting values StartingValuesBVARGROUPPRIORPANEL
    #' @param group_allocation a numeric vector with integer numbers denoting group allocations
    #' @param C a positive integer - the number of countries in the data.
    #' @param G a positive integer specifying the number of country groups.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @return Starting values StartingValuesBVARGROUPPRIORPANEL
    #' @examples 
    #' # starting values for Bayesian Panel VAR 2-country model with 4 lags for a 3-variable system.
    #' sv = specify_starting_values_bvarGroupPriorPANEL$new(C = 2, N = 3, p = 1)
    #' 
    initialize = function(group_allocation = 1:C, C, G = C, N, p, d = 0){
      
      stopifnot("Argument C must be a positive integer number." = C > 0 & C %% 1 == 0)
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)
      stopifnot("Argument group_allocation must be a numeric vector with integer numbers denoting groups." = 
                  is.numeric(group_allocation) & all(group_allocation %% 1 == 0))
      stopifnot("Argument group_allocation must be of length C." = 
                  (length(group_allocation) == C))
      
      K               = N * p + 1 + d
      self$A_c        = array(stats::rnorm(C * K * N, sd = 0.001), c(K, N, C))
      self$Sigma_c    = stats::rWishart(C, N + 1, diag(N))
      self$V          = stats::rWishart(1, K + 1, diag(K))[,,1]
      self$nu         = N + 1 + 0.1
      self$m          = stats::rnorm(1, sd = 0.001)
      self$w          = stats::rgamma(1, 1)
      self$s          = stats::rgamma(1, 1)
      
      self$group_allocation = group_allocation
      self$A_g        = array(stats::rnorm(G * K * N, sd = 0.001), c(K, N, G))
      self$Sigma_g    = stats::rWishart(G, N + 1, diag(N))
      
    }, # END initialize
    
    
    #' @description
    #' Returns the elements of the starting values StartingValuesBVARGROUPPRIORPANEL as 
    #' a \code{list}.
    #' 
    #' @examples 
    #' # starting values for a homoskedastic bsvar with 1 lag for a 3-variable system
    #' sv = specify_starting_values_bvarGroupPriorPANEL$new(rep(1,2), C = 2, N = 3, p = 1)
    #' sv$get_starting_values()   # show starting values as list
    #' 
    get_starting_values   = function(){
      list(
        group_allocation = self$group_allocation,
        A_c           = self$A_c,
        Sigma_c       = self$Sigma_c,
        A_g           = self$A_g,
        Sigma_g       = self$Sigma_g,
        V             = self$V,
        nu            = self$nu,
        m             = self$m,
        w             = self$w,
        s             = self$s
      )
    }, # END get_starting_values
    
    #' @description
    #' Returns the elements of the starting values StartingValuesBVARGROUPPRIORPANEL as a \code{list}.
    #' @param last_draw a list containing the same elements as object StartingValuesBVARGROUPPRIORPANEL
    #' @return An object of class StartingValuesBVARGROUPPRIORPANEL including the last draw 
    #' of the current MCMC as the starting value to be passed to the continuation 
    #' of the MCMC estimation.
    #' 
    #' @examples
    #' sv = specify_starting_values_bvarGroupPriorPANEL$new(rep(1,2), C = 2, N = 3, p = 1)
    #' 
    #' # Modify the starting values by:
    #' sv_list = sv$get_starting_values()   # getting them as list
    #' sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
    #' sv$set_starting_values(sv_list)      # providing to the class object
    #' 
    set_starting_values   = function(last_draw) {
      self$group_allocation = last_draw$group_allocation
      self$A_c            = last_draw$A_c
      self$Sigma_c        = last_draw$Sigma_c
      self$A_g            = last_draw$A_g
      self$Sigma_g        = last_draw$Sigma_g
      self$V              = last_draw$V
      self$nu             = last_draw$nu
      self$m              = last_draw$m
      self$w              = last_draw$w
      self$s              = last_draw$s
    } # END set_starting_values
  ) # END public
)  # END specify_starting_values_bvarGroupPriorPANEL
  


#' R6 Class representing the specification of the BVARGROUPPRIORPANEL model
#'
#' @description
#' The class BVARGROUPPRIORPANEL presents complete specification for the Bayesian Panel
#' Vector Autoregression with county grouping for global prior parameters. The groups can be pre-specified, 
#' which requires the argument \code{group_allocation} to be provided, or estimated,
#' which requires the argument \code{G} for the number of groups to be provided 
#' and the argument \code{group_allocation} to be left empty.
#' 
#' @examples 
#' spec = specify_bvarGroupPriorPANEL$new(
#'    data = ilo_dynamic_panel,
#'    G = 2
#' )
#' 
#' @export
specify_bvarGroupPriorPANEL = R6::R6Class(
  "BVARGROUPPRIORPANEL",
  
  private = list(
    type  = "wozniak" # a value from set, \code{wozniak}, \code{jarocinski}, indicating model specification
  ),
  
  public = list(
    
    #' @field p a non-negative integer specifying the autoregressive lag order of the model. 
    p                      = numeric(),
    
    #' @field G a non-negative integer specifying the number of country groupings. 
    G                      = numeric(),
    
    #' @field estimate_groups a logical value denoting whether the groups are to be estimated. 
    estimate_groups        = logical(),
    
    #' @field prior an object PriorBSVAR with the prior specification. 
    prior                  = list(),
    
    #' @field data_matrices an object DataMatricesBVARPANEL with the data matrices.
    data_matrices          = list(),
    
    #' @field starting_values an object StartingValuesBVARGROUPPRIORPANEL with the starting values.
    starting_values        = list(),
    
    #' @field adaptiveMH a vector of four values setting the adaptive MH sampler 
    #' for nu: adaptive rate, target acceptance rate, the iteration at which to 
    #' start adapting, the initial scaling rate
    adaptiveMH             = numeric(),
    
    #' @description
    #' Create a new specification of the Bayesian Panel VAR model with country 
    #' grouping for global prior parameters BVARGROUPPRIORPANEL. The groups can be pre-specified, 
    #' which requires the argument \code{group_allocation} to be provided, or estimated,
    #' which requires the argument \code{G} for the number of groups to be provided 
    #' and the argument \code{group_allocation} to be left empty.
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
    #' @param G a positive integer specifying the number of country groups. Its 
    #' specification is required if \code{group_allocation} is not provided and 
    #' the country groups to be estimated.
    #' @param group_allocation an argument that can be provided as a numeric 
    #' vector with integer numbers denoting group allocations to pre-specify the
    #' the country groups, in which case they are not estimated, or left empty
    #' if the country groups are to be estimated.
    #' @return A new complete specification for the Bayesian Panel VAR model BVARPANEL.
    initialize = function(
      data,
      p = 1L,
      exogenous = NULL,
      stationary = rep(FALSE, ncol(data[[1]])),
      type = rep("real", ncol(data[[1]])),
      G = NULL,
      group_allocation = NULL
    ) {
      stopifnot("Argument data has to contain matrices with the same number of columns." 
                = length(unique(simplify2array(lapply(data, ncol)))) == 1)
      stopifnot("Argument p has to be a positive integer." 
                = ((p %% 1) == 0 & p > 0))
      
      if (!is.null(G)) {
        
        stopifnot("Argument G has to be a positive integer." = ((G %% 1) == 0 & G > 0))
        
        if (!is.null(group_allocation)) {
          stopifnot("Number of groups in argument group_allocation must be equal to the value of argument G." 
                    = length(unique(group_allocation)) == G)
          message("Country groupings have been pre-specified and will not be estimated.")
          self$estimate_groups = FALSE
        } else {
          message("Country groupings will be estimated.")
          group_allocation = rep(1:G, length.out = length(data))
          self$estimate_groups = TRUE
        }
        
      } else {
        
        if (!is.null(group_allocation)) {
          G = length(unique(group_allocation))
          message("Country groupings have been pre-specified and will not be estimated.")
          self$estimate_groups = FALSE
        } else {
          stop("Either the number of groups G or the group_allocation has to be provided.")
        }
        
      }
      
      self$p    = p
      self$G    = G
      C         = length(data)
      N         = unique(simplify2array(lapply(data, ncol)))
      d             = 0
      if (!is.null(exogenous)) {
        d           = ncol(exogenous[[1]])
      }
      
      self$data_matrices   = specify_panel_data_matrices$new(data, self$p, exogenous, type)
      self$prior           = specify_prior_bvarPANEL$new(C, N, self$p, d, stationary)
      self$starting_values = specify_starting_values_bvarGroupPriorPANEL$new(group_allocation, C, self$G, N, self$p, d)
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
      stopifnot("Argument x has to be a numeric vector of length 4." 
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
    }, # END set_global2pooled
    
    #' @description
    #' Sets the parameters of adaptive Metropolis-Hastings sampler for the parameter nu.
    #' 
    #' @param x a vector of four values setting the adaptive MH sampler for nu:
    #' adaptive rate, target acceptance rate, the iteration at which to 
    #' start adapting, the initial scaling rate
    #' 
    #' @examples 
    #' spec = specify_bvarPANEL$new(
    #'    data = ilo_dynamic_panel[1:5]
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
) # END specify_bvarGroupPriorPANEL



#' R6 Class Representing PosteriorBVARGROUPPRIORPANEL
#'
#' @description
#' The class PosteriorBVARGROUPPRIORPANEL contains posterior output and the specification 
#' including the last MCMC draw for the Bayesian Panel VAR model with country grouping. 
#' Note that due to the thinning of the MCMC output the starting value in element 
#' \code{last_draw} might not be equal to the last draw provided in 
#' element \code{posterior}.
#' 
#' @seealso \code{\link{specify_bvarGroupPriorPANEL}}
#' 
#' @examples 
#' spec = specify_bvarGroupPriorPANEL$new(
#'    data = ilo_dynamic_panel[1:5],
#'    G = 2
#' )
#' #posterior       = estimate(specification, 5)
#' #class(posterior)
#' 
#' @export
specify_posterior_bvarGroupPriorPANEL = R6::R6Class(
  "PosteriorBVARGROUPPRIORPANEL",
  
  inherit = specify_posterior_bvarGroupPANEL,
  
  public = list(
    
    #' @field last_draw an object of class BVARGROUPPRIORPANEL with the last draw of the 
    #' current MCMC run as the starting value to be passed to the continuation 
    #' of the MCMC estimation using \code{estimate()}. 
    last_draw = list(),
    
    #' @field posterior a list containing Bayesian estimation output.
    posterior = list(),
    
    #' @description
    #' Create a new posterior output PosteriorBVARGROUPPRIORPANEL.
    #' @param specification an object of class BVARGROUPPRIORPANEL with the last 
    #' draw of the current MCMC run as the starting value.
    #' @param posterior a list containing Bayesian estimation output.
    #' @return A posterior output PosteriorBVARGROUPPRIORPANEL.
    initialize = function(specification, posterior) {
      
      stopifnot("Argument specification must be of class BVARGROUPPRIORPANEL." = any(class(specification) == "BVARGROUPPRIORPANEL"))
      stopifnot("Argument posterior must must contain MCMC output." = is.list(posterior) & is.array(posterior$V))
      
      self$last_draw    = specification
      self$posterior    = posterior
      
      N = dim(specification$starting_values$A_c)[2]
      K = dim(specification$starting_values$A_c)[1]
      C = dim(specification$starting_values$A_c)[3]
      G = dim(specification$starting_values$A_g)[3]
      S = dim(posterior$V)[3]
      
      Sigma_c           = array(NA, c(N, N, C, S))
      A_c               = array(NA, c(K, N, C, S))
      Sigma_g           = array(NA, c(N, N, G, S))
      A_g               = array(NA, c(K, N, G, S))
      
      for (s in 1:S) {
        A_c[,,,s]       = posterior$A_c_cpp[s,1][[1]]
        Sigma_c[,,,s]   = posterior$Sigma_c_cpp[s,1][[1]]
        A_g[,,,s]       = posterior$A_g_cpp[s,1][[1]]
        Sigma_g[,,,s]   = posterior$Sigma_g_cpp[s,1][[1]]
      }
      self$posterior$Sigma_c   = Sigma_c
      self$posterior$A_c       = A_c
      self$posterior$Sigma_g   = Sigma_g
      self$posterior$A_g       = A_g

    } # END initialize
  ) # END public
) # END specify_posterior_bvarGroupPriorPANEL

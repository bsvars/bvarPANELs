
#' R6 Class Representing PriorBVARPANEL
#'
#' @description
#' The class PriorBVARPANEL presents a prior specification for the Bayesian
#' hierarchical panel VAR model.
#' 
#' @examples 
#' prior = specify_prior_bvarPANEL$new(C = 2, N = 3, p = 1)
#' prior$M
#' 
#' @export
specify_prior_bvarPANEL = R6::R6Class(
  "PriorBVARPANEL",
  
  public = list(
    
    #' @field M an \code{KxN} matrix, the mean of the second-level MNIW prior
    #' distribution for the global parameter matrices \eqn{\undelline{\mathbf{A}}} 
    #' and \eqn{\undelline{\mathbf{V}}}
    M           = matrix(),
    
    #' @field W a \code{KxK} column-specific covariance matrix of the second-level
    #' MNIW prior distribution for the global parameter matrices \eqn{\undelline{\mathbf{A}}}
    #' and \eqn{\undelline{\mathbf{V}}}
    W           = matrix(),
    
    #' @field S_inv an \code{NxN} row-specific precision matrix of the second-level
    #' MNIW prior distribution for the global parameter matrices \eqn{\undelline{\mathbf{A}}}
    #' and \eqn{\undelline{\mathbf{V}}}
    S_inv       = matrix(),
    
    #' @field S_Sigma_inv an \code{NxN} precision matrix of the second-level 
    #' Wishart prior distribution for the global parameter matrix \eqn{\undelline{\mathbf{\Sigma}}}.
    S_Sigma_inv = matrix(),
    
    #' @field eta a positive shape parameter of the second-level MNIW prior distribution
    #' for the global parameter matrices \eqn{\undelline{\mathbf{A}}}
    #' and \eqn{\undelline{\mathbf{V}}}
    eta = NA,
    
    #' @field mu a positive shape parameter of the second-level Wishart prior distribution
    #'  for the global parameter matrix \eqn{\undelline{\mathbf{\Sigma}}}.
    mu  = NA,
    
    #' @field lambda a positive shape of the second-level exp prior distribution 
    #' for the shape parameter \eqn{\undelline{\nu}}.
    lambda  = NA,
    
    #' @field mu_m a scalar mean of the third-level normal prior distribution
    #' for the global average persistence parameter \eqn{\undelline{m}}.
    mu_m   = NA,
    
    #' @field sigma2_m a positive scalar variance of the third-level normal prior distribution
    #' for the global average persistence parameter \eqn{\undelline{m}}.
    sigma2_m  = NA,
    
    #' @field s_w a positive scalar scale of the third-level gamma prior 
    #' distribution for parameter \eqn{\underline{w}}.
    s_w  = NA,
    
    #' @field a_w a positive scalar shape of the third-level gamma prior 
    #' distribution for parameter \eqn{\underline{w}}.
    a_w  = NA,
    
    #' @field s_s a positive scalar scale parameter of the third-level 
    #' inverted-gamma 2 prior distribution for parameter \eqn{\underline{s}}.
    s_s  = NA,
    
    #' @field nu_s a positive scalar shape parameter of the third-level 
    #' inverted-gamma 2 prior distribution for parameter \eqn{\underline{s}}.
    nu_s  = NA,
    
    #' @description
    #' Create a new prior specification PriorBVARPANEL.
    #' @param C a positive integer - the number of countries in the data.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @return A new prior specification PriorBVARPANEL.
    #' @examples 
    #' # a prior for 2-country, 3-variable example with one lag and stationary data
    #' prior = specify_prior_bvarPANEL$new(C = 2, N = 3, p = 1)
    #' prior$M
    #' 
    initialize = function(C, N, p, d = 0){
      stopifnot("Argument C must be a positive integer number." = C > 0 & C %% 1 == 0)
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)

      K                 = N * p + 1 + d
      self$M            = cbind(diag(N), matrix(0, N, K - N))
      self$W            = diag(c(kronecker((1:p)^2, rep(1, N) ), rep(10, d + 1)))
      self$S_inv        = diag(N)
      self$S_Sigma_inv  = diag(N)
      self$eta          = N + 1
      self$mu           = N + 1
      self$lambda       = 0.1
      self$mu_m         = 1
      self$sigma2_m     = 1
      self$s_w          = 1
      self$a_w          = 1
      self$s_s          = 1
      self$nu_s         = 3
    }, # END initialize
    
    #' @description
    #' Returns the elements of the prior specification PriorBSVAR as a \code{list}.
    #' 
    #' @examples 
    #' # a prior for 2-coutnry, 3-variable example with four lags
    #' prior = specify_prior_bvarPANEL$new(C = 2, N = 3, p = 4)
    #' prior$get_prior() # show the prior as list
    #' 
    get_prior = function(){
      list(
        M        = self$M,
        W        = self$W,
        S_inv    = self$S_inv,
        S_Sigma_inv = self$S_Sigma_inv,
        eta      = self$eta,
        mu       = self$mu,
        lambda   = self$lambda,
        mu_m     = self$mu_m,
        sigma2_m = self$sigma2_m,
        s_w      = self$s_w,
        a_w      = self$a_w,
        s_s      = self$s_s,
        nu_s     = self$nu_s
      )
    } # END get_prior
    
  ) # END public
) # END specify_prior_bvarPANEL

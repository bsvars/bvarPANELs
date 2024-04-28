
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
    #' distribution for the global parameter matrices \eqn{\mathbf{A}} 
    #' and \eqn{\mathbf{V}}
    M           = matrix(),
    
    #' @field W a \code{KxK} column-specific covariance matrix of the second-level
    #' MNIW prior distribution for the global parameter matrices \eqn{\mathbf{A}}
    #' and \eqn{\mathbf{V}}
    W           = matrix(),
    
    #' @field S_inv an \code{NxN} row-specific precision matrix of the second-level
    #' MNIW prior distribution for the global parameter matrices \eqn{\mathbf{A}}
    #' and \eqn{\mathbf{V}}
    S_inv       = matrix(),
    
    #' @field S_Sigma_inv an \code{NxN} precision matrix of the second-level 
    #' Wishart prior distribution for the global parameter matrix \eqn{\mathbf{\Sigma}}.
    S_Sigma_inv = matrix(),
    
    #' @field eta a positive shape parameter of the second-level MNIW prior distribution
    #' for the global parameter matrices \eqn{\mathbf{A}}
    #' and \eqn{\mathbf{V}}
    eta = NA,
    
    #' @field mu_Sigma a positive shape parameter of the second-level Wishart prior 
    #' distribution  for the global parameter matrix \eqn{\mathbf{\Sigma}}.
    mu_Sigma  = NA,
    
    #' @field lambda a positive shape of the second-level exp prior distribution 
    #' for the shape parameter \eqn{\nu}.
    lambda  = NA,
    
    #' @field mu_m a scalar mean of the third-level normal prior distribution
    #' for the global average persistence parameter \eqn{m}.
    mu_m   = NA,
    
    #' @field sigma2_m a positive scalar variance of the third-level normal prior 
    #' distribution for the global average persistence parameter \eqn{m}.
    sigma2_m  = NA,
    
    #' @field s_w a positive scalar scale of the third-level gamma prior 
    #' distribution for parameter \eqn{w}.
    s_w  = NA,
    
    #' @field a_w a positive scalar shape of the third-level gamma prior 
    #' distribution for parameter \eqn{w}.
    a_w  = NA,
    
    #' @field s_s a positive scalar scale parameter of the third-level 
    #' inverted-gamma 2 prior distribution for parameter \eqn{s}.
    s_s  = NA,
    
    #' @field nu_s a positive scalar shape parameter of the third-level 
    #' inverted-gamma 2 prior distribution for parameter \eqn{s}.
    nu_s  = NA,
    
    #' @description
    #' Create a new prior specification PriorBVARPANEL.
    #' @param C a positive integer - the number of countries in the data.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @return A new prior specification PriorBVARPANEL.
    #' @examples 
    #' # a prior for 2-country, 3-variable example with one lag and stationary data
    #' prior = specify_prior_bvarPANEL$new(C = 2, N = 3, p = 1)
    #' prior$M
    #' 
    initialize = function(C, N, p){
      stopifnot("Argument C must be a positive integer number." = C > 0 & C %% 1 == 0)
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
    
      K                 = N * p + 1
      self$M            = cbind(diag(N), matrix(0, N, K - N))
      self$W            = diag(c(kronecker((1:p)^2, rep(1, N) ), rep(10, 1)))
      self$S_inv        = diag(N)
      self$S_Sigma_inv  = diag(N)
      self$eta          = N + 1
      self$mu_Sigma     = N + 1
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
        mu_Sigma = self$mu_Sigma,
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







#' R6 Class Representing StartingValuesBVARPANEL
#'
#' @description
#' The class StartingValuesBVARPANEL presents starting values for the Bayesian
#' hierarchical panel VAR model.
#' 
#' @examples 
#' # starting values for a Bayesian Panel VAR
#' sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)
#' 
#' @export
specify_starting_values_bvarPANEL = R6::R6Class(
  "StartingValuesBVARPANEL",
  
  public = list(
    
    #' @field A_c an \code{KxNxC} array of starting values for the local parameter 
    #' \eqn{\mathbf{A}_c}. 
    A_c           = array(),
    
    #' @field Sigma_c an \code{NxNxC} array of starting values for the local
    #' parameter \eqn{\mathbf{\Sigma}_c}. 
    Sigma_c       = array(),
    
    #' @field A an \code{KxN} matrix of starting values for the global parameter 
    #' \eqn{\mathbf{A}}. 
    A             = matrix(),
    
    #' @field V an \code{KxK} matrix of starting values for the global parameter 
    #' \eqn{\mathbf{V}}. 
    V             = matrix(),
    
    #' @field Sigma an \code{NxN} matrix of starting values for the global parameter 
    #' \eqn{\mathbf{\Sigma}}. 
    Sigma         = matrix(),
    
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
    #' Create new starting values StartingValuesBVARPANEL
    #' @param C a positive integer - the number of countries in the data.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @return Starting values StartingValuesBVARPANEL
    #' @examples 
    #' # starting values for Bayesian Panel VAR 2-country model with 4 lags for a 3-variable system.
    #' sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 4)
    #' 
    initialize = function(C, N, p){
      stopifnot("Argument C must be a positive integer number." = C > 0 & C %% 1 == 0)
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      
      K               = N * p + 1
      self$A_c        = matrix(stats::rnorm(C * K * N, sd = 0.001), c(K, N, C))
      self$Sigma_c    = stats::rWishart(C, N + 1, diag(N))
      self$A          = matrix(stats::rnorm(K * N, sd = 0.001), K, N) + diag(K)[,1:N]
      self$V          = stats::rWishart(1, K + 1, diag(K))[,,1]
      self$Sigma      = stats::rWishart(1, N + 1, diag(N))[,,1]
      self$nu         = stats::rgamma(1, 3)
      self$m          = stats::rnorm(1, sd = 0.001)
      self$w          = stats::rgamma(1, 1)
      self$s          = stats::rgamma(1, 1)
    }, # END initialize
    
    #' @description
    #' Returns the elements of the starting values StartingValuesBVARPANEL as 
    #' a \code{list}.
    #' 
    #' @examples 
    #' # starting values for a homoskedastic bsvar with 1 lag for a 3-variable system
    #' sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)
    #' sv$get_starting_values()   # show starting values as list
    #' 
    get_starting_values   = function(){
      list(
        A_c           = self$A_c,
        Sigma_c       = self$Sigma_c,
        A             = self$A,
        V             = self$V,
        Sigma         = self$Sigma,
        nu            = self$nu,
        m             = self$m,
        w             = self$w,
        s             = self$s
      )
    }, # END get_starting_values
    
    #' @description
    #' Returns the elements of the starting values StartingValuesBVARPANEL as a \code{list}.
    #' @param last_draw a list containing the same elements as object StartingValuesBVARPANEL.
    #' @return An object of class StartingValuesBVARPANEL including the last draw 
    #' of the current MCMC as the starting value to be passed to the continuation 
    #' of the MCMC estimation.
    #' 
    #' @examples
    #' sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)
    #' 
    #' # Modify the starting values by:
    #' sv_list = sv$get_starting_values()   # getting them as list
    #' sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
    #' sv$set_starting_values(sv_list)      # providing to the class object
    #' 
    set_starting_values   = function(last_draw) {
      self$A_c            = last_draw$A_c
      self$Sigma_c        = last_draw$Sigma_c
      self$A              = last_draw$A
      self$V              = last_draw$V
      self$Sigma          = last_draw$Sigma
      self$nu             = last_draw$nu
      self$m              = last_draw$m
      self$w              = last_draw$w
      self$s              = last_draw$s
    } # END set_starting_values
  ) # END public
) # END specify_starting_values_bvarPANEL






#' R6 Class Representing DataMatricesBVARPANEL
#'
#' @description
#' The class DataMatricesBVARPANEL presents the data matrices of dependent 
#' variables, \eqn{\mathbf{Y}_c}, and regressors, \eqn{\mathbf{X}_c}, for the 
#' Bayesian Panel VAR model for all countries \eqn{c = 1, ..., C}.
#' 
#' @export
specify_panel_data_matrices = R6::R6Class(
  "DataMatricesBVARPANEL",
  
  public = list(
    
    #' @field Y a list with \code{C} elements with \code{T_c x N} matrices of 
    #' dependent variables, \eqn{\mathbf{Y}_c}. 
    Y     = list(),
    
    #' @field X a list with \code{C} elements with \code{T_c x K} matrices of 
    #' regressors, \eqn{\mathbf{X}_c}. 
    X     = list(),
    
    #' @description
    #' Create new data matrices DataMatricesBVARPANEL
    #' @param data a list containing \code{(T_c+p)xN} matrices with country-specific
    #' time series data.
    #' @param p a positive integer providing model's autoregressive lag order.
    #' @return New data matrices DataMatricesBVARPANEL
    initialize = function(data, p = 1L) {
      if (missing(data)) {
        stop("Argument data has to be specified")
      } else {
        stopifnot("Argument data has to be a list of matrices." = is.list(data) & all(simplify2array(lapply(data, function(x){is.matrix(x) & is.numeric(x)}))))
        stopifnot("Argument data has to contain matrices with the same number of columns." = length(unique(simplify2array(lapply(data, ncol)))) == 1)
        stopifnot("Argument data cannot include missing values." = all(simplify2array(lapply(dd, function(x){!any(is.na(x))}))))
      }
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      
      C             = length(data)
      for (c in 1:C) {
        TT            = nrow(data[[c]])
        T_c           = TT - p
        self$Y[[c]]   = data[[c]][(p + 1):TT,]
        
        X             = matrix(0, T_c, 0)
        for (i in 1:p) {
          X           = cbind(X, data[[c]][(p + 1):TT - i,])
        }
        X             = cbind(X, rep(1, T_c))
        self$X        = X
      } # END c loop
      names(self$Y)   = names(self$X) = names(data)
    }, # END initialize
    
    #' @description
    #' Returns the data matrices DataMatricesBVARPANEL as a \code{list}.
    get_data_matrices = function() {
      list(
        Y = self$Y,
        X = self$X
      )
    } # END get_data_matrices
  ) # END public
) # END specify_panel_data_matrices







#' R6 Class representing the specification of the BVARPANEL model
#'
#' @description
#' The class BVARPANEL presents complete specification for the Bayesian Panel
#' Vector Autoregression.
#' 
#' @examples 
#' \dontrun{
#' data(us_fiscal_lsuw)
#' spec = specify_bvarPANEL$new(
#'    data = us_fiscal_lsuw,
#'    p = 4
#' )
#' }
#' 
#' @export
specify_bvarPANEL = R6::R6Class(
  "BVARPANEL",
  
  public = list(
    
    #' @field p a non-negative integer specifying the autoregressive lag order of the model. 
    p                      = numeric(),
    
    #' @field prior an object PriorBSVAR with the prior specification. 
    prior                  = list(),
    
    #' @field data_matrices an object DataMatricesBVARPANEL with the data matrices.
    data_matrices          = list(),
    
    #' @field starting_values an object StartingValuesBVARPANEL with the starting values.
    starting_values        = list(),
    
    #' @description
    #' Create a new specification of the Bayesian Panel VAR model BVARPANEL.
    #' @param data a list with \code{C} elements of \code{(T_c+p)xN} matrices 
    #' with time series data.
    #' @param p a positive integer providing model's autoregressive lag order.
    #' @return A new complete specification for the Bayesian Panel VAR model BVARPANEL.
    initialize = function(
    data,
    p = 1L
    ) {
      stopifnot("Argument data has to contain matrices with the same number of columns." = length(unique(simplify2array(lapply(data, ncol)))) == 1)
      stopifnot("Argument p has to be a positive integer." = ((p %% 1) == 0 & p > 0))
      
      self$p    = p
      C         = length(data)
      N         = unique(simplify2array(lapply(data, ncol)))
      
      self$data_matrices   = specify_panel_data_matrices$new(data, self$p)
      self$prior           = specify_prior_bvarPANEL$new(C, N, self$p)
      self$starting_values = specify_starting_values_bvarPANEL$new(C, N, self$p)
    }, # END initialize
    
    #' @description
    #' Returns the data matrices as the DataMatricesBVARPANEL object.
    #' 
    #' @examples
    #' \dontrun{ 
    #' data(us_fiscal_lsuw)
    #' spec = specify_bsvar$new(
    #'    data = us_fiscal_lsuw,
    #'    p = 4
    #' )
    #' spec$get_data_matrices()
    #' }
    get_data_matrices = function() {
      self$data_matrices$clone()
    }, # END get_data_matrices
    
    #' @description
    #' Returns the prior specification as the PriorBVARPANEL object.
    #' 
    #' @examples 
    #' \dontrun{
    #' data(us_fiscal_lsuw)
    #' spec = specify_bsvar$new(
    #'    data = us_fiscal_lsuw,
    #'    p = 4
    #' )
    #' spec$get_prior()
    #' }
    get_prior = function() {
      self$prior$clone()
    }, # END get_prior
    
    #' @description
    #' Returns the starting values as the StartingValuesBVARPANEL object.
    #' 
    #' @examples 
    #' \dontrun{
    #' data(us_fiscal_lsuw)
    #' spec = specify_bsvar$new(
    #'    data = us_fiscal_lsuw,
    #'    p = 4
    #' )
    #' spec$get_starting_values()
    #' }
    get_starting_values = function() {
      self$starting_values$clone()
    } # END get_starting_values
  ) # END public
) # END specify_bvarPANEL


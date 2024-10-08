% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/specify_bvarpanel.R
\name{specify_starting_values_bvarPANEL}
\alias{specify_starting_values_bvarPANEL}
\title{R6 Class Representing StartingValuesBVARPANEL}
\description{
The class StartingValuesBVARPANEL presents starting values for the Bayesian
hierarchical panel VAR model.
}
\examples{
# starting values for a Bayesian Panel VAR
sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)


## ------------------------------------------------
## Method `specify_starting_values_bvarPANEL$new`
## ------------------------------------------------

# starting values for Bayesian Panel VAR 2-country model with 4 lags for a 3-variable system.
sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 4)


## ------------------------------------------------
## Method `specify_starting_values_bvarPANEL$get_starting_values`
## ------------------------------------------------

# starting values for a homoskedastic bsvar with 1 lag for a 3-variable system
sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)
sv$get_starting_values()   # show starting values as list


## ------------------------------------------------
## Method `specify_starting_values_bvarPANEL$set_starting_values`
## ------------------------------------------------

sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)

# Modify the starting values by:
sv_list = sv$get_starting_values()   # getting them as list
sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
sv$set_starting_values(sv_list)      # providing to the class object

}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{A_c}}{an \code{KxNxC} array of starting values for the local parameter 
\eqn{\mathbf{A}_c}.}

\item{\code{Sigma_c}}{an \code{NxNxC} array of starting values for the local
parameter \eqn{\mathbf{\Sigma}_c}.}

\item{\code{A}}{an \code{KxN} matrix of starting values for the global parameter 
\eqn{\mathbf{A}}.}

\item{\code{V}}{an \code{KxK} matrix of starting values for the global parameter 
\eqn{\mathbf{V}}.}

\item{\code{Sigma}}{an \code{NxN} matrix of starting values for the global parameter 
\eqn{\mathbf{\Sigma}}.}

\item{\code{nu}}{a positive scalar with starting values for the global parameter
\eqn{\nu}.}

\item{\code{m}}{a positive scalar with starting values for the global hyper-parameter
\eqn{m}.}

\item{\code{w}}{a positive scalar with starting values for the global hyper-parameter
\eqn{w}.}

\item{\code{s}}{a positive scalar with starting values for the global hyper-parameter
\eqn{s}.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-StartingValuesBVARPANEL-new}{\code{specify_starting_values_bvarPANEL$new()}}
\item \href{#method-StartingValuesBVARPANEL-get_starting_values}{\code{specify_starting_values_bvarPANEL$get_starting_values()}}
\item \href{#method-StartingValuesBVARPANEL-set_starting_values}{\code{specify_starting_values_bvarPANEL$set_starting_values()}}
\item \href{#method-StartingValuesBVARPANEL-clone}{\code{specify_starting_values_bvarPANEL$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-StartingValuesBVARPANEL-new"></a>}}
\if{latex}{\out{\hypertarget{method-StartingValuesBVARPANEL-new}{}}}
\subsection{Method \code{new()}}{
Create new starting values StartingValuesBVARPANEL
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specify_starting_values_bvarPANEL$new(C, N, p, d = 0)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{C}}{a positive integer - the number of countries in the data.}

\item{\code{N}}{a positive integer - the number of dependent variables in the model.}

\item{\code{p}}{a positive integer - the autoregressive lag order of the SVAR model.}

\item{\code{d}}{a positive integer - the number of \code{exogenous} variables in the model.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Starting values StartingValuesBVARPANEL
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# starting values for Bayesian Panel VAR 2-country model with 4 lags for a 3-variable system.
sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 4)

}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-StartingValuesBVARPANEL-get_starting_values"></a>}}
\if{latex}{\out{\hypertarget{method-StartingValuesBVARPANEL-get_starting_values}{}}}
\subsection{Method \code{get_starting_values()}}{
Returns the elements of the starting values StartingValuesBVARPANEL as 
a \code{list}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specify_starting_values_bvarPANEL$get_starting_values()}\if{html}{\out{</div>}}
}

\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# starting values for a homoskedastic bsvar with 1 lag for a 3-variable system
sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)
sv$get_starting_values()   # show starting values as list

}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-StartingValuesBVARPANEL-set_starting_values"></a>}}
\if{latex}{\out{\hypertarget{method-StartingValuesBVARPANEL-set_starting_values}{}}}
\subsection{Method \code{set_starting_values()}}{
Returns the elements of the starting values StartingValuesBVARPANEL as a \code{list}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specify_starting_values_bvarPANEL$set_starting_values(last_draw)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{last_draw}}{a list containing the same elements as object StartingValuesBVARPANEL.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
An object of class StartingValuesBVARPANEL including the last draw 
of the current MCMC as the starting value to be passed to the continuation 
of the MCMC estimation.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)

# Modify the starting values by:
sv_list = sv$get_starting_values()   # getting them as list
sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
sv$set_starting_values(sv_list)      # providing to the class object

}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-StartingValuesBVARPANEL-clone"></a>}}
\if{latex}{\out{\hypertarget{method-StartingValuesBVARPANEL-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specify_starting_values_bvarPANEL$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

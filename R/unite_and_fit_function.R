#' Combine Fitted and Non-fitted Parameters
#'
#' Combines fitted and non-fitted parameters and calls the fitting function. Serves as a wrapping function for the user-specified fitting function \code{fit.fn} (see \code{\link{create.profile}}).
#' @param par A named vector containing all parameters that are supposed to be fitted.
#' @param no.fit A named parameter vector containing all parameters that are not supposed to be fitted.
#' @param par.names The names of all parameters
#' @param fit.fn The cost function (see \code{\link{create.profile}} for more details).
#' @param ... Other arguments.
#' @export
#' @return Returns the negative log-likelihood as calculated by the specified cost function
#' @examples
#' #set parameters and cost function
#' fit.par <- c(p1 = 2, p2 = 4)
#' name.par <- c("p1", "p2", "p3")
#' defaults <- list(p1 = 0, p2 = 2, p3 = 4)
#' cost.function <- function(parms){
#'     parms[1] + parms[2] + parms[3]
#' }
#'
#' #call unite.and.fit
#' unite.and.fit(par = fit.par, no.fit = c(p3 = 1),par.names = name.par, fit.fn = cost.function)
unite.and.fit <- function(par, no.fit, par.names, fit.fn, ...) {

  dots <- list(...)
  # combine fitted and non-fitted parameters
  total.par  <- rep(0, length(par.names))
  names(total.par) <- par.names
  total.par[names(par)] <- par
  total.par[names(no.fit)] <- no.fit


  # call user defined function while passing on user defined variables
  diff <- do.call(fit.fn, c(list(parms = total.par), dots))

  return(as.numeric(diff))
}

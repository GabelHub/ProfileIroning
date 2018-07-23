#' Create Profile Likelihood Data Frame
#'
#' Given previously calculated profile values by \code{\link{create.profile}}, \code{get.profile} merges the results together into a data frame.
#' @param which.par The names of the parameter that the profile should be returned for.
#' @param range A numeric vector determining the values for which the profile should be calculated.
#' @param homedir The directory in which the folder \code{"Profile-Results"} is found. Default to \code{\link{getwd}}().
#' @param wait Logical. If true, \code{get.profile} will wait if a certain profile value has not been fitted yet (i.e. the corresponding status file is not equal to "done").
#'
#' @return A data frame containing the negative log-likelihood and the fitted parameter values for a specific range created by \code{\link{create.profile}}.
#' @export
#'
#' @examples
#' #define cost function
#' cost_function <- function(parms){
#' y <- parms[1] + parms[2]*c(1:3)
#' LL <- sum((y - c(1:3))^2)
#' }
#'
#' #create profile values
#' create.profile(which.par = "p1",
#'                par.names = c(p1 = 1, p2 = 3),
#'                range = list(seq(0,5,1)),
#'                fit.fn = cost_function,
#'                future.off = TRUE)
#'
#' #retrieve the calculated profile
#' res <- get.profile(which.par = "p1", range = seq(0,5,1))

get.profile <- function(which.par, range, homedir = getwd(), wait = FALSE){
  res <- NULL
  for(i in 1:length(range)){
    if(wait == FALSE){
      res <- rbind(res, readRDS(paste0(homedir, "/Profile-Results/Fits/", which.par,"_", range[i], ".rds")))
    }else{
      while(readRDS(paste0(homedir, "/Profile-Results/Status/status", which.par,"_", range[i], ".rds")) != "done" || file.exists(paste0(homedir, "/Profile-Results/Fits/", which.par,"_", range[i], ".rds")) == FALSE){
        Sys.sleep(5)
      }
      res <- rbind(res, readRDS(paste0(homedir, "/Profile-Results/Fits/", which.par,"_", range[i], ".rds")))
    }

  }
  return(res)
}

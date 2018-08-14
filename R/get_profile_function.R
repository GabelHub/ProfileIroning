#' Create Profile Likelihood Data Frame
#'
#' Given previously calculated profile values by \code{\link{create.profile}}, \code{get.profile} merges the results together into a data frame.
#' @param which.par The names of the parameter that the profile should be returned for.
#' @param range A numeric vector determining the values for which the profile should be calculated. Alternatively, setting it to "get.all" retrieves all available profile values and binds them together.
#' @param homedir The directory in which the folder \code{"Profile-Results"} is found. Default to \code{\link{getwd}}().
#' @param wait Logical. If true, \code{get.profile} will wait if a certain profile value has not been fitted yet (i.e. the corresponding status file is not equal to "done").
#' @param delete.old Logical. If TRUE, the individual point-wise fits created by \code{\link{point.profile}} will be deleted after using them. Default to FALSE.
#' @param save.it Logical. If TRUE (default), the gathered data will be bound and saved in the "Profile-Results/Tables/" folder.
#'
#' @return A data frame containing the negative log-likelihood and the fitted parameter values for a specific range created by \code{\link{create.profile}}.
#' @export
#'
#' @examples
#' #define cost function
#' cost_function <- function(parms){
#' y <- parms[1] + parms[2]*c(1:3) + parms[3]*sin(c(1:3))
#' LL <- sum((y - c(1:3))^2)
#' }
#'
#' #create profile values
#' create.profile(which.par = "get.p1",
#'                par.names = c(get.p1 = 1, get.p2 = 3, get.p3 = -2),
#'                range = list(seq(0,5,1)),
#'                fit.fn = cost_function,
#'                future.off = TRUE)
#'
#' #retrieve the calculated profile
#' res <- get.profile(which.par = "get.p1", range = seq(0,5,1))

get.profile <- function(which.par, range, homedir = getwd(), wait = FALSE, delete.old = FALSE, save.it = TRUE){

  if(range[1] == "get.all"){
    range <- c()
    all.files <- list.files(paste0(homedir, "/Profile-Results/Fits"), pattern = which.par)
    for(i in 1:length(all.files)){
      val <- gsub(pattern = paste0(which.par, "_"),replacement = "", all.files[i])
      val <- gsub(pattern = paste0(".rds"),replacement = "", val)
      range <- c(range, as.numeric(val))
      range <- unique(range[order(range)])
    }
  }

  res <- c()

  for(i in 1:length(range)){
    if(wait == FALSE){
      res <- rbind(res, readRDS(paste0(homedir, "/Profile-Results/Fits/", which.par,"_", range[i], ".rds")))
    }else{
      while(readRDS(paste0(homedir, "/Profile-Results/Status/status", which.par,"_", range[i], ".rds")) != "done" || file.exists(paste0(homedir, "/Profile-Results/Fits/", which.par,"_", range[i], ".rds")) == FALSE){
        print(paste0("Waiting for ", which.par,"_", range[i], "..."))
        Sys.sleep(10)
      }
      res <- rbind(res, readRDS(paste0(homedir, "/Profile-Results/Fits/", which.par,"_", range[i], ".rds")))
      if(delete.old == TRUE){
        file.remove(paste0(homedir, "/Profile-Results/Fits/", which.par,"_", range[i], ".rds"))
      }

    }

  }

  res <- as.data.frame(res)
  if(save.it == TRUE){
    saveRDS(res, paste0(homedir, "/Profile-Results/Tables/", which.par, ".rds"))
  }

  return(res)
}

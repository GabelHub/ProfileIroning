#' Calculating the Profile Likelihood
#'
#' @param which.par Name of the parameter that a profile should be obtained for. If multiple profiles are supposed to be calculated, a vector can be passed along as well. Alternatively, supplying "all.par" calculates the profiles for all available parameters.
#' @param par.names A vector containing the names and initial values for all available parameters.
#' @param range A list containing the respective ranges for which the profile should be calculated.
#' @param fit.fn A cost function. Has to take the complete parameter vector as an input (needs to be names \code{parms}) and must return the corresponding negative log-likelihood (-2LL, see Burnham and Anderson 2002).
#' @param bind.old Logical. If TRUE, previously calculated values will also be added to the profile if available. Default to FALSE.
#' @param delete.old Logical. If TRUE, the individual point-wise fits created by \code{\link{point.profile}} will be deleted after using them. Default to FALSE.
#' @param do.not.fit A named vector containing the values of the parameters that should not be fitted. Default to NULL.
#' @param homedir The directory to which the results should be saved to.
#' @param optim.runs The number of times that each model will be fitted by \code{\link{optim}}. Default to 5.
#' @param random.borders The ranges from which the random initial parameter conditions for all \code{optim.runs} larger than one are sampled. Can be either given as a vector containing the relative deviations for all parameters or as a matrix containing in its first column the lower and in its second column the upper border values. Parameters are uniformly sampled based on \code{\link{runif}}. Default to 1 (100\% deviation of all parameters). Alternatively, functions such as \code{\link{rnorm}}, \code{\link{rchisq}}, etc. can be used if the additional arguments are passed along as well.
#' @param refit If TRUE, previously fitted ranges will be fitted again and results will be overwritten according to the value set in \code{save.rel.diff}. Default to FALSE. Works only if \code{delete.old} is set to FALSE.
#' @param save.rel.diff A numeric value indicating a relative threshold when to overwrite a pre-existing result. Default to 0, which means that results get overwritten if an improvement is made.
#' @param con.tol The absolute convergence tolerance of each fitting run (see Details). Default is set to 0.1.
#' @param control.optim Control parameters passed along to \code{optim}. For more details, see \code{\link{optim}}.
#' @param parscale.pars Logical. If TRUE (default), the \code{parscale} option will be used when fitting with \code{\link{optim}}. This is helpful, if the parameter values are on different scales.
#' @param future.off Logical. If TRUE, \code{\link{future}} will not be used to calculate the results. Default to FALSE.
#' @param ... Additional parameters that can be passed along to \code{\link{future}} or \code{fit.fn}.
#' @return A list containing the respective profile values for every specified parameter.
#' @export
#'
#' @examples
#' #create data with standard deviation of 1
#' x.values <- 1:7
#' y.values <-  9 * x.values^2 - exp(2 * x.values)
#' sd.y.values <- rep(1,7)
#'
#' #define initial parameter values
#' inits <- c(p1 = 3, p2 = 4, p3 = -2, p4 = 2, p5 = 0)
#'
#' #define cost function that returns the negative log-likelihood
#' cost_function <- function(parms, x.vals, y.vals, sd.y){
#'   # restrict the search range to -5 to +5
#'   if(max(abs(parms)) > 5){
#'     return(NA)
#'   }
#'   with(as.list(c(parms)), {
#'     res <- p1*4 + p2*x.vals + p3^2*x.vals^2 + p4*sin(x.vals)  - exp(p5*x.vals)
#'     diff <- sum((res - y.vals)^2/sd.y)
#'   })
#' }
#'
#' #perform model selection
#' res <- create.profile(which.par = "all.par",
#'                       par.names = inits,
#'                       range = list(seq(0, 2, 0.1),
#'                                  seq(0, 5, 1),
#'                                  seq(2.9, 3.1, 0.01),
#'                                  seq(0, 3, 0.1),
#'                                  seq(1.999999, 2.000001, 0.0000001)),
#'                       fit.fn = cost_function,
#'                       optim.runs = 1,
#'                       delete.old = TRUE,
#'                       x.vals = x.values,
#'                       y.vals = y.values,
#'                       sd.y = sd.y.values)
create.profile <- function(which.par, par.names, range, fit.fn, bind.old = FALSE, delete.old = FALSE, do.not.fit = NULL, homedir = getwd(), optim.runs = 5, random.borders = 1, refit = FALSE, save.rel.diff = 0, con.tol = 0.1, control.optim = list(maxit = 1000), parscale.pars = TRUE, future.off = FALSE, ...){
  #create storage directories
  create.directories(homedir = homedir)

  #check if refitting makes sense
  if(refit == TRUE && delete.old == TRUE){
    stop("Set 'delete.old = FALSE' to allow for refitting.")
  }

  #get all parameters
  if(which.par[1] == "all.par"){
    index <- 1:length(par.names)
    #check if non-fitted parameter is fitted
    if(length(do.not.fit) > 0){
      stop("One of the profile parameters is specified in the do.not.fit vector!")
    }
  }else{
    index <- which(names(par.names) %in% which.par)
    #check if non-fitted parameter is fitted
    if(sum(is.element(names(do.not.fit), which.par) > 0)){
      stop("One of the profile parameters is specified in the do.not.fit vector!")
    }
  }
  overall.min <- Inf
  future.list <- list()
  for(i in 1:length(index)){
    print(paste0("Submitting fits for parameter ",names(par.names)[index[i]]))
    range.x <- range[[i]]
    for(j in 1:length(range.x)){
      if(refit == TRUE || file.exists(paste0(homedir, "/Profile-Results/Fits/",paste0(names(par.names)[index[i]], "_", range.x[j]), ".rds")) == FALSE){

        nofit <- range.x[j]
        names(nofit) <- names(par.names)[index[i]]
        nofit <- c(nofit, do.not.fit)

        if(future.off == TRUE){
          point.profile(no.fit = nofit,
                        parms = par.names,
                        fit.fn = fit.fn,
                        homedir = homedir,
                        optim.runs = optim.runs,
                        random.borders = random.borders,
                        con.tol = con.tol,
                        control.optim = control.optim,
                        parscale.pars = parscale.pars,
                        save.rel.diff = save.rel.diff,
                        ...)
        }else{
          future.list <- c(future.list, list(future::future(point.profile(no.fit = nofit,
                                                                          parms = par.names,
                                                                          fit.fn = fit.fn,
                                                                          homedir = homedir,
                                                                          optim.runs = optim.runs,
                                                                          random.borders = random.borders,
                                                                          con.tol = con.tol,
                                                                          control.optim = control.optim,
                                                                          parscale.pars = parscale.pars,
                                                                          save.rel.diff = save.rel.diff,
                                                                          ...),
                                                            label = paste0(names(par.names)[index[i]], "_", range.x[j]),
                                                            ...)))
        }
      }
    }
  }


  if(length(future.list) > 0){
    #check if futures are done
    for(wff in 1:length(future.list)){
      print.wait <- TRUE
      while(!future::resolved(future.list[[wff]])){
        if(print.wait){
          print(paste0("Waiting for future ", future.list[[wff]]$label, " ..."))
          print.wait <- FALSE
        }
        Sys.sleep(5)
      }
    }

    #check if all files exist
    for(fex in 1:length(future.list)){
      if(!file.exists(paste0(homedir, "/Profile-Results/Fits/", future.list[[wff]]$label, ".rds"))){
        #print error message
        if(class(try(future::value(future.list[[wff]]))) == "try-error"){
          try(future::value(future.list[[wff]]))
        }

        stop(paste0("Future is done but no output file was generated. The corresponding error message of job ",
                    future.list[[wff]]$label,
                    " is shown above. If no output is shown, use 'future.off = TRUE' to debug."))
      }
    }
  }

  #plot the outcome

  nrows = round(sqrt(length(index)))
  ncols = ifelse(length(index) == 1, 1, ceiling(sqrt(length(index))))
  all.res <- list()

  for(i in 1:length(index)){
    if(bind.old == TRUE){
      extra.range <- c()
      all.files <- list.files(paste0(homedir, "/Profile-Results/Fits"), pattern = names(par.names)[index[i]])
      if(length(all.files) >= 1){
        for(j in 1:length(all.files)){
          val <- gsub(pattern = paste0(names(par.names)[index[i]], "_"), replacement = "", all.files[j])
          val <- gsub(pattern = paste0(".rds"), replacement = "", val)
          extra.range <- c(extra.range, as.numeric(val))
        }
      }

      get.range <- unique(c(extra.range, as.vector(range[[i]])))
      get.range <- get.range[order(get.range)]

    }else{
      get.range <- range[[i]]
    }
    table.x <- get.profile(which.par = names(par.names)[index[i]],
                           range = get.range,
                           homedir = homedir,
                           delete.old = delete.old)
    if(i == 1){
      overall.min <- min(table.x[,1])
    }else{
      if(overall.min > min(table.x[,1])){
        overall.min <- min(table.x[,1])
      }
    }

    if(length(index) == 1){
      all.res <- table.x
    }else{
      all.res[[i]] <- table.x
    }

  }
  if(length(index) > 1){
    names(all.res) <- names(par.names)[index]
  }

  old.par <- graphics::par("mfrow")
  on.exit(graphics::par(mfrow = old.par))

  graphics::par(mfrow = c(nrows, ncols))

  for(i in 1:length(index)){
    if(length(index) > 1){
      table.x <- all.res[[i]]
    }

    graphics::plot(table.x[, names(par.names)[index[i]]],
                   table.x$LL,
                   type = "l",
                   lwd = 3,
                   xlab = names(par.names)[index[i]],
                   ylab = "-2LL",
                   main = paste0("Profile likelihood of ", names(par.names)[index[i]]))

    graphics::abline(h = min(table.x$LL) + 3.84, lwd = 3, lty = 2, col = "red")
    graphics::abline(h = overall.min, lty = 3, col = "grey")

  }

  #save plots
  grDevices::pdf(file = paste0(homedir, "/Profile-Results/Figures/RawProfiles.pdf"),
                 width  = 4*ncols,
                 height = 4*nrows,
                 useDingbats = F)

  graphics::par(mfrow = c(nrows, ncols))

  for(i in 1:length(index)){
    if(length(index) > 1){
      table.x <- all.res[[i]]
    }

    graphics::plot(table.x[, names(par.names)[index[i]]],
                   table.x$LL,
                   type = "l",
                   lwd = 3,
                   xlab = names(par.names)[index[i]],
                   ylab = "-2LL",
                   main = paste0("Profile likelihood of ", names(par.names)[index[i]]))

    graphics::abline(h = min(table.x$LL) + 3.84, lwd = 3, lty = 2, col = "red")
    graphics::abline(h = overall.min, lty = 3, col = "grey")

  }

  grDevices::dev.off()


  return(all.res)
}

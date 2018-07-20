#' Fit a Single Profile Point
#'
#' For a fixed value of one of the parameters, \code{point.profile} fits the remaining parameters and stores the results in the folder "Profile-Results/Fits/" to be accessed by \code{\link{create.profile}} later.
#' @param no.fit A named vector containing the values of all parameters that are not to be fitted.
#' @param parms A named vector containing the values of all parameters.
#' @param fit.fn A cost function. Has to take the complete parameter vector as an input (needs to be names \code{parms}) and must return the corresponding negative log-likelihood (-2LL, see Burnham and Anderson 2002).
#' @param homedir The directory to which the results should be saved to.
#' @param optim.runs The number of times that each model will be fitted by \code{\link{optim}}. Default to 5.
#' @param random.borders The ranges from which the random initial parameter conditions for all \code{optim.runs} larger than one are sampled. Can be either given as a vector containing the relative deviations for all parameters or as a matrix containing in its first column the lower and in its second column the upper border values. Parameters are uniformly sampled based on \code{\link{runif}}. Default to 1 (100\% deviation of all parameters).
#' @param con.tol The absolute convergence tolerance of each fitting run (see Details). Default is set to 0.1.
#' @param control.optim Control parameters passed along to \code{optim}. For more details, see \code{\link{optim}}.
#' @param save.rel.diff A numeric value indicating when to overwrite a pre-existing result. Default to 0.01, which means that results get overwritten only if an improvement larger than 1\% of the pre-existing value is made.
#' @param ... Other parameters to be passed on to optim.
#'
#' @return Returns the fitted parameter set and the corresponding log-likelihood.
#' @export
#'
#' @examples
#' #define cost function
#' cost_function <- function(parms){
#'   y <- parms[1] + parms[2]*c(1:3) + parms[3]^2 *c(1:3)
#'   LL <- sum((y - c(1:3))^2)
#' }
#'
#' #create profile values
#' point.profile(no.fit = c(p1 = 1),
#'               parms = c(p1 = 1, p2 = 3, p3 = 2),
#'               fit.fn = cost_function,
#'               optim.runs = 1)

point.profile <- function(no.fit,
                          parms,
                          fit.fn,
                          homedir = getwd(),
                          optim.runs = 5,
                          random.borders = 1,
                          con.tol = 0.1,
                          control.optim = list(maxit = 1000),
                          save.rel.diff = 0.01,
                          ...) {

  #sink output to a log file
  sink(paste0(homedir, "/Profile-Results/LogFiles/Log", paste0(names(no.fit[1]),"_", no.fit[1]), ".txt"), split = TRUE)
  fit.vector <- parms[-which(names(parms) %in% names(no.fit))]

  all.par <- c(no.fit, fit.vector)
  #number of succesful runs
  k <- 1
  #number of insuccessful runs
  total.tries <- 0

  #try k = optim.runs different combinations
  while(k <= optim.runs && total.tries < (4*optim.runs)){
    #set failing status to FALSE
    abort <- F

    #print number of successful runs
    cat(paste0("\nFitting run # ", k, "\n"))

    # update status file to number of current fit run
    saveRDS(object = k, file = paste0(homedir, "/Profile-Results/Status/status", paste0(names(no.fit[1]),"_", no.fit[1]), ".rds"))

    #check if initial parameters set is working
    if(k == 1){
      #take the parameter combination from the currently best model
      ran.par <- fit.vector
      #check if it works in the ls2 function
      print(ran.par)

      works <- is.finite(unite.and.fit(par = ran.par,
                                       no.fit = no.fit,
                                       par.names = names(parms),
                                       fit.fn = fit.fn,
                                       ...))
      #if not, skip the current set
      if(works == FALSE){
        cat("Inherited parameters do not work and are being skipped.\n")
        k <- k + 1
        total.tries <- total.tries + 1
        next
      }
    }else{
      #get random initial conditions and test if they work
      works <- FALSE
      while(works == FALSE){
        if(is.vector(random.borders)){
          if(length(random.borders) == 1){
            random.min <- fit.vector - random.borders*abs(fit.vector)
            random.max <- fit.vector + random.borders*abs(fit.vector)
          }else{
            random.min <- fit.vector - random.borders[-no.fit]*abs(fit.vector)
            random.max <- fit.vector + random.borders[-no.fit]*abs(fit.vector)
          }

        }else if(is.matrix(random.borders)){
          if(nrow(random.borders) == 1){
            random.min <- rep(random.borders[1,1], length(fit.vector))
            random.max <- rep(random.borders[1,2], length(fit.vector))
          }else{
            random.min <- random.borders[-no.fit,1]
            random.max <- random.borders[-no.fit,2]
          }
        }else{
          stop("random.borders must be a number, a vector or a matrix!")
        }
        ran.par <- stats::runif(n = length(fit.vector),
                                min = random.min,
                                max = random.max)
        names(ran.par) <- names(fit.vector)
        works <- is.finite(unite.and.fit(par = ran.par,
                                         no.fit = no.fit,
                                         par.names = names(parms),
                                         fit.fn = fit.fn,
                                         ...))
      }
      #print parameter sets for the log files
      print(ran.par)
    }
    #specify parameters (only for the first run)
    opt.run <- 10
    opt.previous <- 100 # the difference of opt.run and opt.previous has to be bigger than 0.1, but the values are not important
    runs = 1

    while(abs(opt.run - opt.previous) > con.tol){#test if the current run yielded better results than the previous. If yes keep optimising
      # get initial parameter sets for optim
      if(runs == 1){
        opt.par <- ran.par
      }
      if(runs > 1){
        opt.previous <- opt$value
      }
      #update run
      runs = runs + 1

      #run optim
      opt <- stats::optim(par = opt.par,
                          fn = unite.and.fit,
                          no.fit = no.fit,
                          par.names = names(parms),
                          fit.fn = fit.fn,
                          ...,
                          control = control.optim)
      opt.run <- opt$value
      opt.par <- opt$par

      if (opt.run > 1e34) {
        abort <- T
        cat("optim failed. Run skipped.\n")
        total.tries <- total.tries + 1
        break
      }

      cat(paste(opt.run, "\n"))
    }
    #options(warn = 0)
    #test if current run was better
    if(k == 1 || opt.run < opt.min){
      opt.min <- opt.run
      min.par <- opt.par

      #get corresponding parameter values
      out.par  <- rep(0, length(parms))
      names(out.par) <- names(parms)
      out.par[names(min.par)] <- min.par
      out.par[names(no.fit)] <- no.fit


      result <- c(LL = opt.run, out.par)

      #saves progress if the recent fit is the first or better than any previously saved one
      #check if this model has already been tested
      if(file.exists(paste0(homedir, "/Profile-Results/Fits/",paste0(names(no.fit[1]), "_", no.fit[1]), ".rds")) == TRUE){
        result_old <- readRDS(file = paste0(homedir,
                                            "/Profile-Results/Fits/",
                                            paste0(names(no.fit[1]),"_", no.fit[1]),
                                            ".rds"))
        result_new <- result[which(names(result) == 'LL')]
        #overwrite output file only if current fit is better than the previous one
        if(result_old[which(names(result_old) == 'LL')] > (1+ sign(result_new)*save.rel.diff)*result_new){
          cat("Current fit better than previous best fit. Results overwritten.\n")
          saveRDS(object = result, file = paste0(homedir,
                                                 "/Profile-Results/Fits/",
                                                 paste0(names(no.fit[1]),"_", no.fit[1]),
                                                 ".rds"))
        }

      }else{
        saveRDS(object = result, file = paste0(homedir,
                                               "/Profile-Results/Fits/",
                                               paste0(names(no.fit[1]), "_",no.fit[1]),
                                               ".rds"))
        cat("No results file present. Current parameters saved.\n")
      }


    }


    if(!abort || k==1) {
      k <- k + 1
      total.tries <- total.tries + 1
    }


  }

  #update the status file that corresponds with this model for the progression of the main routine
  status <-"done"
  saveRDS(object = status,
          file = paste0(homedir,
                        "/Profile-Results/Status/status",
                        paste0(names(no.fit[1]),"_", no.fit[1]),
                        ".rds"))
  cat("\nStatus file updated. Fitting done.\n")
  sink()

  return(result)
}


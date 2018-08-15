#' Smooth existing profiles
#'
#' Smooths existing profiles by refitting bad estimates. To this end, \code{smooth.profile} checks for neighbouring points that have an unusual large difference in log-likelihood and spikes, i.e. points that have higher values than both of the neighbouring points
#' @param which.par A vector containing the names of all parameters for which the respective profiles should be smoothed. Alternatively, supplying "all.par" smooths all existing profiles.
#' @param fit.fn A cost function. Has to take the complete parameter vector as an input (needs to be names \code{parms}) and must return the corresponding negative log-likelihood (-2LL, see Burnham and Anderson 2002).
#' @param threshold A numeric value determining the minimal difference between two neighbouring points that leads to refitting. Alternatively, threshold can be set to "auto" (default), which chooses a minimal difference automatically (this is calculated by dividing the difference between maximal and minimal values by the number of profile points).
#' @param spike.min A numeric value determining the minimal difference for detecting spikes. Default to 0.01.
#' @param do.not.fit The names of the parameter that are not to be fitted. Can only be supplied if a single profile is smoothed
#' @param homedir The directory to which the results should be saved to. Default to \code{\link{getwd}}().
#' @param optim.runs The number of times that each model will be fitted by \code{\link{optim}}. Default to 5.
#' @param random.borders The ranges from which the random initial parameter conditions for all \code{optim.runs} larger than one are sampled. Can be either given as a vector containing the relative deviations for all parameters or as a matrix containing in its first column the lower and in its second column the upper border values. Parameters are uniformly sampled based on \code{\link{runif}}. Default to 1 (100\% deviation of all parameters). Alternatively, functions such as \code{\link{rnorm}}, \code{\link{rchisq}}, etc. can be used if the additional arguments are passed along as well.
#' @param refit If TRUE, previously fitted ranges will be fitted again and results will be overwritten according to the value set in \code{save.rel.diff}. Default to FALSE.
#' @param con.tol The absolute convergence tolerance of each fitting run (see Details). Default is set to 0.1.
#' @param control.optim Control parameters passed along to \code{optim}. For more details, see \code{\link{optim}}.
#' @param save.rel.diff A numeric value indicating a relative threshold when to overwrite a pre-existing result. Default to 0, which means that results get overwritten if an improvement is made.
#' @param future.off Logical. If TRUE, \code{\link{future}} will not be used to calculate the results. Default to FALSE.
#' @param ... Additional parameters that can be passed along to \code{\link{future}} or \code{fit.fn}.
#'
#' @return NULL. Saves the refitted profiles in the folder "Profiles-Result/Tables/".
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
#' #create profiles
#' res <- create.profile(which.par = "p1",
#'                       par.names = inits,
#'                       range = list(seq(0, 2, 0.2)),
#'                       fit.fn = cost_function,
#'                       homedir = getwd(),
#'                       delete.old = TRUE,
#'                       x.vals = x.values,
#'                       y.vals = y.values,
#'                       sd.y = sd.y.values)
#'
#' #add noise to profile
#' res[,1] <-  res[,1] + runif(n = nrow(res), min = 0, max = 5)
#' saveRDS(res, paste0(getwd(), "/Profile-Results/Tables/p1.rds"))
#'
#' #smooth profile
#' smooth.profile(which.par = "p1",
#'                fit.fn = cost_function,
#'                homedir = getwd(),
#'                optim.runs = 1,
#'                future.off = TRUE,
#'                x.vals = x.values,
#'                y.vals = y.values,
#'                sd.y = sd.y.values)

smooth.profile <- function(which.par, fit.fn, threshold = "auto", spike.min = 0.01, do.not.fit = NULL, homedir = getwd(), optim.runs = 5, random.borders = 1, refit = F, con.tol = 0.1, control.optim = list(maxit = 1000), save.rel.diff = 0, future.off = F, ...){

  if( (length(which.par) > 1 || which.par[1] == "all.par") && is.null(do.not.fit) == FALSE){
    stop("Use 'do.not.fit' only with a single entry of 'which.par'!")
  }

  if(which.par[1] == "all.par"){
    filenames <- list.files(paste0(homedir,"/Profile-Results/Tables"), pattern="*.rds", full.names=FALSE)
    which.par <- as.vector(gsub(".rds", "", filenames))
  }

  #set graphical parameters
  nrows = round(sqrt(length(which.par)))
  ncols = ifelse(length(which.par) == 1, 1, ceiling(sqrt(length(which.par))))
  graphics::par(mfrow = c(nrows, ncols))

  #storage parameters
  save.data <- list()

  #plot the unsmoothed profiles
  for(s in 1:length(which.par)){
    if(file.exists(paste0(homedir, "/Profile-Results/Tables/", which.par[s], ".rds")) == FALSE){
      stop(paste0("The results table for parameter ", which.par[s], " is missing."))
    }
    save.data[[s]] <- readRDS(paste0(homedir, "/Profile-Results/Tables/", which.par[s], ".rds"))

    #plot data

    graphics::plot(save.data[[s]][, which.par[s]],
                   save.data[[s]][,1],
                   type = "l",
                   lwd = 3,
                   xlab = which.par[s],
                   ylab = "-2LL",
                   main = paste0("Unsmoothed profile of parameter ", which.par[s]))
    graphics::abline(h = min(save.data[[s]][,1]) + 3.84, lwd = 3, lty = 2, col = "red")
  }

  #set iteration index
  iteration <- 1
  improvement <- rep(1, length(which.par))
  which.improved <- list(rep(1, nrow(save.data[[1]])))
  if(length(which.par) > 1){
    for(s in 2:length(which.par)){
      which.improved[[s]] <- rep(1, nrow(save.data[[s]]))
    }
  }
  save.steepcliff <- list()
  save.col.pl <- list()

  while(any(improvement != 0)){
    print(paste0("ITERATION ", iteration))
    for(s in which(improvement != 0)){
      #get parameter range
      print(paste0("Analysing profile table for parameter ", which.par[s]))
      data <- save.data[[s]]
      col.pl <- which(names(data)[2:ncol(data)] == which.par[s])
      range <- data[, col.pl + 1]
      if(is.null(do.not.fit) == FALSE){
        col.pl <- c(col.pl, which(names(data)[2:ncol(data)] %in% do.not.fit))
      }
      save.col.pl[[s]] <- col.pl
      par.range <- 2:ncol(data)
      par.range <- par.range[-c(col.pl-1)]



      if(threshold == "auto"){
        cliff.min <- (max(data[,1]) - min(data[,1]))/nrow(data)
      }else{
        cliff.min <- threshold
      }
      #run until refitting of cliffs does not yield any improvement and no more spikes are present
      ###MAIN LOOP####
      ###SMOOTHING####
      #find steep cliffs and specify the type: 1 - cliff to the left, 2- cliff to the right, 12 - spiked value
      sel.steepcliff <- rep(0, length(range))
      #cycle through generated profile
      for(i in 1:length(range)){
        #check left border
        if(i==1){
          if((data[i,1] - data[i+1,1]) > cliff.min){
            sel.steepcliff[i] <- 2
          }
        }else if(i==length(range)){#check right border
          if((data[i,1] - data[i-1,1]) > cliff.min){
            sel.steepcliff[i] <- 1
          }
        }else{#check all values inbetween
          if((data[i,1] - max(data[c(i-1,i+1),1])) > spike.min){#if value is a spike
            sel.steepcliff[i] <- 12
          }else if ((data[i,1] - data[i-1,1]) > cliff.min){#if cliff is to the left
            sel.steepcliff[i] <- 1
          }else if ((data[i,1] - data[i+1,1]) > cliff.min){#if cliff is to the right
            sel.steepcliff[i] <- 2
          }

        }

      }
      #get problematic entries
      steepcliff <- which(sel.steepcliff > 0)
      save.steepcliff[[s]] <- steepcliff
      print(paste0("Found ", length(steepcliff)," unsuitable value(s) in the profile (spike or difference to neighbouring points larger than ", format(cliff.min, digits = 2), ")."))
      print(data[steepcliff, col.pl + 1])

      #fit only if cliffs are present
      if(length(steepcliff) > 0){
        #cycle through cliffs
        for(i in 1:length(steepcliff)){
          k <- steepcliff[i]
          cliff <- sel.steepcliff[k]
          fixed.par <- data[k, col.pl + 1]
          names(fixed.par) <- names(data)[col.pl + 1]
          #set seed for random fitting
          set.seed(unclass(Sys.time()))

          #use the parameter values of the better side of the cliff
          if(cliff == 1){
            params.test <- as.numeric(data[k - 1,-1])
            names(params.test) <- names(data[k - 1,-1])
          }else if(cliff == 2){
            params.test <- as.numeric(data[k + 1,-1])
            names(params.test) <- names(data[k + 1,-1])
          }else{#cliff == 12
            params.test <- as.numeric(data[k - 1,-1])
            names(params.test) <- names(data[k - 1,-1])

            params.test2 <- as.numeric(data[k + 1,-1])
            names(params.test2) <- names(data[k + 1,-1])
          }

          works <- is.finite(unite.and.fit(par = params.test[-col.pl],
                                           no.fit = fixed.par,
                                           par.names = names(data)[2:ncol(data)],
                                           fit.fn = fit.fn,
                                           ...))
          #if parameter combination failed, try 100 times to create a working one by jittering
          if(works == FALSE){
            #save old parameter set
            params.base <- params.test
            #define count variable
            count.test <- 0
            while(works == FALSE & count.test < 100){

              if(is.vector(random.borders)){
                random.min <- params.test - random.borders*abs(params.test)
                random.max <- params.test + random.borders*abs(params.test)
              }else if(is.matrix(random.borders)){
                if(nrow(random.borders) == 1){
                  random.min <- rep(random.borders[1,1], length(params.test))
                  random.max <- rep(random.borders[1,2], length(params.test))
                }else{
                  random.min <- random.borders[,1]
                  random.max <- random.borders[,2]
                }
              }else{
                stop("random.borders must be a number, a vector or a matrix!")
              }

              #create new parameter set for testing
              params.test <- stats::runif(n = length(params.test),
                                          min = random.min,
                                          max = random.max)
              names(params.test) <- names(params.base)
              #test the set
              works <- is.finite(unite.and.fit(par = params.test[-col.pl],
                                               no.fit = fixed.par,
                                               par.names = names(data)[2:ncol(data)],
                                               fit.fn = fit.fn,
                                               ...))
              #update count variable
              count.test <- count.test + 1
            }
          }
          #run optimization only if parameter combination works
          if(works == TRUE){

            # #optimise
            # p.scale <- par.scales(par = abs(params.test),
            #                       scale = abs(params.test),
            #                       fix = length(params.test))

            if(future.off == TRUE){
              point.profile(no.fit = fixed.par,
                            parms = params.test,
                            fit.fn = fit.fn,
                            homedir = homedir,
                            optim.runs = optim.runs,
                            random.borders = random.borders,
                            con.tol = con.tol,
                            control.optim = control.optim,
                            save.rel.diff = save.rel.diff,
                            ...)
            }else{
              future::future(point.profile(no.fit = fixed.par,
                                           parms = params.test,
                                           fit.fn = fit.fn,
                                           homedir = homedir,
                                           optim.runs = optim.runs,
                                           random.borders = random.borders,
                                           con.tol = con.tol,
                                           control.optim = control.optim,
                                           save.rel.diff = save.rel.diff,
                                           ...),
                             label = paste0(which.par[s], "_", fixed.par[1]),
                             ...)
            }
          }
          #if the value is spiked, try parameter combination of the right side
          if(cliff == 12){

            #try and see if parameter combination work
            works <- is.finite(unite.and.fit(par = params.test2[-col.pl],
                                             no.fit = fixed.par,
                                             par.names = names(data)[2:ncol(data)],
                                             fit.fn = fit.fn,
                                             ...))
            #if parameter combination failed, try 100 times to create a working one by jittering
            if(works == FALSE){
              #save old parameter set
              params.base <- params.test2
              #define count variable
              count.test <- 0
              while(works == FALSE & count.test < 100){
                #create new parameter set for testing
                if(is.vector(random.borders)){
                  random.min <- params.test2 - random.borders*abs(params.test2)
                  random.max <- params.test2 + random.borders*abs(params.test2)
                }else if(is.matrix(random.borders)){
                  if(nrow(random.borders) == 1){
                    random.min <- rep(random.borders[1,1], length(params.test2))
                    random.max <- rep(random.borders[1,2], length(params.test2))
                  }else{
                    random.min <- random.borders[,1]
                    random.max <- random.borders[,2]
                  }
                }else{
                  stop("random.borders must be a number, a vector or a matrix!")
                }

                #create new parameter set for testing
                params.test <- stats::runif(n = length(params.test2),
                                            min = random.min,
                                            max = random.max)
                names(params.test2) <- names(params.base)
                #test the set
                works <- is.finite(unite.and.fit(par = params.test2[-col.pl],
                                                 no.fit = fixed.par,
                                                 par.names = names(data)[2:ncol(data)],
                                                 fit.fn = fit.fn,
                                                 ...))
                #update count variable
                count.test <- count.test + 1
              }
            }
            #run optimization only if parameter combination works
            if(works == TRUE){

              if(future.off == TRUE){
                point.profile(no.fit = fixed.par,
                              parms = params.test2,
                              fit.fn = fit.fn,
                              homedir = homedir,
                              optim.runs = optim.runs,
                              random.borders = random.borders,
                              con.tol = con.tol,
                              control.optim = control.optim,
                              save.rel.diff = save.rel.diff,
                              ...)
              }else{
                future::future(point.profile(no.fit = fixed.par,
                                             parms = params.test2,
                                             fit.fn = fit.fn,
                                             homedir = homedir,
                                             optim.runs = optim.runs,
                                             random.borders = random.borders,
                                             con.tol = con.tol,
                                             control.optim = control.optim,
                                             save.rel.diff = save.rel.diff,
                                             ...),
                               label = paste0(which.par[s], "_", fixed.par[1]),
                               ...)
              }

            }
          }

        }
      }

    }
    #read in data
    for(k in which(improvement != 0)){
      #define improvement variable
      improvement[k] <- 0
      for(j in 1:length(save.steepcliff[[k]])){
        while(file.exists(paste0(homedir, "/Profile-Results/Status/status", which.par[k],"_", save.data[[k]][save.steepcliff[[k]][j], save.col.pl[[k]] + 1], ".rds")) == FALSE){
          Sys.sleep(5)
        }
        res <- readRDS(paste0(homedir, "/Profile-Results/Fits/", which.par[k],"_", save.data[[k]][save.steepcliff[[k]][j], save.col.pl[[k]] + 1], ".rds"))
        if(res[1] < save.data[[k]][save.steepcliff[[k]][j], 1]){
          save.data[[k]][save.steepcliff[[k]][j], ] <- res
          improvement[k] <- improvement[k] + 1
          if(j == 1){
            slct <- c(1,2)
          }else if(j == nrow(save.data[[k]])){
            slct <- c(j-1, j)
          }else{
            slct <- c(j-1, j, j +1)
          }
          which.improved[[k]][slct] <- 1
        }else{
          which.improved[[k]][j] <- 0
        }
      }

      if(improvement[k] > 0){
        print(paste0("Improved ", improvement[k], " out of ", length(steepcliff), " profile values for parameter ",which.par[s],"."))
        saveRDS(save.data[[k]], paste0(homedir, "/Profile-Results/Tables/", which.par[k], ".rds"))
      }else{
        print(paste0("Improved ", improvement[k], " out of ", length(steepcliff), " profile values for parameter ",which.par[s],". Smoothing terminated for this parameter."))
      }

    }

    graphics::par(mfrow = c(nrows, ncols))
    for(s in 1:length(which.par)){
      graphics::plot(save.data[[s]][, which.par[s]],
                     save.data[[s]][,1],
                     type = "l",
                     lwd = 3,
                     xlab = which.par[s],
                     ylab = "-2LL",
                     main = paste0("Smoothed profile, iteration ", iteration))
      graphics::abline(h = min(save.data[[s]][,1]) + 3.84, lwd = 3, lty = 2, col = "red")
    }

    grDevices::pdf(file = paste0(homedir, "/Profile-Results/Figures/Profiles", ".pdf"),
                   width  = 4*ncols,
                   height = 4*nrows,
                   useDingbats = F)
    graphics::par(mfrow = c(nrows, ncols))
    for(s in 1:length(which.par)){
      graphics::plot(save.data[[s]][, which.par[s]],
                     save.data[[s]][,1],
                     type = "l",
                     lwd = 3,
                     xlab = which.par[s],
                     ylab = "-2LL",
                     main = paste0("Smoothed profile, iteration ", iteration))
      graphics::abline(h = min(save.data[[s]][,1]) + 3.84, lwd = 3, lty = 2, col = "red")
    }
    grDevices::dev.off()

    iteration <- iteration + 1

    if(any(improvement > 0) == FALSE){
      print("No further improvement in any profiles was achieved. Terminating.")
    }
  }
}

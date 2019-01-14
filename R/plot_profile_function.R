#' Plot Calculated Profiles
#'
#'Plots the previously calculated profiles of \code{\link{create.profile}}.
#' @param which.par Name of the parameter that a profile should be obtained for. If multiple profiles are supposed to be calculated, a vector can be passed along as well. Alternatively, supplying "all.par" calculates the profiles for all available parameters.
#' @param homedir The directory to which the results should be saved to. Default to \code{\link{getwd}}().
#' @param conf.level The confidence level used for calculating the confidence intervals. Default to 0.95.
#' @param ... Additional arguments that can be passed on to \code{\link{plot}}.
#'
#' @return Plots the profiles and returns the confidence intervals
#' @export
#'
#' @examples
#' #' #create data with standard deviation of 1
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
#'
#' plotting.profile("all.par")
plotting.profile <- function(which.par, homedir = getwd(), conf.level = 0.95, ...){

  old.par <- graphics::par("mfrow")
  on.exit(graphics::par(mfrow = old.par))

  if(which.par[1] == "all.par"){
    filenames <- list.files(paste0(homedir,"/Profile-Results/Tables"), pattern="*.rds", full.names=FALSE)
    which.par <- as.vector(gsub(".rds", "", filenames))
  }

  overall.min <- Inf
  store <- list()
  for(s in 1:length(which.par)){
    if(file.exists(paste0(homedir, "/Profile-Results/Tables/", which.par[s], ".rds")) == FALSE){
      stop(paste0("The results table for parameter ", which.par[s], " is missing."))
    }
    data <- readRDS(paste0(homedir, "/Profile-Results/Tables/", which.par[s], ".rds"))
    if(overall.min > min(data$LL)){
      overall.min <- min(data$LL)
    }
    store[[s]] <- data

  }

  cut.off <- overall.min + stats::qchisq(p = conf.level, df = 1)

  nrows = sqrt(length(which.par))
  ncols = ifelse(length(which.par) == 1, 1, ceiling(sqrt(length(which.par))))

  store.conf <- c()
  dots <- list(...)

  graphics::par(mfrow = c(nrows, ncols))

  for(i in 1:length(which.par)){
    table.x <- store[[i]]
    conf <- c(min(which(table.x$LL <= cut.off)), max(which(table.x$LL <= cut.off)))
    conf.interval <- c(min = min(table.x[which(table.x$LL == min(table.x$LL)), which.par[i]]),
                       lower.ci = table.x[conf[1], which.par[i]],
                       upper.ci = table.x[conf[2], which.par[i]])
    store.conf <- rbind(store.conf, conf.interval)

    do.call(graphics::plot, c(list(x = table.x[, which.par[i]],
                                  y = table.x$LL,
                                  type = "l",
                                  lwd = 3,
                                  xlab = which.par[i],
                                  ylab = "-2LL",
                                  main = paste0("Profile likelihood of ", which.par[i])),
                              dots))

    graphics::abline(h = cut.off, lwd = 3, lty = 2, col = "red")
    graphics::abline(h = overall.min, lty = 3, col = "grey")

  }
  rownames(store.conf) <- which.par
  return(store.conf)

}

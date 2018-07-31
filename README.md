
<!-- README.md is generated from README.Rmd. Please edit that file -->
ProfileIroning
==============

The ProfileIroning package is intended to calculate confidence intervals for model parameters based on a profile likelihood approach. Additionally, it is able to smooth existing profiles by refitting specific entries in the profiles

Installation
------------

You can install ProfileIroning from github with:

``` r
# install.packages("devtools")
devtools::install_github("GabelHub/ProfileIroning")
# alternative installation command
devtools::install_git("git://github.com/GabelHub/ProfileIroning.git", branch = "master")
```

Main functions
--------------

This function creates the profiles for specified parameters, based on a user-defined cost-function. Profiles are calculated by using the *future* package which allows to parallelise the whole process. Depending on the preference, *create.profile* can already extensively try to achieve a good profile. Alternatively, the generated profiles can be smoothed afterwards by calling *smooth.profile*, which attempts to reduce the number of badly fitted parameter sets. In the end, *plotting.profile* can be called to display the profiles and to return the respective confidence intervals. An example code is given here

``` r
library(ProfileIroning)

#create data with standard deviation of 1
x.values <- 1:7
y.values <-  9 * x.values^2 - exp(2 * x.values)
sd.y.values <- rep(1,7)

#define initial parameter values
inits <- c(p1 = 3, p2 = 4, p3 = -2, p4 = 2, p5 = 0)

#define cost function that returns the negative log-likelihood
cost_function <- function(parms, x.vals, y.vals, sd.y){
  # restrict the search range to -5 to +5
  if(max(abs(parms)) > 5){
    return(NA)
  }
  with(as.list(c(parms)), {
    res <- p1*4 + p2*x.vals + p3^2*x.vals^2 + p4*sin(x.vals)  - exp(p5*x.vals)
    diff <- sum((res - y.vals)^2/sd.y)
  })
}

dont.fit <- c(p1 = 3)
#perform model selection
#perform model selection
res <- create.profile(which.par = "p1",
                      par.names = inits,
                      range = list(seq(0, 2, 0.2)),
                      fit.fn = cost_function,
                      homedir = getwd(),
                      x.vals = x.values,
                      y.vals = y.values,
                      sd.y = sd.y.values)

res[,1] <-  res[,1] + runif(n = nrow(res), min = 0, max = 5)
saveRDS(res, paste0(getwd(), "/Profile-Results/Tables/p1.rds"))

smooth.profile(which.par = "p1",
               fit.fn = cost_function,
               homedir = getwd(),
               optim.runs = 1,
               future.off = TRUE,
               x.vals = x.values,
               y.vals = y.values,
               sd.y = sd.y.values)

plotting.profile(which.par = "p1")
```

context("Test the optimising function point.profile")

system("rm -r tests/testthat/Profile-Results")
system("rm -r Profile-Results")

library(ProfileIroning)

set.seed(1)

inits <- c(p1 = 3, p2 = 4, p3 = -2, p4 = 2, p5 = 0)

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
res <- create.profile(which.par = "all.par",
                      par.names = inits,
                      range = list(seq(0, 2, 0.1),
                                   seq(0, 5, 1),
                                   seq(2.9, 3.1, 0.01),
                                   seq(0, 3, 0.1),
                                   seq(1.999999, 2.000001, 0.0000001)),
                      fit.fn = cost_function,
                      homedir = getwd(),
                      optim.runs = 5,
                      random.borders = 1,
                      con.tol = 0.1,
                      control.optim = list(maxit = 1000),
                      x.vals = x.values,
                      y.vals = y.values,
                      sd.y = sd.y.values)

# test_that("Gives the correct output", {
#   expect_equal(as.numeric(round(res["LL"])), 52)
#
# })


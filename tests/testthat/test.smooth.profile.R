context("Test the optimising function smooth.profile")

system("rm -r tests/testthat/Profile-Results")
system("rm -r Profile-Results")

library(ProfileIroning)

set.seed(1)

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

#perform model selection
res <- create.profile(which.par = c("p1","p2","p3"),
               par.names = inits,
               range = list(seq(0, 2, 0.2), seq(0, 2, 0.2), seq(0, 2, 0.2)),
               fit.fn = cost_function,
               optim.runs = 1,
               homedir = getwd(),
               x.vals = x.values,
               y.vals = y.values,
               sd.y = sd.y.values)

res[[1]][,1] <-  res[[1]][,1] + runif(n = nrow(res[[1]]), min = 0, max = 5)
saveRDS(res[[1]], paste0(getwd(), "/Profile-Results/Tables/p1.rds"))

res[[2]][,1] <-  res[[2]][,1] + runif(n = nrow(res[[2]]), min = 0, max = 5)
saveRDS(res[[2]], paste0(getwd(), "/Profile-Results/Tables/p2.rds"))

res[[3]][,1] <-  res[[3]][,1] + runif(n = nrow(res[[3]]), min = 0, max = 5)
saveRDS(res[[3]], paste0(getwd(), "/Profile-Results/Tables/p3.rds"))

smooth.profile(which.par = "all.par",
               fit.fn = cost_function,
               homedir = getwd(),
               optim.runs = 1,
               future.off = FALSE,
               x.vals = x.values,
               y.vals = y.values,
               sd.y = sd.y.values)

test_that("Gives the correct output", {
  expect_equal(as.numeric(round(res[[1]][1,"LL"])), 1)

})


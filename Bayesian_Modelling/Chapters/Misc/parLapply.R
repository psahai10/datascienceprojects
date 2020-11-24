library(parallel)
library(MASS)

x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- seq(1, 10000)
boot_fx <- function(trial) {
  ind <- sample(100, 100, replace=TRUE)
  result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
  r <- coefficients(result1)
  res <- rbind(data.frame(), r)
}
system.time({
  results <- lapply(trials, boot_fx)
})

library(parallel)
cl <- makeCluster(mc <- getOption("cl.cores", 4))
#https://stackoverflow.com/questions/12019638/using-parallels-parlapply-unable-to-access-variables-within-parallel-code
clusterExport(cl=cl, varlist=c("x"))

system.time({
  parLapply(cl, trials, boot_fx)
})

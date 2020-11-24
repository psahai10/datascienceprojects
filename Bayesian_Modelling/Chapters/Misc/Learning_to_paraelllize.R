library(parallel)
library(MASS)
library("rqdatatable")
library("microbenchmark")
library("ggplot2")
library("WVPlots")
library(furrr)
library(purrr)
library(data.table)
library(brms)

d <- data.table(yc = rnorm(100, 0, 1),
                yt = rnorm(100, .5, 1),
                id = 1:100) %>% 
  melt(id.vars = "id", measure.vars = c("yc","yt"), variable.name = "group")

fit <- brm(data = d,
           family = gaussian,
           value ~ 0 + intercept + group,
           prior = c(prior(normal(0, 10), class = b),
                     prior(student_t(3, 1, 10), class = sigma)),
           seed = 1)

sim_d <- function(seed, n) {
  set.seed(seed)
  
  data.table(group = rep(c("control", "treatment"), n/2)) %>% 
    .[, value := ifelse(group == "control", rnorm(n), rnorm(n, 0.5))]
}

n_sim <- 5
n <- 4

system.time({
sims_norm <- data.table(seed = 1:n_sim) %>% 
  .[, d := map(seed, sim_d, n = n)] %>% 
  .[, fit := map2(d, seed, ~ update(fit, newdata = .x, seed = .y))]
})

system.time({
  sims_norm <- data.table(seed = 1:n_sim) %>% 
    .[, d := future_map(seed, sim_d, n = n)] %>% 
    .[, fit := future_map2(d, seed, ~ update(fit, newdata = .x, seed = .y))]
})


sims_wrapr <- wrapr::execute_parallel(
  tables = list(data = data, 
                annotation = annotation),
  f = dt_f,
  partition_column = "key_group",
  cl = cl)



timings <- microbenchmark(
  data_table_parallel = 
    nrow(data_table_parallel_f(data, annotation)),
  data_table = nrow(data_table_f(data, annotation)),
  times = 10L)

# autoplot(timings)
timings <- as.data.frame(timings)
timings$seconds <- timings$time/1e+9

ScatterBoxPlotH(timings, 
                xvar = "seconds", yvar = "expr", 
                title="task duration by method")

head(sims)

sims[1, fit]

head(sims[1, d][[1]])
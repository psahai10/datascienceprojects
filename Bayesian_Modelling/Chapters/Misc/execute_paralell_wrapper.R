# https://davisvaughan.github.io/furrr/
# https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html
# http://www.win-vector.com/blog/2018/07/speed-up-your-r-work/

library(furrr)
library(purrr)
library(parallel)
library(MASS)
library("rqdatatable")
library("microbenchmark")
library("ggplot2")
library("WVPlots")

ncore <- parallel::detectCores()
cl <- parallel::makeCluster(ncore)

set.seed(2362)
mk_example <- function(nkey, nrep, ngroup = 20) {
  keys <- paste0("key_", seq_len(nkey))
  key_group <- sample(as.character(seq_len(ngroup)), 
                      length(keys), replace = TRUE)
  names(key_group) <- keys
  key_table <- data.frame(
    key = rep(keys, nrep),
    stringsAsFactors = FALSE)
  key_table$data <- runif(nrow(key_table))
  instance_table <- data.frame(
    key = rep(keys, nrep),
    stringsAsFactors = FALSE)
  instance_table$id <- seq_len(nrow(instance_table))
  instance_table$info <- runif(nrow(instance_table))
  # groups should be no finer than keys
  key_table$key_group <- key_group[key_table$key]
  instance_table$key_group <- key_group[instance_table$key]
  list(key_table = key_table,
       instance_table = instance_table)
}

dlist <- mk_example(10, 10)
data <- dlist$instance_table
annotation <- dlist$key_table

data_table_f <- function(data, annotation) {
  data <- data.table::as.data.table(data)
  annotation <- data.table::as.data.table(annotation)
  joined <- merge(data, annotation, 
                  by = "key", 
                  all=FALSE, 
                  allow.cartesian=TRUE)
  joined <- joined[joined$data <= joined$info, ]
  data.table::setorderv(joined, cols = "data")
  joined <- joined[, .SD[.N], id]
  data.table::setorderv(joined, cols = "id")
}
resdt <- data_table_f(data, annotation)
head(resdt)

parallel::clusterEvalQ(cl, library("data.table"))

parallel::clusterExport(cl, "data_table_f")

dt_f <- function(tables_list) {
  data <- tables_list$data
  annotation <- tables_list$annotation
  data_table_f(data, annotation)
}

data_table_parallel_f <- function(data, annotation) {
  respdt <- wrapr::execute_parallel(
    tables = list(data = data, 
                  annotation = annotation),
    f = dt_f,
    partition_column = "key_group",
    cl = cl) %.>%
    data.table::rbindlist(.)
  data.table::setorderv(respdt, cols = "id")
  respdt
}
respdt <- data_table_parallel_f(data, annotation)
head(respdt)


dlist <- mk_example(300, 300)
data <- dlist$instance_table
annotation <- dlist$key_table

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

#https://rpubs.com/kaz_yos/stan_count2

library(tidyverse)
library(dplyr)
## devtools::install_github('jburos/biostan', build_vignettes = TRUE, dependencies = TRUE)
library(rstan)
## To avoid recompilation of unchanged Stan programs, we recommend calling
rstan_options(auto_write = TRUE)
##
library(bayesplot)
##
library(loo)
## Seed from random.org
set.seed(673788956)

data1 <- haven::read_dta("https://stats.idre.ucla.edu/stat/stata/dae/nb_data.dta")
head(data1)

data <- data1 %>%
  mutate(prog = factor(prog, levels=1:3, labels = c("General", "Academic", "Vocational")),
         id = factor(id))

## For counting
data_count <- data %>%
  rename(y = daysabs) %>%
  count(y)
data_count

data %>%
  ggplot(mapping = aes(x = daysabs, fill = prog)) +
  geom_bar(stat = "count") +
  facet_grid(prog ~ ., margin = TRUE, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())

data %>%
  group_by(prog) %>%
  summarize(mean = mean(daysabs),
            var = var(daysabs))

X <- model.matrix(object = daysabs ~ math + prog, data = data)
head(X, n = 10)

##
extract_post_pred <- function(stan_fit) {
  tidybayes::tidy_draws(stan_fit) %>%
    select(.chain, .iteration, .draw, starts_with("y_new")) %>%
    gather(key = key, value = value, starts_with("y_new")) %>%
    mutate(key = gsub("y_new|\\[|\\]", "", key) %>% as.integer())
}
##
plot_draws <- function(stan_fit, n_sample, data_count = data_count) {
  draw_data <- extract_post_pred(stan_fit)
  sample_indices <- sample(seq_len(max(draw_data$.iteration)), size = n_sample)
  draw_data %>%
    group_by(.chain, .iteration, .draw) %>%
    count(value) %>%
    filter(.iteration %in% sample_indices) %>%
    ggplot(mapping = aes(x = value, y = n, group = interaction(.chain, .iteration, .draw))) +
    ## Plot random draws from posterior
    geom_line(alpha = 0.5, size = 0.1) +
    ## Include actual data distribution
    geom_line(data = data_count, color = "gray", alpha = 0.7, size = 2,
              mapping = aes(x = y, group = NULL)) +
    facet_wrap( ~ .chain) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())
}

stan_code_poisson <- "
data { 
  real<lower=0> a; 
  real<lower=0> b; 
  int<lower=0> N; 
  int<lower=0> M; 
  matrix[N, M] X; 
  int<lower=0> y[N]; 
} 
  
parameters { 
  real<lower=0> lambda; 
} 

model { 
  target += gamma_lpdf(lambda | a, b); 
 
     for (i in 1:N) {
         target += poisson_lpmf(y[i] | lambda); 
  } 
} 

generated quantities { 
    int<lower=0> y_new[N]; 
    vector[N] log_lik; 
      
     for (i in 1:N) { 
        y_new[i] = poisson_rng(lambda); 
        log_lik[i] = poisson_lpmf(y[i] | lambda); 
  } 
} 
"

dat_list <-list(a = 10^(-3), b = 10^(-3),
                N = nrow(X), M = ncol(X),
                X = X, y = y)

stan_model_poisson <- stan(model_code = stan_code_poisson,
                          data=dat_list, chains = 4, cores=4)

check_hmc_diagnostics(stan_model_poisson)

pars <- c("lambda","lp__")
print(stan_model_poisson, pars = pars)

pairs(stan_model_poisson, pars = pars)

traceplot(stan_model_poisson, inc_warmup = TRUE, pars = pars)

plot_draws(stan_model_poisson, n_sample = 20)

loo(stan_model_poisson)

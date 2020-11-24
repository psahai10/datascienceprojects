library(lme4)
library(brms)
data('cbpp')
d <- cbpp
head(cbpp)


fit1 <- brm(incidence | trials(size) ~ period + (1|herd), 
            data = d, family = binomial())

fit1 <- brm(data=d,
            family=binomial(),
            incidence | trials(size) ~ period + (1|herd),
            iter = 2500, warmup = 500, cores = 4, chains = 4,
            seed = 10,)

beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]"
)

stanvars <- stanvar(scode = stan_funs, block = "functions")

fit2 <- brm(
  incidence | vint(size) ~ period + (1|herd), data = cbpp, 
  family = beta_binomial2, stanvars = stanvars
)

expose_functions(fit2, vectorize = TRUE)

log_lik_beta_binomial2 <- function(i, prep) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  beta_binomial2_lpmf(y, mu, phi, trials)
}

loo(fit1, fit2)

```{r posterior_predict}
posterior_predict_beta_binomial2 <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  trials <- prep$data$vint1[i]
  beta_binomial2_rng(mu, phi, trials)
}

pp_check(fit2)

posterior_epred_beta_binomial2 <- function(prep) {
  mu <- prep$dpars$mu
  trials <- prep$data$vint1
  trials <- matrix(trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  mu * trials
}

conditional_effects(fit2, conditions = data.frame(size = 1))
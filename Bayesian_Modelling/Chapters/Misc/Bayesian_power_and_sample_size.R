# https://www4.stat.ncsu.edu/~reich/ABA/Code/Power

N           <- 100  # Number of MC reps
n           <- 100  # Sample size
sigma       <- 1    # Error standard devation

true_mu     <- 0.2  # True value of mu

pri.mn.anal <- 0    # Prior for mu
pri.sd.anal <- Inf

L <- rep(0,N)
U <- rep(0,N)

set.seed(0820)

for(rep in 1:N){
  
  #Generate data:
  Y  <- rnorm(n,true_mu,sigma)
  
  #Compute posterior 95% interval:
  post.var <- 1/(n/sigma^2+1/pri.sd.anal^2)
  post.mn  <- post.var*(pri.mn.anal/pri.sd.anal^2+sum(Y)/sigma^2)
  L[rep]   <- post.mn-1.96*sqrt(post.var)
  U[rep]   <- post.mn+1.96*sqrt(post.var)
}

plot(NA,xlim=c(-1,1),ylim=c(1,N),xlab="95% interval",ylab="Replication")
abline(v=0)

for(rep in 1:N){
  reject <- L[rep]>0 | U[rep]<0
  lines(c(L[rep],U[rep]),c(rep,rep),col=ifelse(reject>0,2,1))
}

mean(L>0 | U<0)

#https://tysonbarrett.com/jekyll/update/2019/07/21/BayesianSims/

library(tidyverse)
library(data.table)
library(brms)


theme_set(theme_dark() +
            theme(legend.position = "none",
                  panel.grid = element_blank(),
                  panel.background = element_rect(fill = "grey20"),
                  plot.background = element_rect(fill = "grey20"),
                  text = element_text(color = "white"),
                  axis.text = element_text(color = "white")))



d <- data.table(yc = rnorm(100, 0, 1),
                yt = rnorm(100, .5, 1),
                id = 1:100) %>% 
  melt(id.vars = "id", measure.vars = c("yc","yt"), variable.name = "group")


library(plyr)
mu <- ddply(d, "group", summarise, grp.mean=mean(value))
head(mu)

p <- ggplot(d, aes(x=value, fill=group)) + labs(x='', y='') + ggtitle("Sample VS Normal Distribution") +
  geom_density(alpha=.2) + xlim(-4,4) +
  geom_rug(aes(x = value, y = 0, color=group, alpha=0.7), position = position_jitter(height = 0))


p + geom_vline(data=mu, aes(xintercept=grp.mean, color=group, alpha=0.7), linetype="dashed", size=1)



head(d)

fit <- brm(data = d,
           family = gaussian,
           value ~ 0 + intercept + group,
           prior = c(prior(normal(0, 10), class = b),
                     prior(student_t(3, 1, 10), class = sigma)),
           seed = 1)
plot(fit)

fit

# devtools::install_github("mjskay/tidybayes")
library(tidybayes)
posts <- posterior_samples(fit)
posts %>% 
  ggplot(aes(b_groupyt)) +
  geom_halfeyeh(aes(y = 0), color = "grey90") +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(y = "",
       x = "Treatment Effect")

# https://davisvaughan.github.io/furrr/
# https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html
# http://www.win-vector.com/blog/2018/07/speed-up-your-r-work/
library(parallel)
library(MASS)
library("rqdatatable")
library("microbenchmark")
library("ggplot2")
library("WVPlots")


sim_d <- function(seed, n) {
  set.seed(seed)

  mu_t <- .5
  mu_c <- 0
  
  data.table(group = rep(c("control", "treatment"), n/2)) %>% 
    .[, treatment := ifelse(group == "control", 0, 1), 
        y := ifelse(group == "control",
                    rnorm(n, mean = mu_c, sd = 1), 
                    rnorm(n, mean = mu_t, sd = 1))]
}


library(furrr)
library(purrr)

n_sim <- 100
n <- 80
plan(multisession)
sims <- data.table(seed = 1:n_sim) %>% 
  .[, d := map(seed, sim_d, n = n)] %>% 
  .[, fit := map2(d, seed, ~ update(fit, newdata = .x, seed = .y))]

head(sims)

sims[1, d]
sims[1, fit]
sims[1, effects]

head(sims[1, d][[1]])

library(broom)
n_sim <- 1
n <- 10
sims <- data.table(seed = 1:n_sim) %>% 
  .[, d := map(seed, sim_d, n = n)] %>% 
  .[, fit := map2(d, seed, ~update(fit, newdata = .x, seed = .y))] %>% 
  .[, treatment := map(fit, ~tidy(.x, prob = .95))] %>% 
  .[, unnest(.SD, cols=c(treatment))] %>% 
  filter(term == "b_treatment") %>% 
  mutate(powered = case_when(lower > 0 ~ 1, TRUE ~ 0))

sims[, .N, by = powered]

sims %>% 
  ggplot(aes(x = seed, y = estimate, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_pointrange(size = .1, aes(color = factor(powered))) +
  labs(y = "Estimate (Credibility Interval)",
       x = "Index") +
  scale_color_manual(values = c("dodgerblue1", "#F1C40F"))

n_sim <- 100
n <- 140
plan(multisession)
sims140 <- data.table(seed = 1:n_sim) %>% 
  .[, d := map(seed, sim_d, n = n)] %>% 
  .[, fit := map2(d, seed, ~update(fit, newdata = .x, seed = .y))] %>% 
  .[, effects := map(fit, ~tidy(.x, prob = .95))] %>% 
  .[, unnest(.SD, effects)] %>% 
  .[term == "b_grouptreatment"] %>% 
  .[, powered := case_when(lower > 0 ~ 1, TRUE ~ 0)] %>% 
  .[, sample := 140]

n <- 120
plan(multisession)
sims120 <- data.table(seed = 1:n_sim) %>% 
  .[, d := future_map(seed, sim_d, n = n)] %>% 
  .[, fit := future_map2(d, seed, ~update(fit, newdata = .x, seed = .y))] %>% 
  .[, effects := future_map(fit, ~tidy(.x, prob = .95))] %>% 
  .[, unnest(.SD, effects)] %>% 
  .[term == "b_grouptreatment"] %>% 
  .[, powered := case_when(lower > 0 ~ 1, TRUE ~ 0)] %>% 
  .[, sample := 120]

n <- 100
plan(multisession)
sims100 <- data.table(seed = 1:n_sim) %>% 
  .[, d := future_map(seed, sim_d, n = n)] %>% 
  .[, fit := future_map2(d, seed, ~update(fit, newdata = .x, seed = .y))] %>% 
  .[, effects := future_map(fit, ~tidy(.x, prob = .95))] %>% 
  .[, unnest(.SD, effects)] %>% 
  .[term == "b_grouptreatment"] %>% 
  .[, powered := case_when(lower > 0 ~ 1, TRUE ~ 0)] %>% 
  .[, sample := 100]

sims80 <- sims[, sample := 80]

n <- 60
plan(multisession)
sims60 <- data.table(seed = 1:n_sim) %>% 
  .[, d := future_map(seed, sim_d, n = n)] %>% 
  .[, fit := future_map2(d, seed, ~update(fit, newdata = .x, seed = .y))] %>% 
  .[, effects := future_map(fit, ~tidy(.x, prob = .95))] %>% 
  .[, unnest(.SD, effects)] %>% 
  .[term == "b_grouptreatment"] %>% 
  .[, powered := case_when(lower > 0 ~ 1, TRUE ~ 0)] %>% 
  .[, sample := 60]

n <- 40
plan(multisession)
sims40 <- data.table(seed = 1:n_sim) %>% 
  .[, d := future_map(seed, sim_d, n = n)] %>% 
  .[, fit := future_map2(d, seed, ~update(fit, newdata = .x, seed = .y))] %>% 
  .[, effects := future_map(fit, ~tidy(.x, prob = .95))] %>% 
  .[, unnest(.SD, effects)] %>% 
  .[term == "b_grouptreatment"] %>% 
  .[, powered := case_when(lower > 0 ~ 1, TRUE ~ 0)] %>% 
  .[, sample := 40]

cur <- rbindlist(list(sims40, sims60, sims80, sims100, sims120, sims140))
cur %>% 
  .[, mean(powered), by = sample] %>% 
  ggplot(aes(sample, V1)) +
  geom_point(color = "grey80") +
  geom_line(color = "grey80") +
  geom_hline(yintercept = .8, color = "grey80", linetype = "dashed") +
  labs(x = "Sample Size", 
       y = "Power")

ggplot(cur, aes(seed, estimate, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_pointrange(size = .1, aes(color = factor(powered))) +
  labs(y = "Estimate (Credibility Interval)",
       x = "Index") +
  scale_color_manual(values = c("dodgerblue1", "#F1C40F")) +
  facet_wrap(~sample)

# https://cran.r-project.org/web/packages/BayesMAMS/vignettes/BayesFreq.pdf

library("BayesMAMS")
ssbayes(k=2, nu=1, q0=c(0, 0, 0), eta=0.95, zeta=0.90, deltastar=0.5, prec="known",
        crit="1")

k <- 2
alloc <- sqrt(k)
nu <- 1
deltastar <- 0.5
alpha <- 0.05
power <- 0.90

ssfreq_bon <- ((qnorm(1 - alpha/k) + qnorm(power)) / (sqrt(nu) * deltastar))^2 *
  (1 + 1/sqrt(k))
ceiling(c(sqrt(k) * ssfreq_bon, rep(ssfreq_bon, k)))

library("mvtnorm")
rho <- 1 / (1 + alloc)
corr <- matrix(rho, k, k) + diag(1 - rho, k)
quan <- qmvnorm(0.95, mean=rep(0, k), corr=corr)$quantile
ssfreq_dun <- ((quan + qnorm(power)) / (sqrt(nu) * deltastar))^2 * (1 + 1/alloc)
ceiling(c(sqrt(k) * ssfreq_dun, rep(ssfreq_dun, k)))

library("MAMS")
pstar <- pnorm(deltastar / sqrt(2 * 1/nu))
mams(K=k, J=1, r=1, r0=alloc, p=pstar, p0=0.5)

ssbayes(k=2, nu=1, q0=c(0, 0, 0), eta=0.95, zeta=0.90, deltastar=0.5, prec="known",
        crit="2")

ssfreq_unadj <- ((qnorm(1 - alpha) + qnorm(power)) / (sqrt(nu) * deltastar))^2 *
  (1 + 1/sqrt(k))
ceiling(c(sqrt(k) * ssfreq_unadj, rep(ssfreq_unadj, k)))

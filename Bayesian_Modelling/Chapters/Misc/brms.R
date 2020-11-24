# https://www.barelysignificant.com/slides/RGUG2019#51

library('brms')
library(lme4)
library(dplyr)
library(ggplot2)
data(sleepstudy)
head(sleepstudy, 10)

sleepstudy %>%
  ggplot(aes(x = Days, y = Reaction) ) +
  geom_smooth(method = "lm", colour = "black") +
  geom_point() +
  facet_wrap(~Subject, nrow = 2) +
  theme_bw(base_size = 10)

get_prior(Reaction ~ Days + (1 + Days|Subject), sleepstudy)

prior1 <- c(
  prior(normal(200, 10), class = Intercept),
  prior(normal(0, 10), class = b, coef = Days),
  prior(cauchy(0, 10), class = sigma),
  prior(lkj(2), class = cor)
)

mod1 <- brm(
  Reaction ~ Days + (1 + Days | Subject),
  data = sleepstudy,
  family = gaussian(),
  prior = prior1,
  warmup = 2000, iter = 1e4,
  chains = 2
)

posterior_summary(mod1, pars = c("^b_", "^sd_", "sigma"), probs = c(0.025, 0.975) )

ex1 <- brm(Reaction ~ Days + (1 + Days | Subject))
ex2 <- brm(mvbind(Reaction, Memory) ~ Days + (1 + Days | Subject))
ex3 <- brm(mvbind(Reaction, Memory) ~ Days + (1 + Days | Subject))
ex4 <- brm(mvbind(Reaction, Memory) ~ 1 + Days + (1 + Days | Subject))
ex5 <- brm(mvbind(Reaction, Memory) ~ 0 + Days + (1 + Days | Subject))

ex6 <- brm(mvbind(Reaction, Memory) ~ Days + (1 | Subject))
ex7 <- brm(mvbind(Reaction, Memory) ~ Days + (Days | Subject))

ex8 <- brm(mvbind(Reaction, Memory) ~ Days + (1 + Days || Subject))

ex9 <- brm(Reaction ~ 1 + Days + (1 + Days | Subject), family = lognormal() )


library(tidyverse)
(data <-
    read.csv(
      "https://raw.githubusercontent.com/lnalborczyk/lnalborczyk.github.io/master/post/absence2.csv",
      stringsAsFactors = FALSE) )

data %>% sample_frac %>% head(10)


prior2 <- c(
  prior(normal(0, 10), class = Intercept, coef = ""),
  prior(cauchy(0, 10), class = sd),
  prior(normal(0, 10), class = b),
  prior(lkj(2), class = cor)
)

mod2 <- brm(
  presence | trials(total) ~ 1 + reminder + (1 + reminder|researcher), 
  family = binomial(link = "logit"),
  prior = prior2,
  data = data,
  sample_prior = TRUE,
  warmup = 2000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = 0.95)
)

mod2 %>%
  plot(
    combo = c("hist", "trace"), widths = c(1, 1.5),
    theme = theme_bw(base_size = 16)
  )

posterior_summary(mod2, pars = c("^b_", "^sd_"), probs = c(0.025, 0.975) )

# retrieving the intercept
a <- fixef(mod2)[1]
# transforming it back to the probability scale (equivalent to plogis(a))
exp(a) / (1 + exp(a) )

fixef(mod2)[2, c(1, 3, 4)] %>% exp

data %>%
  group_by(researcher, total) %>%
  data_grid(reminder = seq_range(reminder, n = 1e2) ) %>%
  add_fitted_samples(mod2, newdata = ., n = 100, scale = "linear") %>%
  mutate(estimate = plogis(estimate) ) %>%
  ggplot(aes(x = reminder, y = estimate, group = .iteration) ) +
  geom_hline(yintercept = 0.5, lty = 2) +
  geom_line(aes(y = estimate, group = .iteration), size = 0.5, alpha = 0.1) +
  facet_wrap(~researcher, nrow = 2) +
  theme_bw(base_size = 20) + labs(x = "Reminder", y = "Estimate")\


(hyp1 <- hypothesis(mod2, "reminder = 0") )

1 / hyp1$hypothesis$Evid.Ratio

plot(hyp1, theme = theme_bw(base_size = 20) )

data.frame(prior = hyp1$prior_samples$H1, posterior = hyp1$samples$H1) %>%
  gather(type, value) %>%
  ggplot(aes(x = value) ) +
  geom_histogram(bins = 50, alpha = 0.8) +
  geom_vline(xintercept = 0, lty = 2, size = 1) +
  facet_wrap(~type, scales = "free") +
  xlab(expression(beta[reminder]) ) +
  theme_bw(base_size = 20)

prior3 <- c(
  prior(normal(0, 10), class = Intercept, coef = ""),
  prior(cauchy(0, 10), class = sd),
  prior(lkj(2), class = cor) )
mod2 <- brm(presence | trials(total) ~ 1 + reminder + (1 + reminder|researcher), 
            family = binomial(link = "logit"),
            prior = prior2,
            data = data,
            # this line is important for bridgesampling
            save_all_pars = TRUE,
            warmup = 2000, iter = 1e4,
            cores = parallel::detectCores(),
            control = list(adapt_delta = 0.95) )
mod3 <- brm(presence | trials(total) ~ 1 + (1 + reminder|researcher), 
            family = binomial(link = "logit"),
            prior = prior3,
            data = data,
            save_all_pars = TRUE,
            warmup = 2000, iter = 1e4,
            cores = parallel::detectCores(),
            control = list(adapt_delta = 0.95) )

bayes_factor(mod2, mod3)

data %>%
  ggplot(aes(x = presence / total) ) +
  geom_density(fill = "grey20") +
  theme_bw(base_size = 20)

pp_check(mod2, nsamples = 1e2) + theme_bw(base_size = 20)

pp_check(mod2, nsamples = 1e3, type = "stat_2d") + theme_bw(base_size = 20)

mod2 <- brm(
  presence | trials(total) ~ 1 + reminder + (1 + reminder|researcher), 
  family = binomial(link = "logit"),
  prior = prior2,
  data = data,
  warmup = 2000, iter = 1e4,
  cores = parallel::detectCores(), # using all availables cores
  control = list(adapt_delta = 0.95) # adjusting the delta step size
)








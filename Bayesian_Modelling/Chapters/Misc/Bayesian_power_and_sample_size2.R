# https://solomonkurz.netlify.app/post/bayesian-power-analysis-part-i/
  
library(tidyverse)
library(brms)
library(broom)

theme_set(theme_dark() +
            theme(legend.position = "none",
                  panel.grid = element_blank(),
                  panel.background = element_rect(fill = "grey20"),
                  plot.background = element_rect(fill = "grey20"),
                  text = element_text(color = "white"),
                  axis.text = element_text(color = "white")))

# define the means
mu_c <- 0
mu_t <- 0.5

# set up the data

tibble(x = seq(from = -4, to = 5, by = .01)) %>%
  mutate(c = dnorm(x, mean = mu_c, sd = 1),
         t = dnorm(x, mean = mu_t, sd = 1)) %>% 
  
  # plot
  ggplot(aes(x = x, ymin = 0)) +
  geom_ribbon(aes(ymax = c),
              size = 0, alpha = 0.25, fill = "red") +
  geom_ribbon(aes(ymax = t),
              size = 0, alpha = 1/3, fill = "blue2") +
  geom_text(data = tibble(x = c(-.5, 1),
                          y = .385,
                          label = c("control", "treatment"),
                          hjust = 1:0),
            aes(y = y, label = label, color = label, hjust = hjust),
            size = 5, show.legend = F) +
  scale_x_continuous(NULL, breaks = -4:5) +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_color_manual(values = c("red", "blue2"))

n <- 50

set.seed(1)

d <-
  tibble(group     = rep(c("control", "treatment"), each = n)) %>% 
  mutate(treatment = ifelse(group == "control", 0, 1),
         y         = ifelse(group == "control", 
                            rnorm(n, mean = mu_c, sd = 1),
                            rnorm(n, mean = mu_t, sd = 1)))

glimpse(d)

get_prior(data = d,
          family = gaussian,
          y ~ 0 + Intercept + treatment)

fit <-
  brm(data = d,
      family = gaussian,
      y ~ 0 + Intercept + treatment,
      prior = c(prior(normal(0, 2), class = b),
                prior(student_t(3, 1, 1), class = sigma)),
      seed = 1)

plot(fit)

tidy(fit, prob = .95)

# set a new seed
set.seed(2)

# simulate new data based on that new seed
d <- tibble(group     = rep(c("control", "treatment"), each = n)) %>% 
     mutate(treatment = ifelse(group == "control", 0, 1),
      y = ifelse(group == "control", 
                            rnorm(n, mean = mu_c, sd = 1),
                            rnorm(n, mean = mu_t, sd = 1)))

updated_fit <-update(fit, newdata = d,seed = 2)

tidy(updated_fit, prob = .95)


sim_d <- function(seed, n) {
  
  mu_t <- .5
  mu_c <- 0
  
  set.seed(seed)
  
  tibble(group = rep(c("control", "treatment"), n/2)) %>% 
    mutate(treatment = ifelse(group == "control", 0, 1),
           y  = ifelse(group == "control", 
                rnorm(n, mean = mu_c, sd = 1),
                rnorm(n, mean = mu_t, sd = 1)))
}

# how many simulations would you like?
n_sim <- 100
n <- 80
# this will help us track time
t1 <- Sys.time()

# here's the main event!
sims <- tibble(seed = 1:n_sim) %>% 
  mutate(d = map(seed, sim_d, n = n)) %>% 
  mutate(fit = map2(d, seed, ~ update(fit, newdata = .x, seed = .y)))

t2 <- Sys.time()
t2 - t1

sims <- sims %>% 
  mutate(treatment = map(fit, tidy, prob = .95)) %>% 
  unnest(treatment) %>% 
  filter(term == "b_treatment") %>% 
  mutate(powered = case_when(lower > 0 ~ 1, TRUE ~ 0))

sims %>% 
  ggplot(aes(x = seed, y = estimate, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_pointrange(size = .3, aes(color = factor(powered))) +
  labs(y = "Estimate (Credibility Interval)",
       x = "Index") +
  scale_color_manual(values = c("dodgerblue1", "#F1C40F"))

n_sim <- 100
n <- 140
sims140 <- tibble(seed = 1:n_sim) %>% 
  mutate(d = map(seed, sim_d, n = n)) %>% 
  mutate(fit = map2(d, seed, ~ update(fit, newdata = .x, seed = .y))) %>%
  mutate(treatment = map(fit, tidy, prob = .95)) %>% 
  unnest(treatment) %>% 
  filter(term == "b_treatment") %>% 
  mutate(powered = case_when(lower > 0 ~ 1, TRUE ~ 0)) %>%
  mutate(sample = n)

n <- 120
sims120 <- tibble(seed = 1:n_sim) %>% 
  mutate(d = map(seed, sim_d, n = n)) %>% 
  mutate(fit = map2(d, seed, ~ update(fit, newdata = .x, seed = .y))) %>%
  mutate(treatment = map(fit, tidy, prob = .95)) %>% 
  unnest(treatment) %>% 
  filter(term == "b_treatment") %>% 
  mutate(powered = case_when(lower > 0 ~ 1, TRUE ~ 0)) %>%
  mutate(sample = n)

n <- 100
sims100 <- tibble(seed = 1:n_sim) %>% 
  mutate(d = map(seed, sim_d, n = n)) %>% 
  mutate(fit = map2(d, seed, ~ update(fit, newdata = .x, seed = .y))) %>%
  mutate(treatment = map(fit, tidy, prob = .95)) %>% 
  unnest(treatment) %>% 
  filter(term == "b_treatment") %>% 
  mutate(powered = case_when(lower > 0 ~ 1, TRUE ~ 0)) %>%
  mutate(sample = n)

n <- 80
sims80 <- tibble(seed = 1:n_sim) %>% 
  mutate(d = map(seed, sim_d, n = n)) %>% 
  mutate(fit = map2(d, seed, ~ update(fit, newdata = .x, seed = .y))) %>%
  mutate(treatment = map(fit, tidy, prob = .95)) %>% 
  unnest(treatment) %>% 
  filter(term == "b_treatment") %>% 
  mutate(powered = case_when(lower > 0 ~ 1, TRUE ~ 0)) %>%
  mutate(sample = n)

n <- 60
sims60 <- tibble(seed = 1:n_sim) %>% 
  mutate(d = map(seed, sim_d, n = n)) %>% 
  mutate(fit = map2(d, seed, ~ update(fit, newdata = .x, seed = .y))) %>%
  mutate(treatment = map(fit, tidy, prob = .95)) %>% 
  unnest(treatment) %>% 
  filter(term == "b_treatment") %>% 
  mutate(powered = case_when(lower > 0 ~ 1, TRUE ~ 0)) %>%
  mutate(sample = n)

n <- 40
sims40 <- tibble(seed = 1:n_sim) %>% 
  mutate(d = map(seed, sim_d, n = n)) %>% 
  mutate(fit = map2(d, seed, ~ update(fit, newdata = .x, seed = .y))) %>%
  mutate(treatment = map(fit, tidy, prob = .95)) %>% 
  unnest(treatment) %>% 
  filter(term == "b_treatment") %>% 
  mutate(powered = case_when(lower > 0 ~ 1, TRUE ~ 0)) %>%
  mutate(sample = n)

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
  geom_pointrange(size = .2, aes(color = factor(powered))) +
  labs(y = "Estimate (Credibility Interval)",
       x = "Index") +
  scale_color_manual(values = c("dodgerblue1", "#F1C40F")) +
  facet_wrap(~sample)

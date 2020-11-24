library(rethinking)
library(brms)
library(tidyverse)
library(ggrepel)
data(WaffleDivorce)
d <- WaffleDivorce

d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)


priors <- c(prior(normal(0, 0.2), class = Intercept),
            prior(normal(0, 0.5), class = b),
            prior(exponential(1), class = sigma))

b5.1 <- brm(data = d, family = gaussian,
      D ~ 1 + A,
      prior = priors,
      iter = 2000, warmup = 500, chains = 4, cores = 4,
      seed = 5)

print(b5.1)

# define the range of `MedianAgeMarriage_s` values we'd like to feed into `fitted()`
A_seq <- tibble(A = seq(from = -3, to = 3.2, length.out = 30))

f <- fitted(b5.1, newdata = A_seq) %>%
  as_tibble() %>%
  # tack the `nd` data onto the `fitted()` results
  bind_cols(A_seq)

ggplot(data = f, 
       aes(x = A, y = Estimate)) +
  geom_smooth(aes(ymin = Q2.5, ymax = Q97.5),
              stat = "identity",
              fill = "blue", color = "blue", alpha = 1/5, size = 1/4) +
  geom_point(data = d, 
             aes(y = D), 
             size = 2.5, color = "blue", alpha=0.4) +
  ylab("Divorce") +
  coord_cartesian(xlim = range(d$A), 
                  ylim = range(d$D)) +
  theme_bw() +
  theme(panel.grid = element_blank())   


priors <- c(prior(normal(0, 0.2), class = Intercept),
            prior(normal(0, 0.5), class = b),
            prior(exponential(1), class = sigma))

b5.2 <- brm(data = d, family = gaussian,
      D ~ 1 + M,
      prior=priors,
      iter = 2000, warmup = 500, chains = 4, cores = 4,
      seed = 5)

A_seq <- tibble(M = seq(from = -2.5, to = 3.5, length.out = 30))

f <- 
  fitted(b5.2, newdata = A_seq) %>%
  as_tibble() %>%
  bind_cols(A_seq)

ggplot(data = f, 
       aes(x = M, y = Estimate)) +
  geom_smooth(aes(ymin = Q2.5, ymax = Q97.5),
              stat = "identity",
              fill = "blue", color = "blue", alpha = 1/5, size = 1/4) +
  geom_point(data = d, 
             aes(y = D), 
             size = 2.5, color = "blue", alpha=0.4) +
  coord_cartesian(xlim = range(d$M), 
                  ylim = range(d$D)) +
  ylab("Divorce") +
  theme_bw() +
  theme(panel.grid = element_blank())  

b5.3 <- brm(data = d, family = gaussian,
      D ~ 1 + M + A,
      prior=priors,
      iter = 2000, warmup = 500, chains = 4, cores = 4,
      seed = 5)

stanplot(b5.3)


library(tidybayes)


post <- posterior_samples(b5.3)

post %>% 
  select(-lp__) %>% 
  gather() %>% 
  
  ggplot(aes(x = value, y = reorder(key, value))) +  # note how we used `reorder()` to arrange the coefficients
  geom_vline(xintercept = 0, color = "firebrick4", alpha = 3/10) +
  stat_pointintervalh(point_interval = mode_hdi, .width = .95, 
                      size = 2, color = "firebrick4") +
  labs(title = "My tidybayes-based coefficient plot",
       x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid   = element_blank(),
        panel.grid.major.y = element_line(color = alpha("firebrick4", 2), linetype = 3),
        axis.text.y  = element_text(hjust = 0, size=11),
        axis.ticks.y = element_blank())

b5.4 <- brm(data = d, family = gaussian,
      M ~ 1 + A,
      prior = priors,
      iter = 2000, warmup = 500, chains = 4, cores = 4,
      seed = 5)


f <- fitted(b5.4) %>%
    as_tibble() %>%
    bind_cols(d)

head(f)

f %>% 
  
  ggplot(aes(x = A, y = M)) +
  geom_point(size = 3.5, shape = 1, color = "blue") +
  geom_segment(aes(xend = A, yend = Estimate), 
               size = 1/4) +
  geom_line(aes(y = Estimate), 
            color = "firebrick4") +
  coord_cartesian(ylim = range(d$M)) +
  theme_bw() +
  theme(panel.grid = element_blank())  

r <- residuals(b5.4) %>%
  # to use this in ggplot2, we need to make it a tibble or data frame
    as_tibble() %>% 
    bind_cols(d)

# for the annotation at the top
text <- tibble(Estimate = c(- 0.5, 0.5),
         Divorce = 14.1,
         label = c("slower", "faster"))

r %>% 
  ggplot(aes(x = Estimate, y = Divorce)) +
  stat_smooth(method = "lm", fullrange = T,
              color = "firebrick4", fill = "firebrick4", 
              alpha = 1/5, size = 1/2) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 3, color = "firebrick4", alpha = 2/3) +
  geom_text(data = text,
            aes(label = label)) +
  scale_x_continuous("Marriage rate residuals", limits = c(-2, 2)) +
  coord_cartesian(xlim = range(r$Estimate),
                  ylim = c(6, 14.1)) +
  theme_bw() +
  theme(panel.grid = element_blank())


b5.4b <- brm(data = d, family = gaussian,
      A ~ 1 + M,
      prior = priors,
      iter = 2000, warmup = 500, chains = 4, cores = 4,
      seed = 5)


text <-
  tibble(Estimate = c(- 0.7, 0.5),
         Divorce  = 14.1,
         label    = c("younger", "older"))

residuals(b5.4b) %>%
  as_tibble() %>%
  bind_cols(d) %>% 
  
  ggplot(aes(x = Estimate, y = Divorce)) +
  stat_smooth(method = "lm", fullrange = T,
              color = "firebrick4", fill = "firebrick4", 
              alpha = 1/5, size = 1/2) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 3, color = "firebrick4", alpha = 2/3) +
  geom_text(data = text,
            aes(label = label)) +
  scale_x_continuous("Age of marriage residuals", limits = c(-2, 3)) +
  coord_cartesian(xlim = range(r$Estimate),
                  ylim = c(6, 14.1)) +
  theme_bw() +
  theme(panel.grid = element_blank())  


library(rethinking)
library(brms)


a <-  1.4
sigma <-  1.5
n_ponds <- 60

set.seed(12)

  dsim <- 
    tibble(pond   = 1:n_ponds,
           ni     = rep(c(5, 10, 25, 35), each = n_ponds / 4) %>% as.integer(),
           true_a = rnorm(n = n_ponds, mean = a, sd = sigma))
  set.seed(12)
  (
    dsim <-
      dsim %>%
      mutate(si = rbinom(n = n(), prob = inv_logit_scaled(true_a), size = ni))
  )
  
  (
    dsim <-
      dsim %>%
      mutate(p_nopool = si / ni)
  )
  
priors <- c(prior(normal(0, 1), class = Intercept),
          prior(cauchy(0, 1), class = sd))

b12.3 <-brm(data = dsim, family = binomial,
        si | trials(ni) ~ 1 + (1 | pond),
        iter = 10000, warmup = 1000, chains = 1, cores = 1,
        seed = 12)  
  
  
data(chimpanzees)
d <- chimpanzees
d$treatment <- d$prosoc_left+2*d$condition

priors <- prior(normal(0, 10), class = Intercept)

b10.1 <-  brm(data = d, family = binomial,
          pulled_left | trials(1) ~ 1,
          prior = priors,
          seed = 10)

priors <- c(prior(normal(0, 10), class = Intercept),
          prior(normal(0, 10), class = b))

b10.2 <-  brm(data = d, family = binomial,
          pulled_left | trials(1) ~ 1 + prosoc_left,
          prior = priors,
          seed = 10)

b10.3 <- update(b10.2,
         newdata = d,
         formula = pulled_left | trials(1) ~ 1 + treatment)


priors <- c(prior(normal(0, 10), class = Intercept),
            prior(normal(0, 10), class = b),
            prior(cauchy(0, 1), class = sd))


d_aggregated <-
  d %>%
  select(-recipient, -block, -trial, -chose_prosoc) %>%
  group_by(actor, treatment) %>%
  summarise(x = sum(pulled_left))

priors = c(prior(normal(0, 10), class = Intercept),
          prior(normal(0, 10), class = b))

b10.5 <-
  brm(data = d_aggregated, family = binomial,
      x | trials(18) ~ 1 + treatment,
      prior = priors, 
      iter = 2500, warmup = 500, cores = 2, chains = 2, 
      seed = 10)

d_aggregated %>%
  filter(actor %in% c(1, 2))

b12.4 <- 
  brm(data = d, family = binomial,
      pulled_left | trials(1) ~ 1 + treatment + (1 | actor),
      prior=priors,
      # I'm using 4 cores, instead of the `cores=3` in McElreath's code
      iter = 5000, warmup = 1000, chains = 4, cores = 4,  
      control = list(adapt_delta = 0.95),
      seed = 12)

b12.5 <- 
  update(b12.4,
         newdata = d,
         formula = pulled_left | trials(1) ~ 1 + treatment + (1 | actor) + (1 | block),
         iter = 6000, warmup = 1000, cores = 4, chains = 4, 
         control = list(adapt_delta = 0.99),
         seed = 12)

library(rethinking)
library(dplyr)
library(wesanderson)

# simulate career choices among 500 individuals
n      <- 500           # number of individuals
income <- 1:3           # expected income of each career
score  <- 0.5 * income  # scores for each career, based on income

# next line converts scores to probabilities
p <- softmax(score[1], score[2], score[3])

# now simulate choice
# outcome career holds event type values, not counts
career <- rep(NA, n)  # empty vector of choices for each individual

set.seed(10)
# sample chosen career for each individual
for(i in 1:n) career[i] <- sample(1:3, size = 1, prob = p)

tibble(career = career) %>%
  ggplot(aes(x = career)) +
  geom_bar(size = 0, fill = wes_palette("Moonrise2")[2])

tibble(career) %>% 
  count(career) %>% 
  mutate(percent     = (100 * n / sum(n)),
         probability =        n / sum(n))

tibble(income = 1:3) %>% 
  mutate(score = 0.5 * income) %>% 
  mutate(p = exp(score) / sum(exp(score)))

tibble(score = 104 + c(0.5, 1, 1.5)) %>% 
  mutate(p = exp(score) / sum(exp(score)))

tibble(income = 1:3,
       score  = c(0, 0.5, 1)) %>% 
  mutate(p = exp(score) / sum(exp(score)))

library(brms)

get_prior(data = list(career = career), 
          family = categorical(link = logit),
          career ~ 1)

priors <- c(prior(normal(0, 5), class = Intercept, dpar = mu2),
            prior(normal(0, 5), class = Intercept, dpar = mu3))
  
b10.16 <-
  brm(data = list(career = career), 
      family = categorical(link = logit),
      prior = priors,
      career ~ 1,
      iter = 2500, warmup = 500, cores = 2, chains = 2,
      seed = 10,
      file = "fits/b10.16")  

print(b10.16)

fixef(b10.16) %>% inv_logit_scaled()

fitted(b10.16) %>% str()

fitted(b10.16)[1, , ] %>% 
  round(digits = 2)

fitted(b10.16)[1, , ] %>% 
  t() %>% 
  round(digits = 2)

tibble(career) %>% 
  count(career) %>% 
  mutate(percent     = (100 * n / sum(n)),
         probability =        n / sum(n))

# this helps us set our custom color scheme
color_scheme_set(wes_palette("Moonrise2")[c(1, 3, 2, 2, 2, 3)])

pp_check(b10.16, type = "hist", binwidth = 1) +
  theme(legend.position = c(.91, .125),
        legend.key = element_rect(color = "transparent"))

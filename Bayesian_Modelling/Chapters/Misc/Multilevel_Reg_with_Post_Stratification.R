# https://cran.r-project.org/web/packages/rstanarm/vignettes/mrp.html
# https://github.com/stan-dev/rstanarm/blob/master/vignettes/mrp.Rmd

library(rstanarm)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default())
options(mc.cores = 4) 
library(dplyr)
library(tidyr)

mrp_sim <- simulate_mrp_data(n=1200)
str(mrp_sim)

sample <- mrp_sim[["sample"]]
rbind(head(sample), tail(sample))

poststrat <- mrp_sim[["poststrat"]]
rbind(head(poststrat), tail(poststrat))

true_popn <- mrp_sim[["true_popn"]]
rbind(head(true_popn), tail(true_popn))

sample$state <- factor(sample$state, levels=1:50)
sample$state <- with(sample, factor(state, levels=order(table(state))))
true_popn$state <- factor(true_popn$state,levels = levels(sample$state))
poststrat$state <- factor(poststrat$state,levels = levels(sample$state))

fit <- stan_glmer(
  cat_pref ~ factor(male) + factor(male) * factor(age) + 
    (1 | state) + (1 | age) + (1 | eth) + (1 | income),
  family = binomial(link = "logit"),
  data = sample
)


posterior_prob <- posterior_linpred(fit, transform = TRUE, newdata = poststrat)
poststrat_prob <- posterior_prob %*% poststrat$N / sum(poststrat$N)
model_popn_pref <- c(mean = mean(poststrat_prob), sd = sd(poststrat_prob))
round(model_popn_pref, 3)


sample_popn_pref <- mean(sample$cat_pref)
round(sample_popn_pref, 3)

compare2 <- compare2 +
  geom_hline(yintercept = model_popn_pref[1], colour = '#2ca25f', size = 1) +
  geom_text(aes(x = 5.2, y = model_popn_pref[1] + .025), label = "MRP", colour = '#2ca25f')
bayesplot_grid(compare, compare2, 
               grid_args = list(nrow = 1, widths = c(8, 1)))

true_popn_pref <- sum(true_popn$cat_pref * poststrat$N) / sum(poststrat$N)
round(true_popn_pref, 3)

state_df <- data.frame(
  State = 1:50,
  model_state_sd = rep(-1, 50),
  model_state_pref = rep(-1, 50),
  sample_state_pref = rep(-1, 50),
  true_state_pref = rep(-1, 50),
  N = rep(-1, 50)
)

for(i in 1:length(levels(as.factor(poststrat$state)))) {
  poststrat_state <- poststrat[poststrat$state == i, ]
  posterior_prob_state <- posterior_linpred(
    fit,
    transform = TRUE,
    draws = 1000,
    newdata = as.data.frame(poststrat_state)
  )
  poststrat_prob_state <- (posterior_prob_state %*% poststrat_state$N) / sum(poststrat_state$N)
  #This is the estimate for popn in state:
  state_df$model_state_pref[i] <- round(mean(poststrat_prob_state), 4)
  state_df$model_state_sd[i] <- round(sd(poststrat_prob_state), 4)
  #This is the estimate for sample
  state_df$sample_state_pref[i] <- round(mean(sample$cat_pref[sample$state == i]), 4)
  #And what is the actual popn?
  state_df$true_state_pref[i] <-
    round(sum(true_popn$cat_pref[true_popn$state == i] * poststrat_state$N) /
            sum(poststrat_state$N), digits = 4)
  state_df$N[i] <- length(sample$cat_pref[sample$state == i])
}

state_df[c(1,3:6)]

state_df$State <- factor(state_df$State, levels = levels(sample$state))

round(100 * c(
  mean = mean(abs(state_df$sample_state_pref-state_df$true_state_pref), na.rm = TRUE),
  max = max(abs(state_df$sample_state_pref-state_df$true_state_pref), na.rm = TRUE)
))

round(100 * c(
  mean = mean(abs(state_df$model_state_pref-state_df$true_state_pref)),
  max = max(abs(state_df$model_state_pref-state_df$true_state_pref))
))



# not evaluated to avoid dependency on tidyverse
sample_alt <- sample %>%
  group_by(male, age, income, state, eth) %>%
  summarise(N_cat_pref = sum(cat_pref), N = n()) %>%
  ungroup()

fit2 <- stan_glmer(
  cbind(N_cat_pref, N - N_cat_pref) ~ factor(male) + factor(male) * factor(age) + 
    (1 | state) + (1 | age) + (1 | eth) + (1 | income),
  family = binomial("logit"),
  data = sample_alt,
  refresh = 0
)

posterior_prob_alt <- posterior_linpred(fit2, transform = TRUE, newdata = poststrat)
poststrat_prob_alt <- posterior_prob_alt %*% poststrat$N / sum(poststrat$N)
model_popn_pref_alt <- c(mean = mean(poststrat_prob_alt), sd = sd(poststrat_prob_alt))
round(model_popn_pref_alt, 3)




############################################################################################################

simulate_mrp_data <- function(n) {
  J <- c(2, 3, 7, 3, 50) # male or not, eth, age, income level, state
  poststrat <- as.data.frame(array(NA, c(prod(J), length(J)+1))) # Columns of post-strat matrix, plus one for size
  colnames(poststrat) <- c("male", "eth", "age","income", "state",'N')
  count <- 0
  for (i1 in 1:J[1]){
    for (i2 in 1:J[2]){
      for (i3 in 1:J[3]){
        for (i4 in 1:J[4]){
          for (i5 in 1:J[5]){
            count <- count + 1
            # Fill them in so we know what category we are referring to
            poststrat[count, 1:5] <- c(i1-1, i2, i3,i4,i5) 
          }
        }
      }
    }
  }
  # Proportion in each sample in the population
  p_male <- c(0.52, 0.48)
  p_eth <- c(0.5, 0.2, 0.3)
  p_age <- c(0.2,.1,0.2,0.2, 0.10, 0.1, 0.1)
  p_income<-c(.50,.35,.15)
  p_state_tmp<-runif(50,10,20)
  p_state<-p_state_tmp/sum(p_state_tmp)
  poststrat$N<-0
  for (j in 1:prod(J)){
    poststrat$N[j] <- round(250e6 * p_male[poststrat[j,1]+1] * p_eth[poststrat[j,2]] *
                              p_age[poststrat[j,3]]*p_income[poststrat[j,4]]*p_state[poststrat[j,5]]) #Adjust the N to be the number observed in each category in each group
  }
  
  # Now let's adjust for the probability of response
  p_response_baseline <- 0.01
  p_response_male <- c(2, 0.8) / 2.8
  p_response_eth <- c(1, 1.2, 2.5) / 4.7
  p_response_age <- c(1, 0.4, 1, 1.5,  3, 5, 7) / 18.9
  p_response_inc <- c(1, 0.9, 0.8) / 2.7
  p_response_state <- rbeta(50, 1, 1)
  p_response_state <- p_response_state / sum(p_response_state)
  p_response <- rep(NA, prod(J))
  for (j in 1:prod(J)) {
    p_response[j] <-
      p_response_baseline * p_response_male[poststrat[j, 1] + 1] *
      p_response_eth[poststrat[j, 2]] * p_response_age[poststrat[j, 3]] *
      p_response_inc[poststrat[j, 4]] * p_response_state[poststrat[j, 5]]
  }
  people <- sample(prod(J), n, replace = TRUE, prob = poststrat$N * p_response)
  
  ## For respondent i, people[i] is that person's poststrat cell,
  ## some number between 1 and 32
  n_cell <- rep(NA, prod(J))
  for (j in 1:prod(J)) {
    n_cell[j] <- sum(people == j)
  }
  
  coef_male <- c(0,-0.3)
  coef_eth <- c(0, 0.6, 0.9)
  coef_age <- c(0,-0.2,-0.3, 0.4, 0.5, 0.7, 0.8, 0.9)
  coef_income <- c(0,-0.2, 0.6)
  coef_state <- c(0, round(rnorm(49, 0, 1), 1))
  coef_age_male <- t(cbind(c(0, .1, .23, .3, .43, .5, .6),
                           c(0, -.1, -.23, -.5, -.43, -.5, -.6)))
  true_popn <- data.frame(poststrat[, 1:5], cat_pref = rep(NA, prod(J)))
  for (j in 1:prod(J)) {
    true_popn$cat_pref[j] <- plogis(
      coef_male[poststrat[j, 1] + 1] +
        coef_eth[poststrat[j, 2]] + coef_age[poststrat[j, 3]] +
        coef_income[poststrat[j, 4]] + coef_state[poststrat[j, 5]] +
        coef_age_male[poststrat[j, 1] + 1, poststrat[j, 3]]
    )
  }
  
  #male or not, eth, age, income level, state, city
  y <- rbinom(n, 1, true_popn$cat_pref[people])
  male <- poststrat[people, 1]
  eth <- poststrat[people, 2]
  age <- poststrat[people, 3]
  income <- poststrat[people, 4]
  state <- poststrat[people, 5]
  
  sample <- data.frame(cat_pref = y, 
                       male, age, eth, income, state, 
                       id = 1:length(people))
  
  #Make all numeric:
  for (i in 1:ncol(poststrat)) {
    poststrat[, i] <- as.numeric(poststrat[, i])
  }
  for (i in 1:ncol(true_popn)) {
    true_popn[, i] <- as.numeric(true_popn[, i])
  }
  for (i in 1:ncol(sample)) {
    sample[, i] <- as.numeric(sample[, i])
  }
  list(
    sample = sample,
    poststrat = poststrat,
    true_popn = true_popn
  )
}

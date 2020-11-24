# http://www.stat.cmu.edu/~ryurko/post/bayesian-baby-steps-intro/
sample(c("C", "D"), size = 10, replace = TRUE)

library(tidyverse)

juju_games <- read_csv("https://raw.githubusercontent.com/ryurko/nflscrapR-data/master/legacy_data/season_play_by_play/pbp_2017.csv") %>%
  filter(Receiver_ID == "00-0033857",
         PassAttempt == 1) %>%
  select(GameID, PassOutcome) %>%
  group_by(GameID) %>%
  summarise(receptions = length(which(PassOutcome == "Complete")),
            targets = n(),
            catch_rate = receptions / targets) %>%
  mutate(total_receptions = cumsum(receptions),
         total_targets = cumsum(targets),
         total_catch_rate = total_receptions / total_targets,
         index = 1:n(),
         game_index = paste("game_", index, sep = ""),
         game_index = fct_relevel(factor(game_index),
                                  "game_1", "game_2", "game_3",
                                  "game_4", "game_5", "game_6",
                                  "game_7", "game_8", "game_9",
                                  "game_10", "game_11", "game_12",
                                  "game_13"))

p_grid <- seq(0, 1, length.out =50)
prior <- rep(1, 50)

likelihood <- dbinom(x=juju_games$receptions[1],
                    size=juju_games$targets[1],
                    prob = p_grid)

bayes_numerator <- likelihood * prior

posterior <- bayes_numerator/sum(bayes_numerator)

data.frame(p_grid = p_grid, p_posterior = posterior) %>%
  ggplot(aes(x = p_grid, y = p_posterior)) +
  geom_point(size = 3, color = "darkblue") + 
  geom_line(color = "darkblue") +
  # Add a vertical line for JuJu's observed catch rate:
  geom_vline(xintercept = juju_games$catch_rate[1], color = "darkorange",
             linetype = "dashed", size = 3, alpha = .5) +
  # Label!
  labs(x = "Catch rate", y = "Posterior probability",
       title = "Posterior approximation for\nJuJu's catch rate after one game") +
  # Clean it up:
  theme_bw() + theme(axis.text = element_text(size = 10), 
                     title = element_text(size = 10))

game_posteriors <- map_dfc(c(1:nrow(juju_games)),
                           function(x) {
                             p_grid <- seq(from=0,to=1,by=0.05)
                             prior <- rep(1, 21)
                             likelihood <- dbinom(x = juju_games$total_receptions[x],
                                                  size=juju_games$total_targets[x],
                                                  prob=p_grid)
                             bayes_numerator <- likelihood * prior
                             posterior <- bayes_numerator/sum(bayes_numerator)
                             result <- data.frame(posterior)
                             colnames(result) <- paste('game_', x, sep='')
                             return(result)
                           })


# Join these columns with p_grid and column for the prior probability:
 df <- data.frame(p_grid = p_grid, prior = rep(1/21, 21)) %>%
  bind_cols(game_posteriors) %>% 
  # Gather the columns so the data is long, one row for each week and grid value
  gather(key = "game_index", value = "posterior_prob", -p_grid) %>%
  # Relevel the game_index variable:
  mutate(game_index = fct_relevel(factor(game_index),
                                  "prior", "game_1", "game_2", "game_3",
                                  "game_4", "game_5", "game_6",
                                  "game_7", "game_8", "game_9",
                                  "game_10", "game_11", "game_12",
                                  "game_13"))
  # Visualize the posteriors for each game:
ggplot(df, aes(x = p_grid, y = posterior_prob)) + 
  geom_point(size = 2, color = "darkblue") + 
  geom_line(color = "darkblue") +
  facet_wrap(~ game_index) +
  # Add vertical lines for each cumulative observed rate
  geom_vline(data = juju_games, 
             aes(xintercept = total_catch_rate), color = "darkorange",
             linetype = "dashed", size = 1, alpha = .5) +
  geom_text(data = juju_games, size = 3,
            x = .25, y = .3, aes(label = paste("Caught", 
                                               receptions, "of",
                                               targets, sep = " "))) +
  # Label!
  labs(x = "Catch rate", y = "Posterior probability",
       title = "Posterior approximation for JuJu's catch rate after each game") +
  # Clean it up:
  theme_bw() + theme(axis.text.y = element_text(size = 10), 
                     axis.text.x = element_text(size = 6),
                     title = element_text(size = 10)) 

library(nflWAR)

wr_catch_rates <- get_pbp_data(2016) %>%
  add_positions(2016) %>%
  add_model_variables() %>%
  prepare_model_data() %>%
  add_position_tables() %>%
  join_position_statistics() %>%
  pluck("WR_table") %>%
  select(Player_ID_Name, Rec_Perc, Targets) %>%
  filter(Targets >= 25)

wr_catch_rates %>%
  ggplot(aes(x = Rec_Perc)) + 
  geom_density(color = "darkblue") +
  geom_rug(color = "darkblue") + 
  theme_bw() +
  labs(x = "Catch rate", y = "Density", 
       title = "Distribution of WR catch rates in 2016 (minimum 25 targets)")

# Use the fitdistr function from the MASS library:
# Maximum-likelihood fitting of univariate distributions, 
# allowing parameters to be held fixed if desired.
# USAGE fitdistr(x, densfun, start, .)
# x = A numeric vector of length at least one containing only finite values.
# densfun = type of distribution
# start = A named list giving the parameters to be optimized with initial values.
prior_beta_model <- MASS::fitdistr(wr_catch_rates$Rec_Perc, dbeta, 
                                   start = list(shape1 = 10, shape2 = 10))

# Extract the approximated  priors
alpha0 <- prior_beta_model$estimate[1]
beta0 <- prior_beta_model$estimate[2]

# Create a data frame by applying the grid approximation steps to each row
# of juju_games:
new_game_posteriors <- map_dfc(c(1:nrow(juju_games)),
                               function(x) {
                                 p_grid <- seq(from = 0, to = 1, by = .05)
                                 prior <- dbeta(p_grid, alpha0, beta0)
                                 likelihood <- dbinom(x = juju_games$total_receptions[x],
                                                      size = juju_games$total_targets[x],
                                                      prob = p_grid)
                                 bayes_numerator <- likelihood * prior
                                 posterior <- bayes_numerator / sum(bayes_numerator)
                                 # Return this as a data frame:
                                 result <- data.frame(posterior)
                                 colnames(result) <- paste("game_", x, sep = "")
                                 return(result)
                               })

# Join these columns with p_grid and column for the prior probability:
new_df <- data.frame(p_grid = p_grid, 
           prior = dbeta(p_grid, alpha0, beta0) / sum(dbeta(p_grid, alpha0, beta0))) %>%
  bind_cols(new_game_posteriors) %>% 
  # Gather the columns so the data is long, one row for each week and grid value
  gather(key = "game_index", value = "posterior_prob", -p_grid) %>%
  # Relevel the game_index variable:
  mutate(game_index = fct_relevel(factor(game_index),
                                  "prior", "game_1", "game_2", "game_3",
                                  "game_4", "game_5", "game_6",
                                  "game_7", "game_8", "game_9",
                                  "game_10", "game_11", "game_12",
                                  "game_13"))
  # Visualize the posteriors for each game:
ggplot(new_df, aes(x = p_grid, y = posterior_prob)) + 
  geom_point(size = 2, color = "darkblue") + 
  geom_line(color = "darkblue") +
  facet_wrap(~ game_index) +
  # Add vertical lines for each cumulative observed rate
  geom_vline(data = juju_games, 
             aes(xintercept = total_catch_rate), color = "darkorange",
             linetype = "dashed", size = 1, alpha = .5) +
  geom_text(data = juju_games, size = 3,
            x = .25, y = .3, aes(label = paste("Caught", 
                                               receptions, "of",
                                               targets, sep = " "))) +
  # Label!
  labs(x = "Catch rate", y = "Posterior probability",
       title = "Posterior approximation for JuJu's catch rate after each game") +
  # Clean it up:
  theme_bw() + theme(axis.text.y = element_text(size = 10), 
                     axis.text.x = element_text(size = 6),
                     title = element_text(size = 10)) 

kn_post_df <- data.frame(p_density = dbeta(seq(0, 1, length.out = 1000),
                              # Add receptions to prior alpha
                              juju_games$total_receptions[nrow(juju_games)] + alpha0,
                              # Add incompletions to prior beta
                              with(juju_games, 
                              (total_targets[nrow(juju_games)] - 
                              total_receptions[nrow(juju_games)]) + beta0)),
                              p_grid = seq(0, 1, length.out = 1000))

ggplot(kn_post_df, aes(x = p_grid, y = p_density)) +
  geom_line(size=1.5, color = "darkorange") +
  # Label!
  labs(x = "Catch rate", y = "Posterior density",
       title = "Known posterior distribution using beta prior") +
  # Clean it up:
  theme_bw() + theme(axis.text = element_text(size = 10), 
                     title = element_text(size = 10)) 

grid_posterior <- new_game_posteriors %>%
  select(game_13) %>%
  bind_cols(data.frame(p_grid = p_grid))
  
  
ggplot(grid_posterior, aes(x = p_grid, y = game_13)) + 
  geom_point(size = 2, color = "darkblue") + 
  geom_line(color = "darkblue") +
  # Label!
  labs(x = "Catch rate", y = "Posterior probability",
       title = "Grid approximation") +
  # Clean it up:
  theme_bw() + theme(axis.text = element_text(size = 10), 
                     title = element_text(size = 10)) 

# Install cowplot if you don't have it!
install.packages("cowplot")
cowplot::plot_grid(known_posterior, grid_posterior)

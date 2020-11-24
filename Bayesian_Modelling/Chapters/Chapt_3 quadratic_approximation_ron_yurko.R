library(tidyverse)
library(nflWAR)
library(rethinking)
library(mvtnorm)
team_summary_df <- get_season_summary(2009:2017)

score_diff_hist <- team_summary_df %>%
  ggplot(aes(x = Total_Score_Diff)) +
  geom_histogram(aes(y = ..density..), bins = 25, colour="white", fill="grey") + theme_bw() + 
  scale_x_continuous(limits = c(-300, 300)) +  
  geom_density(alpha=.1, fill="blue") + #geom_line(size=1, color='steelblue', aes(x=neuron, y=lambda))
  geom_rug(aes(x = Total_Score_Diff, y = 0), position = position_jitter(height = 0)) +
  # Always label:
  labs(x = "Score differential", y = "Density",
       title = "Distribution of individual team score differential in each season from 2009-17 ",
       caption = "Data accessed with nflscrapR") +
  # Just some cosmetic settings: (really should do a custom theme at the 
  # beginning instead of repeating this)
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14))

score_diff_hist

laplace_approx <- function(log_posterior_fn, init_points, n_samples, ...) {
  fit <- optim(init_points, log_posterior_fn, 
               method = "BFGS",
               control = list(fnscale = -1),
               hessian = TRUE,
               ...)
  
  param_mean <- fit$par
  param_cov_mat <- solve(-fit$hessian)
  
  mvtnorm::rmvnorm(n_samples, param_mean, param_cov_mat) %>%
    data.frame()
}

score_diff_model <- function(params, score_diff_values) {
  sum(dnorm(score_diff_values, params["mu"], params["sigma"], log = TRUE)) +
    dnorm(params["mu"], 0, 100, log = TRUE) + 
    dunif(params["sigma"], 0, 200, log = TRUE)
}

init_samples <- laplace_approx(log_posterior_fn = score_diff_model,
                               init_points = c(mu = 200, sigma = 10),
                               n_samples = 10000,
                               # Specify the data for the optimization!
                               score_diff_values = team_summary_df$Total_Score_Diff)

library(cowplot)

# First the joint distribution plot with points for the values:
joint_plot <- ggplot(init_samples, aes(x = mu, y = sigma)) + 
  geom_point(alpha = .3, color = "darkblue")  +
  labs(title = "Posterior distribution for model parameters with marginals") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16))

# Marginal density for mu along x-axis using the cowplot function axis_canvas:
mu_dens <- axis_canvas(joint_plot, axis = "x") +
  geom_density(data = init_samples, aes(x = mu), fill = "darkblue",
               alpha = 0.5, size = .2)

# Same thing for sigma but along y and with coord_flip = TRUE to make it vertical:
sigma_dens <- axis_canvas(joint_plot, axis = "y", coord_flip = TRUE) +
  geom_density(data = init_samples, aes(x = sigma), fill = "darkblue",
               alpha = 0.5, size = .2) +
  coord_flip()

# Now generate by adding these objects to the main plot:
# Need grid:
# install.packages("grid")
joint_plot_1 <- insert_xaxis_grob(joint_plot, mu_dens, 
                                  grid::unit(.2, "null"), position = "top")
joint_plot_2 <- insert_yaxis_grob(joint_plot_1, sigma_dens, 
                                  grid::unit(.2, "null"), position = "right")
ggdraw(joint_plot_2)


# R is vectorized so can just give it the vector of values for the Gaussian
# distribution parameters to create the simulated score-diff distribution:
sim_score_diff_data <- data.frame(Total_Score_Diff = 
                                    rnorm(10000, mean = init_samples$mu,
                                          sd = init_samples$sigma))

# Create the histogram from before but now add this posterior density on top:
score_diff_hist + geom_density(data = sim_score_diff_data,
                               aes(x = Total_Score_Diff), 
                               fill = NA, color = "red")

# Load in the four datasets summarising team performance for passing and rushing
# both from offensive and defensive perspectives, selecting only the necessary
# columns and renaming them:

team_passing_off <- read_csv("https://raw.githubusercontent.com/ryurko/nflscrapR-data/master/legacy_data/season_team_stats/team_season_passing_df.csv") %>%
  dplyr::select(Season, Team, EPA_per_Att) %>%
  rename(pass_epa_att_off = EPA_per_Att) %>%
  filter(!is.na(Team))

team_passing_def <- read_csv("https://raw.githubusercontent.com/ryurko/nflscrapR-data/master/legacy_data/season_team_stats/team_def_season_passing_df.csv") %>%
  dplyr::select(Season, Team, EPA_per_Att) %>%
  rename(pass_epa_att_def = EPA_per_Att) %>%
  filter(!is.na(Team))

team_rushing_off <- read_csv("https://raw.githubusercontent.com/ryurko/nflscrapR-data/master/legacy_data/season_team_stats/team_season_rushing_df.csv") %>%
  dplyr::select(Season, Team, EPA_per_Car) %>%
  rename(rush_epa_att_off = EPA_per_Car) %>%
  filter(!is.na(Team))

team_rushing_def <- read_csv("https://raw.githubusercontent.com/ryurko/nflscrapR-data/master/legacy_data/season_team_stats/team_def_season_rushing_df.csv") %>%
  dplyr::select(Season, Team, EPA_per_Car) %>%
  rename(rush_epa_att_def = EPA_per_Car) %>%
  filter(!is.na(Team))

# Join the data to the team_summary_df and calculate the differential columns:

team_summary_df <- team_summary_df %>%
  inner_join(team_passing_off, by = c("Team", "Season")) %>%
  inner_join(team_passing_def, by = c("Team", "Season")) %>%
  inner_join(team_rushing_off, by = c("Team", "Season")) %>%
  inner_join(team_rushing_def, by = c("Team", "Season")) %>%
  mutate(pass_epa_diff = pass_epa_att_off - pass_epa_att_def,
         rush_epa_diff = rush_epa_att_off - rush_epa_att_def)

epa_pass_plot <- ggplot(team_summary_df,
                        aes(x = pass_epa_diff, y = Total_Score_Diff)) +
  geom_point(alpha = 0.5, color = "darkblue") +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Passing EPA/Attempt differential",
       y = "Regular season score differential",
       title = "Relationship between score differential\nand passing efficiency differential") +
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 13))

epa_rush_plot <- ggplot(team_summary_df,
                        aes(x = rush_epa_diff, y = Total_Score_Diff)) +
  geom_point(alpha = 0.5, color = "darkblue") +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Rushing EPA/Attempt differential",
       y = "Regular season score differential",
       title = "Relationship between score differential\n and rushing efficiency differential") +
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 13))

plot_grid(epa_pass_plot, epa_rush_plot, rel_widths = c(1, 1))

# First generate the scattterplot with the points
joint_score_plot <- ggplot(team_summary_df,
                           aes(x = pass_epa_diff, y = rush_epa_diff)) +
  geom_point(aes(color = Total_Score_Diff), alpha = 0.5, size = 4) +
  labs(x = "Passing EPA/Attempt differential",
       y = "Rushing EPA/Attempt differential",
       color = "Regular season\nscore differential",
       title = "Joint distribution of passing and rushing efficiency differentials\ncolored by regular season score differential") +
  # Make the color scale go from darkorange (bad) to midpoint gray (usually 
  # associate white with missing) to darkblue (good):
  scale_color_gradient2(low = "darkorange", mid = "gray", high = "darkblue",
                        midpoint = 0) +
  # Add dashed red lines through the origin:
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  # Add some text for the different quadrants for ease of reading:
  annotate("text", label = "Better pass,\nbetter run", x = .30, y = .28,
           size = 5, color = "darkred") +
  annotate("text", label = "Better pass,\nworse run", x = .30, y = -.15,
           size = 5, color = "darkred") +
  annotate("text", label = "Worse pass,\nbetter run", x = -.30, y = .28,
           size = 5, color = "darkred") +
  annotate("text", label = "Worse pass,\nworse run", x = -.30, y = -.15,
           size = 5, color = "darkred") +
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16))

# Make passing marginal along x axis:
pass_dens <- axis_canvas(joint_score_plot, axis = "x") +
  geom_density(data = team_summary_df, aes(x = pass_epa_diff), 
               fill = "darkblue",
               alpha = 0.5, size = .2)

# Same thing for rushing but along y and with coord_flip = TRUE to make it vertical:
rush_dens <- axis_canvas(joint_score_plot, axis = "y", coord_flip = TRUE) +
  geom_density(data = team_summary_df, aes(x = rush_epa_diff), 
               fill = "darkblue",
               alpha = 0.5, size = .2) +
  coord_flip()

# Now generate by adding these objects to the main plot:
joint_score_plot_1 <- insert_xaxis_grob(joint_score_plot, pass_dens, 
                                        grid::unit(.2, "null"), position = "top")
joint_score_plot_2 <- insert_yaxis_grob(joint_score_plot_1, rush_dens, 
                                        grid::unit(.2, "null"), position = "right")
ggdraw(joint_score_plot_2)

score_diff_reg_eff_model <- function(reg_params, hyperparams, 
                                     score_diff_values, pass_eff_values,
                                     rush_eff_values) {
  # Log likelihood that now uses the function form for mu using the two
  # predictors and their respective values:
  sum(dnorm(score_diff_values, 
            reg_params["alpha"] + 
              reg_params["beta_p"]*pass_eff_values +
              reg_params["beta_r"]*rush_eff_values, 
            reg_params["sigma"], log = TRUE)) +
    # plus the log priors for each parameter results in log posterior:
    dnorm(reg_params["alpha"], 
          hyperparams$alpha_mu, hyperparams$alpha_sd, log = TRUE) + 
    dnorm(reg_params["beta_p"], 
          hyperparams$beta_p_mu, hyperparams$beta_p_sd, log = TRUE) + 
    dnorm(reg_params["beta_r"], 
          hyperparams$beta_r_mu, hyperparams$beta_r_sd, log = TRUE) + 
    dunif(reg_params["sigma"], hyperparams$sigma_a, 
          hyperparams$sigma_b, log = TRUE)
  
}

init_reg_samples <- laplace_approx(log_posterior_fn = score_diff_reg_eff_model,
                                   init_points = c(alpha = 0, beta_p = 0, 
                                                   beta_r = 0, sigma = 10),
                                   n_samples = 10000,
                                   # Vector of hyperparameters:
                                   hyperparams = list(alpha_mu = 0, alpha_sd = 100,
                                                      beta_p_mu = 0, beta_p_sd = 100,
                                                      beta_r_mu = 0, beta_r_sd = 100,
                                                      sigma_a = 0, sigma_b = 200),
                                   # Specify the data for the optimization!
                                   score_diff_values = team_summary_df$Total_Score_Diff,
                                   pass_eff_values = team_summary_df$pass_epa_diff,
                                   rush_eff_values = team_summary_df$rush_epa_diff)

# We'll use the latex2exp package here:
library(latex2exp)
# and ggridges:
library(ggridges)

# First use the gather function to make this simpler to work with:
init_reg_samples %>%
  gather(param, value) %>%
  # only use the beta parameters:
  filter(param %in% c("beta_p", "beta_r")) %>%
  # visualize the distributions for each with density curves:
  ggplot(aes(x = value, y = param)) +
  geom_density_ridges(alpha = 0.7, fill = "darkblue",
                      # add the rugs underneath as well:
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_color = "darkblue",
                      point_alpha = 0.7) +
  # properly label and also change the y axis spacing so the Passing density
  # is along the bottom of the axis
  scale_y_discrete(labels = c("Passing", "Rushing"), expand = c(0.01, 0.01)) +
  # label using the latex symbols:
  labs(x = TeX("$\\beta$ value"),
       y = TeX("$\\beta$ parameter"),
       title = TeX("Sampled posterior distributions for $\\beta_P$ and $\\beta_R$")) +
  # theme settings:
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14))

brian_reg_samples <- laplace_approx(log_posterior_fn = score_diff_reg_eff_model,
                                    init_points = c(alpha = 0, beta_p = 0, 
                                                    beta_r = 0, sigma = 10),
                                    n_samples = 10000,
                                    # Vector of hyperparameters:
                                    hyperparams = list(alpha_mu = 0, alpha_sd = 100,
                                                       beta_p_mu = 336, beta_p_sd = 0.2236,
                                                       beta_r_mu = 492, beta_r_sd = 0.2236,
                                                       sigma_a = 0, sigma_b = 200),
                                    # Specify the data for the optimization!
                                    score_diff_values = team_summary_df$Total_Score_Diff,
                                    pass_eff_values = team_summary_df$pass_epa_diff,
                                    rush_eff_values = team_summary_df$rush_epa_diff)

# First use the gather function to make this simpler to work with:
brian_reg_samples %>%
  gather(param, value) %>%
  # only use the beta parameters:
  filter(param %in% c("beta_p", "beta_r")) %>%
  # visualize the distributions for each with density curves:
  ggplot(aes(x = value, y = param)) +
  geom_density_ridges(alpha = 0.7, fill = "darkblue",
                      # add the rugs underneath as well:
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_color = "darkblue",
                      point_alpha = 0.7) +
  # properly label and also change the y axis spacing so the Passing density
  # is along the bottom of the axis
  scale_y_discrete(labels = c("Passing", "Rushing"), expand = c(0.01, 0.01)) +
  # label using the latex symbols:
  labs(x = TeX("$\\beta$ value"),
       y = TeX("$\\beta$ parameter"),
       title = TeX("Sampled posterior distributions for $\\beta_P$ and $\\beta_R$ accounting for Brian's prior")) +
  # theme settings:
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14))

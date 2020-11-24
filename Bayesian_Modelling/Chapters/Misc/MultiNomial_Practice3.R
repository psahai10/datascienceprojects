# https://github.com/paul-buerkner/brms/issues/338

df <- data.frame(
  x = rep(0:1, each = 50),
  y1 = c(sample(1:10, 50, TRUE), sample(6:15, 50, TRUE)),
  y2 = c(sample(11:20, 50, TRUE), sample(6:15, 50, TRUE)),
  y3 = c(sample(21:30, 50, TRUE), sample(6:15, 50, TRUE)),
  obs = 1:100
)




library(tidyverse)
df_long <- df %>%
  gather("key", "y", y1:y3)

library(nnet)
df$Y <- cbind(df$y1, df$y2, df$y3)
fit1 <- multinom(Y ~ x, df)
summary(fit1)

library(brms)
fit2 <- brm(y ~ key + key:x + (1 | obs), 
            df_long, family = poisson())
summary(fit2)


# https://data.princeton.edu/pop510/pop510slides12.pdf

# Load required libraries
library(rethinking)
library(plotly)
library(ggplot2)
library(tidyr)
library(devtools)
library(plotly)
library(rgl)
library(rstan)


# https://cperretti.com/post/2019/05/01/classifying-fish-species-using-categorical-regression-in-stan/

set.seed(321) # for reproducibility

n_obs <- 100 # 100 observations (fish collected)
n_cat <- 3 # 3 outcome categories (species)
n_cov <- 3 # 3 covariates (fin measurements)
X <- matrix(5 + runif(n_cov*n_obs), n_obs, n_cov) +
  matrix(c(1,0,0), n_obs, n_cov)# fin measurements

B <- matrix(c(4, 0, 0, # effects of fin measurements on score of each species
              0, 4, 0, 
              0, 0, 4), nrow = n_cat, ncol = n_cov, byrow = TRUE)
score <- X %*% B

p <- exp(score)/rowSums(exp(score))


X2plot <- cbind(seq(4,7, l=100), rep(5.5, 100), rep(5.5, 100))
s2plot <- X2plot %*% B
p2plot <- exp(s2plot)/rowSums(exp(s2plot))

df2plot <- 
  data.frame(species1 = p2plot[,1],
             species2 = p2plot[,2],
             species3 = p2plot[,3] + 0.005,
             cov1_val = seq(4,7, l=100)) %>%
  tidyr::gather(species, probability, -cov1_val) 

ggplot(data = df2plot, aes(x = cov1_val, y = probability, color = species)) +
  geom_line() +
  theme_bw() +
  ylab("Probability of species") +
  xlab("Dorsal fin length (cm)") +
  theme(legend.title = element_blank())


outcome <- rep(NA, n_obs)
for (i in 1:n_obs) {
  outcome[i] <- sample(1:n_cat, size = 1, prob = p[i,])  
}

#p <- plot_ly(data = data.frame(cbind(X, paste("species", outcome))), 
#       x = ~X1, y = ~X2, z = ~X3, 
#      color = ~X4) %>%
# add_markers() %>%
# layout(scene = list(xaxis = list(title = 'Dorsal fin length'),
#                      yaxis = list(title = 'Pectoral fin length'),
#                      zaxis = list(title = 'Caudal fin length')))

#htmlwidgets::saveWidget(p, "index.html")


m_fin <- "data {
  int n_obs;
  int n_cat;
  int n_cov;
  int outcome[n_obs];
  matrix[n_obs, n_cov] X;
}

transformed data {
  vector[n_cov] zeros = rep_vector(0, n_cov);
}

parameters {
  matrix[n_cov, n_cat - 1] B_sub;
}

transformed parameters {
  matrix[n_cov, n_cat] B = append_col(zeros, B_sub);
}

model {
  matrix[n_obs, n_cat] score = X * B;
  
  target += normal_lpdf(to_vector(B_sub) | 0, 10); // prior on B_sub
  
  for(i in 1:n_obs) {
    vector[n_cat] p = softmax(to_vector(score[i,]));
    target += categorical_lpmf(outcome[i] | p);
  }
}

generated quantities {
  matrix[n_obs, n_cat] p_est;
  for(i in 1:n_obs) {
    p_est[i,] = to_row_vector(softmax(to_vector(X[i,] * B)));
  }
}"

#dat_list <-list(n_obs = n_obs, 
                #n_cat = n_cat,
                #n_cov = n_cov,
                #outcome = outcome,
                #X = X)

#fit <- stan( model_code=m_fin, data= dat_list, chains = 4, cores = 4, iter=5000)

fit <- stan(model_code = m_fin, 
            data = list(n_obs = n_obs, 
                        n_cat = n_cat,
                        n_cov = n_cov,
                        outcome = outcome,
                        X = X),
            chains = 4, cores = 4, iter=5000)

#precis(fit, 3)

summary(fit, probs = c(.1, .5, .9))$summary[1:6,]

fit_summary1 <- 
  summary(fit, pars = "p_est", probs = c(0.05, 0.95))$summary %>% 
  as.data.frame %>%
  dplyr::mutate(obs_cat = substr(row.names(.), 7, 999)) %>%
  tidyr::separate(obs_cat, c("observation", "category"), ",") %>%
  dplyr::mutate(category = gsub(pattern = "]", replacement = "", x = category),
                observation = as.integer(observation))

fit_summary2 <- 
  summary(fit, pars = "p_est", probs = c(0.05, 0.95))$summary %>% 
  as.data.frame %>%
  dplyr::mutate(obs_cat = substr(row.names(.), 7, 999)) %>%
  tidyr::separate(obs_cat, c("observation", "category"), ",") %>%
  dplyr::mutate(category = gsub(pattern = "]", replacement = "", x = category),
                observation = as.integer(observation)) %>%
  dplyr::left_join({p %>% 
      as.data.frame() %>% 
      dplyr::mutate(observation = 1:n_obs) %>%
      tidyr::gather(category, true_p, -observation) %>%
      dplyr::mutate(category = gsub("V", "", category))})

# data.frame()
df1 <- data.frame(matrix(1:12,3,4),1:3)

# as.data.frame()
df2 <- as.data.frame(matrix(1:12,3,4),1:3)

df1
df2

ggplot(fit_summary, aes(x = true_p, y = mean)) +
  geom_point() +
  geom_abline() +
  facet_wrap(~category) +
  theme_bw() +
  xlab("True probablity") +
  ylab("Estimated probability") 

ggplot(fit_summary, aes(x = true_p, y = mean)) +
  geom_point() +
  geom_abline() +
  facet_wrap(~category) +
  theme_bw() +
  xlab("True probablity") +
  ylab("Estimated probability") +
  geom_errorbar(width=.02, alpha = 0.4, aes(ymin=`5%`, ymax=`95%`))


# https://thinkinator.com/2016/01/12/r-users-will-now-inevitably-become-bayesians/

library(brms)
library(rstan)
library(rstanarm)

rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ()) # Run on multiple cores

set.seed (3875)

ir <- data.frame (scale (iris[, -5]), Species=iris[, 5])

### With improper prior it takes about 12 minutes, with about 40% CPU utilization and fans running,
### so you probably don't want to casually run the next line...

b1 <- brm(Species ~ Petal.Length + Petal.Width + Sepal.Length + Sepal.Width, data=ir,
          family="categorical", prior=c(set_prior ("normal (0, 8)")))


stan_plot(b1$fit, show_density=TRUE) + coord_cartesian (xlim=c(-15, 15))

pp_check(mm, check="scatter", nreps=12) + geom_jitter (width=0.1, height=0.1, color=rgb (0, 0, 0, 0.2))

# file:///C:/Users/psahai/Downloads/ordinal_tutorial.pdf
library(rethinking)

survey_rating <- c("1", "2", "3", "4")


f_list <- list(c(rep(1, 40), rep(2, 54), rep(3, 119), rep(4,55)))
fdf <- data.frame(response=(c(rep(1, 40), rep(2, 54), rep(3, 119), rep(4,55))), view='fundamentalist')

m_list <- list(c(rep(1, 25), rep(2, 41), rep(3, 135), rep(4,71)))
mdf <- data.frame(response=(c(rep(1, 25), rep(2, 41), rep(3, 135), rep(4,71))), view='moderate')

l_list <- list(c(rep(1, 23), rep(2, 31), rep(3, 113), rep(4,122)))
ldf <- data.frame(response=(c(rep(1, 23), rep(2, 31), rep(3, 113), rep(4,122))), view='liberal')

df <- rbind(fdf, mdf)
df <- rbind(df, ldf)


simplehist(df$response, xlim=c(1,4), xlab='response')

pr_k <- table(df$response)/nrow(df)

df$views <- ifelse(df$view=='liberal', 1L, ifelse(df$view=='moderate', 2L, 3L))

delta <- rdirichlet(10, alpha=rep(2,4))

dat <- list(
  R=df$response,
  E=as.integer(df$political_views),
  alpha=rep(2,4)
)

mstem <- ulam(
  alist(
  R ~ ordered_logistic(phi, kappa),
  phi <- bE*sum(delta_j[1:E]),
  kappa ~ normal(0, 1.5),
  bE ~ normal(0,1),
  vector[5]: delta_j <<- append_row(0, delta),
  simplex[4]: delta ~ dirichlet(alpha)
  ), data=dat, chains=4, cores=4
)

precis(mstem, depth=2)

plot(precis(mstem))

delta_labels <- c('liberal', 'moderate', 'conservative')
pairs(mstem, pars='delta', lables=delta_labels)


fit_sc1 <- brm(
  formula = response ~ 1 + views,
  data = df,
  family = cumulative("probit")
)

summary(fit_sc1)

conditional_effects(fit_sc1, "view", categorical = TRUE)

fit_sc2 <- brm(
  formula = response ~ 1 + cs(views),
  data = df,
  family = acat("probit")
)
summary(fit_sc2)
conditional_effects(fit_sc2, "view", categorical = TRUE)

fit_sc3 <- brm(
  formula = response ~ 1 + views,
  data = df,
  family = acat("probit")
)
summary(fit_sc3)
conditional_effects(fit_sc3, "view", categorical = TRUE)

fit_sc4 <- brm(
  formula = bf(response ~ 1 + views) +
    lf(disc ~ 0 + views, cmc = FALSE),
  data = df,
  family = cumulative("probit")
)
summary(fit_sc4)
conditional_effects(fit_sc4, "view", categorical = TRUE)


loo(fit_sc1, fit_sc2, fit_sc3, fit_sc4)


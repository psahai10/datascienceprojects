# https://discourse.mc-stan.org/t/predicted-probabilities-for-multinomial-regression/11504

library(dplyr)
library(brms)

key <- expand.grid(sex = c('Male', 'Female'),
                   age = c('18-34', '35-64', '45-64'))
n <- 1000
set.seed(123)
df <- data.frame(sex = sample(c('Male', 'Female'), n, replace = TRUE),
                 age = sample(c('18-34', '35-64', '45-64'), n, replace = TRUE),
                 stringsAsFactors = FALSE)

age <- model.matrix(~ age, data = df)[, -1]
sex <- model.matrix(~ sex, data = df)[, -1]

mm <- list(age, sex)

#Coefficients
beta_cat1 <- list(age = matrix(c(2, 1)), sex = matrix(-1))
beta_cat2 <- list(age = matrix(c(-3, 2)), sex = matrix(2))

xb1 <- Reduce(rowSums, lapply(seq_along(beta_cat1), function(i) mm[[i]] %*% beta_cat1[[i]]))
xb2 <- Reduce(rowSums, lapply(seq_along(beta_cat2), function(i) mm[[i]] %*% beta_cat2[[i]]))

eta1 <- xb1 + rnorm(length(xb1))
eta2 <- xb2 + rnorm(length(xb2))

sum_exp <- rowSums(cbind(exp(eta1), exp(eta2)))

logit <- function(x){
  exp(x)/(1 + sum_exp)
}

#Class probabilities
p3 <- logit(0)
p2 <- logit(eta2)
p1 <- logit(eta1)

df <- cbind(df, p1 = p1, p2 = p2, p3 = p3)
#Generate outcome variables
ys <- lapply(seq_len(nrow(df)), function(i) t(rmultinom(1, 1, prob = c(p1[i], p2[i], p3[i]))))
ys <- do.call(rbind, ys)
y_names <- paste0('y', 1:3)
colnames(ys) <- y_names

df <- cbind(df, ys)

df_agg <- df %>% group_by(age, sex) %>% summarise_at(.vars = y_names, .funs = sum)
df_agg$N <- rowSums(df_agg[, y_names])
df_agg$response <- with(df_agg, cbind(y1, y2, y3))

brms_fit <- brm(response | trials(N) ~ (1|age) + (1|sex), data = df_agg, family = multinomial())

p_grid <- seq(from=0, to=1, length.out = 1000)
prob_p <- rep(1, 1000)
prob_data <- dbinom(6, 9, p_grid)
posterior <- prob_data*prob_p
posterior <- posterior/sum(posterior)
samples <- sample(p_grid, prob=posterior, size=1e4, replace=TRUE)

df <- data.frame(p_grid= p_grid, p_posterior = posterior)
library(rethinking)
library(ggplot2)
library(dplyr)
library(tidyverse)
ggplot(df, aes(x=number, y=sample)) + geom_point(color='navyblue', size=3, alpha=0.2)

ggplot(df, aes(x = p_grid, y = p_posterior)) +
  geom_point(size = 3, color = "darkblue") + 
  geom_line(color = "darkblue") +
  geom_vline(xintercept = 0.667, color = "darkorange",
             linetype = "dashed", size = 3, alpha = .5) +
  labs(x = "Catch rate", y = "Posterior probability",
       title = "Posterior approximation for\nJuJu's catch rate after one game") +
  theme_bw() + theme(axis.text = element_text(size = 10), 
                     title = element_text(size = 10))

sum(posterior[p_grid>0.5])

sum(samples<0.5)/1e4

quantile(samples, 0.95)

quantile(samples, c(0.025, 0.975))

p_grid <- seq(from=0, to=1, length.out = 1000)
prob_p <- rep(1, 1000)
prob_data <- dbinom(3, 3, p_grid)
posterior <- prob_data*prob_p
posterior <- posterior/sum(posterior)
samples <- sample(p_grid, prob=posterior, size=1e4, replace=TRUE)
df <- data.frame(sample= samples, number=seq(1, 10000))

ggplot(df, aes(x=sample)) + labs(x = "", y = "") + geom_density(alpha=.2, fill="#FF6666")
PI(samples, prob=0.05)
HPDI(samples, prob=0.05)

p_grid[which.max(posterior)]

p_grid[which.min(posterior)]

chainmode(samples, adj=0.01)

chainmode(samples, adj=.0001)

mean(samples)
median(samples)

sum(posterior*abs(0.84-p_grid))

loss <-sapply(p_grid, function(x) sum(posterior*abs(x-p_grid)))

p_grid[which.min(loss)]
p_grid[which.max(loss)]

median(samples)
dbinom(0:2, size=2, prob=0.7)
rbinom(1, size=2, prob=0.7)

dummy_w <- rbinom(1e5, size=9, prob=sampl0.6)

table(dummy_w)/1e5

simplehist(dummy_w, xlab='dummy water count')

######################### PRACTICE ######################### PRACTICE ############################

p_grid <- seq(from=0, to=1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(6, 9, p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)
set.seed(100)
simulation <- rbinom(1e5, size=9, prob=0.7)
norm_sim <- sapply(simulation, function(x) x/9)
samples <- sample(p_grid, prob=posterior, size=1e4, replace=TRUE)

## 3E1
sum(posterior[p_grid<0.2])
sum(samples<0.2)/1e4

## 3E2
sum(posterior[p_grid>0.8])
sum(samples>0.8)/1e4

## 3E3
sum(posterior[p_grid>0.8])
sum(samples>0.8)/1e4

# 3E4 & 3E5
HPDI(samples, prob=0.6)

#3E6 
HPDI(samples, prob=0.66) 

# 3E7
PI(samples, prob=0.66)

#3M1
p_grid <- seq(0, 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(8, 15, prob=p_grid)
posterior <- prior*likelihood
posterior <- posterior/sum(posterior)

#3M2
samples <- sample(p_grid, prob=posterior, size=1e4, replace=TRUE)
plot(samples)
dens(samples)
HPDI(samples, prob=0.9)

#3M3
w<-rbinom(1e4, size=15, prob=samples)
n <- length(w)
sum(w==8)/n

#3M4
dummy_w <- dbinom(6, 9, prob=posterior)
sum(dummy_w)

#3M5
p_grid <- seq(0, 1, length.out = 1000)
prior <- rep(0.5, 1000)
likelihood <- dbinom(8, 15, prob=p_grid)
posterior <- prior*likelihood
posterior <- posterior/sum(posterior)

# HARD
data(homeworkch3)
sum(birth1) + sum(birth2)
length(birth1) + length(birth2)
total_birth = table(birth1) + table(birth2)

# 3H1
p_grid <- seq(0, 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(111, 200, prob=p_grid)
posterior <- prior*likelihood
posterior <- posterior/sum(posterior)
samples <- sample(p_grid, prob=posterior, size=1e4, replace=TRUE)
plot(samples)
dens(samples)
max(samples)
HPDI(samples, prob=0.02)

# 3H2
HPDI(samples, prob=0.5)
HPDI(samples, prob=0.89)
HPDI(samples, prob=0.97)

# 3H3
dummy_w <-rbinom(1e5, 200,  prob=0.485)
table(dummy_w)/1e5
simplehist(dummy_w)
dens(dummy_w)

#3H4
table(birth1)
p_grid <- seq(0, 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(51, 100, prob=p_grid)
posterior <- prior*likelihood
posterior <- posterior/sum(posterior)
samples <- sample(p_grid, prob=posterior, size=1e4, replace=TRUE)
plot(samples)
dens(samples)
max(samples)
HPDI(samples, prob=0.02)

# 3H5
table(birth1)
table(birth2)
df <- data.frame(first_birth=birth1, second_birth=birth2)
library(dplyr)
df_first_birth_female <- df %>% filter_at(vars(first_birth), any_vars(. %in% c('0')))   
df_first_birth_male <- df %>% filter_at(vars(first_birth), any_vars(. %in% c('1')))   
table(df_first_birth_female)
table(df_first_birth_male)

# 3H5 Part 2
dummy_g <-rbinom(1e5, 100,  prob=0.49)
table(dummy_w)/1e5
simplehist(dummy_g)
dens(dummy_g)

# 3H5 Part 3
dummy_g <-rbinom(1e5, 100,  prob=0.49)
table(dummy_w)/1e5
simplehist(dummy_g)
dens(dummy_g)

###########
x<-rbinom(1e5, size=49, prob=.49)
sum(x>=39)/1e5

###########
p_grid <- seq(0, 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(39, 49, prob=p_grid)
posterior <- prior*likelihood
posterior <- posterior/sum(posterior)
samples <- sample(p_grid, prob=posterior, size=1e4, replace=TRUE)
plot(samples)
dens(samples)

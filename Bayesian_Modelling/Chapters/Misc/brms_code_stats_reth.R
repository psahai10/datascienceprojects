# https://solomonkurz.netlify.app/post/

# https://bookdown.org/connect/#/apps/1850/access

library(rethinking)
library(brms)
library(loo)
library(tidybayes)
library("wesanderson")

# simulate career choices among 500 individuals
n      <- 500           # number of individuals
income <- c(1,2,5)           # expected income of each career
score  <- 0.5 * income  # scores for each career, based on income

# next line converts scores to probabilities
p <- softmax(score[1], score[2], score[3])

# now simulate choice
# outcome career holds event type values, not counts
career <- rep(NA, n)  # empty vector of choices for each individual

set.seed(34302)
# sample chosen career for each individual
for(i in 1:n) career[i] <- sample(1:3, size = 1, prob = p)

career %>%
  as_tibble() %>%
  ggplot(aes(x = value %>% as.factor())) +
  geom_bar(size = 0, fill = wes_palette("Moonrise2")[2])

dat_list <- list( N=N , K=3 , career=career , 
                  career_income= ifelse(career==3, 5, ifelse(career==2, 2, 1)))

b11.13 <-
  brm(data = list(career = dat_list$career,
                  career_income= dat_list$career_income), 
      family = categorical(link = logit),
      career ~ career_income + (1|p|career_income),
      prior = c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 0.5), class = b)),
      iter = 2500, warmup = 500, cores = 4, chains = 4,
      seed = 34302)

post <- posterior_samples(b11.13)
s1 <- with(post, b_mu2_Intercept + b_mu2_Intercept*2)
s2_orig <- with(post, b_mu3_Intercept + b_mu3_Intercept*5)
s2_new <- with(post, b_mu3_Intercept + b_mu3_Intercept*5*2)

p_orig <- sapply( 1:length(post$b_mu2_Intercept), function(i)
  softmax( c(s1[i], s2_orig[i], 0) ) )

p_new <- sapply( 1:length(post$b_mu2_Intercept), function(i)
  softmax( c(s1[i], s2_new[i], 0) ) )

p_diff <- p_new[2,] - p_orig[2,]
precis(p_diff, hist=FALSE)


# https://cran.r-project.org/web/packages/brms/vignettes/brms_multivariate.html

print(b11.13)

n <- 500

set.seed(34302)
# simulate family incomes for each individual
family_income <- runif(n)

# assign a unique coefficient for each type of event
b      <- c(-2, 0, 2)
career <- rep(NA, n)  # empty vector of choices for each individual

for (i in 1:n) {
  score     <- 0.5 * (1:3) + b * family_income[i]
  p         <- softmax(score[1], score[2], score[3])
  career[i] <- sample(1:3, size = 1, prob = p)
}

fit2 <- brm(career ~ family_income,
  data = list(career = career,  # note how we used a list instead of a tibble
              family_income = family_income), 
  chains = 2, cores = 2, 
)

m11.14 <- quap(
  alist(
    career ~ dcategorical(softmax(s1,s2,s3)),
    s1 <- a1 + b1*family_income,
    s2 <- a2 + b2*family_income,
    s3 <- 0,
    c(a1,a2) ~ dnorm(0,1.5),
    c(b1,b2) ~ dnorm(0,1)),
  data=list(career=career,family_income=family_income))

precis(m11.14)

b11.14 <-
  brm(data = list(career= career,  # note how we used a list instead of a tibble
                  family_income = family_income), 
      family = categorical(link = logit),
      career ~ 1 + family_income,
      prior = c(prior(normal(0, 1.5), class = Intercept),
                prior(normal(0, 1), class = b)),
      iter = 2500, warmup = 500, cores = 3, chains = 3,
      seed = 34302)

print(b11.14)

# binomial model of overall admission probability
b_binom <-
  brm(data = d, family = binomial,
      admit | trials(applications) ~ 1,
      prior(normal(0, 100), class = Intercept),
      iter = 2000, warmup = 1000, cores = 3, chains = 3,
      seed = 10)

# Poisson model of overall admission rate and rejection rate
b_pois <-
  brm(data = d %>%
        mutate(rej = reject),  # 'reject' is a reserved word
      family = poisson,
      mvbind(admit, rej) ~ 1,
      prior(normal(0, 100), class = Intercept),
      iter = 2000, warmup = 1000, cores = 3, chains = 3,
      seed = 10)
# extract the samples
post <- posterior_samples(b_pois)

# wrangle
post %>%
  transmute(admit  = exp(b_admit_Intercept), 
            reject = exp(b_rej_Intercept)) %>% 
  gather() %>% 
  
  # plot
  ggplot(aes(x = value, y = key, fill = key)) +
  geom_halfeyeh(point_interval = median_qi, .width = .95,
                color = wes_palette("Moonrise2")[4]) +
  scale_fill_manual(values = c(wes_palette("Moonrise2")[1],
                               wes_palette("Moonrise2")[2])) +
  labs(title = " Mean admit/reject rates across departments",
       x     = "# applications",
       y     = NULL) +
  theme(legend.position = "none",
        axis.ticks.y    = element_blank())


# simulate
n <- 100
set.seed(10)
x <- runif(n)

set.seed(10)
y <- rgeom(n, prob = inv_logit_scaled(-1 + 2 * x))

list(y = y, x = x) %>%
  as_tibble() %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(size = 3, alpha = 0.5)

b10.18 <-
  brm(data = list(y = y, x = x), 
      family = geometric(link = log),
      y ~ 0 + intercept + x,
      prior = c(prior(exponential(1), class = b, coef = intercept),
                prior(normal(0, 1), class = b)),
      iter = 2500, warmup = 500, chains = 2, cores = 2,
      seed = 10)

plot(marginal_effects(b10.18),
     points = T,
     point_args = c(size = 2, alpha = 2/3),
     line_args  = c(color = wes_palette("Moonrise2")[1],
                    fill  = wes_palette("Moonrise2")[1]))


# estimate
m10.18 <- quap(
  alist(
    y ~ dgeom( p ),
    logit(p) <- a + b*x,
    a ~ dnorm(0,10),
    b ~ dnorm(0,1)
  ),
  data=list(y=y,x=x) )
precis(m10.18)

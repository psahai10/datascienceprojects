# 4.4.3.5 Prediction intervals.
# Much as brms::fitted() was our analogue to rethinking::link(), brms::predict() is our 
# analogue to rethinking::sim()

library(rethinking)
library(brms)
library(dplyr)

data("Howell1"); d <- Howell1; d2 <- d[ d$age >= 18 , ]
xbar <- mean( d2$height )

m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b * (weight - xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ), data=d2
)

precis(m4.3)


b4.3 <- brm(data = d2, family = gaussian,
        height ~ 1 + weight,
        prior = c(prior(normal(156, 100), class = Intercept),
                  prior(normal(0, 10), class = b),
                  prior(uniform(0, 50), class = sigma)),
        iter = 41000, warmup = 40000, chains = 4, cores = 4,
         seed = 4)

plot(b4.3)

posterior_summary(b4.3)[1:3, ]

posterior_samples(b4.3) %>%
  select(-lp__) %>%
  cor() %>%
  round(digits = 2)

d2 <- d2 %>%
  mutate(weight_c = weight - mean(weight))

b4.4 <- 
  brm(data = d2, family = gaussian,
      height ~ 1 + weight_c,
      prior = c(prior(normal(178, 100), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 50), class = sigma)),
      iter = 46000, warmup = 45000, chains = 4, cores = 4,
      seed = 4)

plot(b4.4)

posterior_summary(b4.4)[1:3, ]

posterior_samples(b4.4) %>%
  select(-lp__) %>%
  cor() %>%
  round(digits = 2)

d2 %>%
  ggplot(aes(x = weight, y = height)) +
  geom_abline(intercept = fixef(b4.3)[1], 
              slope     = fixef(b4.3)[2]) +
  geom_point(shape = 1, size = 2, color = "royalblue") +
  theme_bw() +
  theme(panel.grid = element_blank())

post <- posterior_samples(b4.3)

post %>%
  slice(1:5)  # this serves a similar function as `head()`

priors = c(prior(normal(178, 100), class = Intercept),
          prior(normal(0, 10), class = b),
          prior(cauchy(0, 1), class = sigma))

n <- 50

b.10 <- brm(data = d2 %>%
        slice(1:n),  # note our tricky use of `n` and `slice()`
        family = gaussian,
        height ~ 1 + weight,
        prior = priors,
        iter = 2000, warmup = 1000, chains = 4, cores = 4,
        seed = 4)

post10  <- posterior_samples(b.10)


  ggplot(data =  d2[1:50 , ], 
         aes(x = weight, y = height)) +
  geom_abline(intercept = post10[1:20, 1], 
              slope     = post10[1:20, 2],
              size = 1/3, alpha = .3) +
  geom_point(shape = 1, size = 2, color = "royalblue") +
  coord_cartesian(xlim = range(d2$weight),
                  ylim = range(d2$height)) +
  labs(subtitle = "N = 10") +
  theme_bw() +
  theme(panel.grid = element_blank())

  
data("WaffleDivorce")
d <- list()
d$A <- standardize( WaffleDivorce$MedianAgeMarriage)
d$D <- standardize( WaffleDivorce$Divorce)
d$M <- standardize( WaffleDivorce$Marriage)

m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 ),
    ## A -> M
    M ~ dnorm( mu_M , sigma_M ),
    mu_M <- aM + bAM*A,
    aM ~ dnorm( 0 , 0.2 ),
    bAM ~ dnorm( 0 , 0.5 ),
    sigma_M ~ dexp( 1 )
  ) , data = d )

## R code 5.20
A_seq <- seq( from=-2 , to=2 , length.out=30 )

## R code 5.21
# prep data
sim_dat <- data.frame( A=A_seq )

# simulate M and then D, using A_seq
s <- sim( m5.3_A , data=sim_dat , vars=c("M","D") )

## R code 5.22
plot( sim_dat$A , colMeans(s$D) , ylim=c(-2,2) , type="l" ,
      xlab="manipulated A" , ylab="counterfactual D"  )
shade( apply(s$D,2,PI) , sim_dat$A )
mtext( "Total counterfactual effect of A on D" )

## R code 5.23
# new data frame, standardized to mean 26.1 and std dev 1.24
sim2_dat <- data.frame( A=(c(20,30)-26.1)/1.24 )
s2 <- sim( m5.3_A , data=sim2_dat , vars=c("M","D") )
mean( s2$D[,2] - s2$D[,1] )

## R code 5.24
sim_dat <- data.frame( M=seq(from=-2,to=2,length.out=30) , A=0 )
s <- sim( m5.3_A , data=sim_dat , vars="D" )

plot( sim_dat$M , colMeans(s) , ylim=c(-2,2) , type="l" ,
      xlab="manipulated M" , ylab="counterfactual D"  )
shade( apply(s,2,PI) , sim_dat$M )
mtext( "Total counterfactual effect of M on D" )


b5.3 <- 
  brm(data = d, family = gaussian,
      D ~ 1 + M + A,
      prior = c(prior(normal(10, 10), class = Intercept),
                prior(normal(0, 1), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 500, chains = 4, cores = 4,
      seed = 5)

nd <- tibble(M = seq(from=-3, to=3, length.out = 30), A = mean(d$A))


fitted(b5.3, newdata = nd) %>% 
  as_tibble() %>% 
  rename(f_ll = Q2.5,
         f_ul = Q97.5) %>% 
  bind_cols(
    predict(b5.3, newdata = nd) %>% 
      as_tibble() %>% 
      transmute(p_ll = Q2.5,
                p_ul = Q97.5),
    nd) %>% 
  ggplot(aes(x = M, y = Estimate)) +
  geom_ribbon(aes(ymin = p_ll, ymax = p_ul),
              fill = "firebrick", alpha = 1/5) +
  geom_smooth(aes(ymin = f_ll, ymax = f_ul),
              stat = "identity",
              fill = "firebrick", color = "firebrick4", alpha = 1/5, size = 1/4) +
  coord_cartesian(xlim = range(d$M),
                  ylim = c(6, 14)) +
  labs(subtitle = "Counterfactual plot for which\nMedianAgeMarriage_s = 0",
       y = "Divorce") +
  theme_bw() +
  theme(panel.grid = element_blank())    


# new data

nd <- tibble(A = seq(from =-3, to =3.5, length.out = 30), M= mean(d$M))

# `fitted()` + `predict()`
fitted(b5.3, newdata = nd) %>% 
  as_tibble() %>% 
  rename(f_ll = Q2.5,
         f_ul = Q97.5) %>% 
  bind_cols(
    predict(b5.3, newdata = nd) %>% 
      as_tibble() %>% 
      transmute(p_ll = Q2.5,
                p_ul = Q97.5),
    nd
  ) %>% 
  # plot
  ggplot(aes(x = A, y = Estimate)) +
  geom_ribbon(aes(ymin = p_ll, ymax = p_ul),
              fill = "firebrick", alpha = 1/5) +
  geom_smooth(aes(ymin = f_ll, ymax = f_ul),
              stat = "identity",
              fill = "firebrick", color = "firebrick4", alpha = 1/5, size = 1/4) +
  coord_cartesian(xlim = range(d$A),
                  ylim = c(6, 14)) +
  labs(subtitle = "Counterfactual plot for which\nMarriage_s = 0",
       y = "Divorce") +
  theme_bw() +
  theme(panel.grid = element_blank())  

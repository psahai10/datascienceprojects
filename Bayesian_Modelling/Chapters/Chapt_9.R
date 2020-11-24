library(rstan)
library(rethinking)
data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[ complete.cases(d$rgdppc_2000) , ]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
dd$cid <- ifelse( dd$cont_africa==1 , 1 , 2 )

## R code 9.12
m8.3 <- quap(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp( 1 )
  ) , data=dd )
precis( m8.3 , depth=2 )

## R code 9.13
dat_slim <- list(
  log_gdp_std = dd$log_gdp_std,
  rugged_std = dd$rugged_std,
  cid = as.integer( dd$cid )
)
str(dat_slim)

## R code 9.14
m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp( 1 )
  ) , data=dat_slim , chains=4, cores=4 )


show(m9.1)

precis( m9.1 , depth=2 )

pairs(m9.1)

traceplot(m9.1, chains=1)

trankplot(m9.1)

stancode(m9.1)


y <- c(-1,1)

set.seed(11)

m9.2 <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- alpha,
    alpha ~ dnorm(0, 1000),
    sigma ~ dexp(0.0001)
  ), data=list(y=y), chains=3
)

precis(m9.2)
traceplot(m9.2, chains=1)
trankplot(m9.2)
pairs(m9.1@stanfit)
pairs(m9.2@stanfit)


m9.2a <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- alpha,
    alpha ~ dnorm(0, 0.5),
    sigma ~ dexp(0.5)
  ), data=list(y=y), chains=3
)
precis(m9.2a)
traceplot(m9.2a, chains=1)
trankplot(m9.2a)
pairs(m9.2a@stanfit)

set.seed(41)
y <- rnorm(100, mean=0, sd=1)

set.seed(384)

m9.4 <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- alpha1 + alpha2,
    alpha1 ~ dnorm(0, 1000),
    alpha2 ~ dnorm(0, 1000),
    sigma ~ dexp(1)
  ), data=list(y=y), chains=3
)
precis(m9.4)

m9.5 <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- alpha1 + alpha2,
    alpha1 ~ dnorm(0, 1),
    alpha2 ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=list(y=y), chains=3
)
precis(m9.5)

m.9m1 <- ulam(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dunif(0, 1)
  ) , data=dat_slim , chains=4, cores=4 )
precis(m9.1)
precis(m.9m1)


mp <- ulam(
  alist(
    a ~ dnorm(0,1),
    b ~ dcauchy(0,1)
  ), data=list(y=1), chains=3
)
precis(mp)
traceplot(mp, chains=1)
trankplot(mp)

data("WaffleDivorce")
d <- WaffleDivorce

d$D <-standardize(d$Divorce)
d$M <-standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

sd(d$MedianAgeMarriage)


d_slim <- list(
  D = d$D ,
  M = d$M ,
  A = d$A
)

str(d_slim)

m5.1 <- ulam(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data= d_slim, chains=4, log_lik=TRUE
)

m5.2 <- ulam(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data= d_slim, chains=4, log_lik=TRUE
)

m5.3 <- ulam(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M + bA * A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data= d_slim, chains=4, log_lik=TRUE
)

precis(m5.1)
precis(m5.2)
precis(m5.3)


compare(m5.1, m5.2, m5.3, func=WAIC)

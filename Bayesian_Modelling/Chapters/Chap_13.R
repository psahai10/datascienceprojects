library(rethinking)
data("reedfrogs")
d <- reedfrogs
str(d)
head(d, 15)

d$tank <- 1:nrow(d)

dat <- list(
  S = d$surv,
  N = d$density,
  tank = d$tank
)

m13.1 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a[tank],
    a[tank] ~ dnorm(0, 1.5)
  ), data=dat, chains=4, cores=4, log_lik = TRUE
)

m13.1b <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a[tank],
    a[tank] ~ dnorm(1.35, 1.62)
  ), data=dat, chains=4, cores=4, log_lik = TRUE
)

precis(m13.1, 2)
length(precis(m13.1, 2)$mean)

m13.2 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a[tank],
    a[tank] ~ dnorm(a_bar, sigma),
    a_bar ~ dnorm(0, 1.5),
    sigma ~ dexp(1)
  ), data=dat, chains=4, cores=4, log_lik = TRUE
)

precis(m13.2, 2)
length(precis(m13.2, 2)$mean)

compare(m13.1, m13.2)

post <- extract.samples(m13.2)

d$propsurv.est <- logistic(apply(post$a, 2, mean))

plot(d$propsurv, ylim=c(0,1), pch=16, xant='n',
     xlab='tank', yalb='prop survival', col=rangi2)
axis(1,at=c(1,16,32,48), labels=c(1,16,32,48))

points(d$propsurv.est)

abline(h=mean(inv_logit(post$a_bar)), lty=2)

abline( v=16.5, lwd=0.5)
abline( v=32.5, lwd=0.5)

text(8,0,'small tanks')
text(16+8, 0, 'medium tanks')
text(32+8, 0, 'large tanks')

plot(NULL, xlim=c(-3,4), ylim=c(0, 0.35),
     xlab='log-odds survive', ylab='Density')
for (i in 1:200)
  curve( dnorm(x, post$a_bar[i],post$sigma[i]), add=TRUE,
         col=col.alpha('black', 0.2))
sim_tanks <- rnorm(8000, post$a_bar, post$sigma)

dens(inv_logit(sim_tanks), lwd=2, adj=0.1)

library(rethinking)
a_bar <- 1.5
sigma <- 1.5
n_ponds <- 60

Ni <- as.integer(rep(c(5,10,25,35), each=15))

set.seed(5005)
a_pond <- rnorm(n_ponds, mean=a_bar, sd=sigma)
#a_pond <- rnorm(15, mean=1.15, sd=1)
#b_pond <- rnorm(15, mean=1.25, sd=1.4)
#c_pond <- rnorm(25, mean=1.65, sd=1)
#d_pond <- rnorm(25, mean=1.50, sd=1.5)
#e_pond <- rnorm(25, mean=1.15, sd=1.75)
#a_pond <-  do.call(c, list(d_pond, c_pond, e_pond ))
#pond <- ifelse(pond<0, 0, pond)

dsim <- data.frame(pond=1:n_ponds, Ni=Ni, true_a=a_pond)

inv_logit(dsim$true_a)
logistic((dsim$true_a))

dsim$Si <- rbinom(n_ponds, prob=logistic(dsim$true_a), size=dsim$Ni)
dsim$p_nopool <-dsim$Si/dsim$Ni

dat <- list(Si=dsim$Si, Ni=dsim$Ni, pond=dsim$pond)

m13.3 <- ulam(
  alist(
    Si ~ dbinom(Ni, p),
    logit(p) <- a_pond[pond],
    a_pond[pond] ~ dnorm(a_bar, sigma),
    a_bar ~ dnorm(0, 1.5),
    sigma ~ dexp(1)
  ), data=dat, chains=4, cores=4, log_lik = TRUE
)
precis(m13.3, 2)

post <- extract.samples(m13.3)
dsim$p_partpool <- apply(inv_logit(post$a_pond), 2, mean)

dsim$p_true <- inv_logit( dsim$true_a )

nopool_error <- abs(dsim$p_nopool - dsim$p_true)
partpool_error <- abs(dsim$p_partpool - dsim$p_true)

plot(1:60, nopool_error, xlab='pond', ylab='abs error', col=rangi2, pch=16)
points(1:60, partpool_error)


library(rethinking)
library(rstan)

data(chimpanzees)
d <- chimpanzees
d$treatment <- d$prosoc_left+2*d$condition

dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  block_id = d$block,
  treatment = as.integer(d$treatment) )

set.seed(13)

m13.4 <- ulam(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a[actor] + g[block_id] + b[treatment],
    b[treatment] ~ dnorm(0, 0.5),
    ## Adaptive Priors
    a[actor] ~ dnorm(a_bar, sigma_a),
    g[block_id] ~ dnorm(0, sigma_b),
    # Hyper Priors
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1),
    sigma_b ~ dexp(1)
  ), data= dat_list, chains=4, cores=4, log_lik = TRUE
)
precis(m13.4)

dat_list2 <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = as.integer(d$treatment))

m13.5 <- ulam(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a[actor] + b[treatment],
    b[treatment] ~ dnorm(0, 0.5),
    ## Adaptive Priors
    a[actor] ~ dnorm(a_bar, sigma_a),
    # Hyper Priors
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1)
  ), data= dat_list2, chains=4, cores=4, log_lik = TRUE
)

m13.6 <- ulam(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a[actor] + g[block_id] + b[treatment],
    b[treatment] ~ dnorm(0, simga_c),
    ## Adaptive Priors
    a[actor] ~ dnorm(a_bar, sigma_a),
    g[block_id] ~ dnorm(0, sigma_b),
    # Hyper Priors
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1),
    sigma_b ~ dexp(1),
    simga_c ~ dexp(1)
  ), data= dat_list, chains=4, cores=4, log_lik = TRUE
)

m13.4nc <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a_bar + z[actor]*sigma_a + # actor intercepts
      x[block_id]*sigma_g +      # block intercepts
      b[treatment] ,
    b[treatment] ~ dnorm( 0 , 0.5 ),
    z[actor] ~ dnorm( 0 , 1 ),
    x[block_id] ~ dnorm( 0 , 1 ),
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1),
    gq> vector[actor]:a <<- a_bar + z*sigma_a,
    gq> vector[block_id]:g <<- x*sigma_g
  ) , data=dat_list , chains=4 , cores=4 )

m13.6nc <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a_bar + z[actor]*sigma_a + # actor intercepts
      x[block_id]*sigma_g +      # block intercepts
      b[treatment]*sigma_k ,
    b[treatment] ~ dnorm( 0 , 1 ),
    z[actor] ~ dnorm( 0 , 1 ),
    x[block_id] ~ dnorm( 0 , 1 ),
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1),
    sigma_k, dexp(1),
    gq> vector[actor]:a <<- a_bar + z*sigma_a,
    gq> vector[block_id]:g <<- x*sigma_g,
    gq> vector[treatment]:k <<- b*sigma_k
  ) , data=dat_list , chains=4 , cores=4 )

## R code 13.21
library(rethinking)
data(chimpanzees)
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition

dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  block_id = d$block,
  treatment = as.integer(d$treatment) )

set.seed(13)
m13.4 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + g[block_id] + b[treatment] ,
    b[treatment] ~ dnorm( 0 , 0.5 ),
    ## adaptive priors
    a[actor] ~ dnorm( a_bar , sigma_a ),
    g[block_id] ~ dnorm( 0 , sigma_g ),
    ## hyper-priors
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1)
  ) , data=dat_list , chains=4 , cores=4 , log_lik=TRUE )


chimp <- 2
d_pred <- list(
  actor = rep(chimp, 4),
  treatment = 1:4,
  block_id = rep(1,4)
)

p <- link( m13.4, data=d_pred)
p_mu <- apply(p,2,mean)
p_ci <- apply(p, 2, PI)


post <- extract.samples(m13.4)
str(post)
head(post)

dens(post$a[,4])

p_link <- function(treatment, actor=1, block_id=1) {
  logodds <- with( post ,
                   a[,actor] + g[,block_id] + b[,treatment])
  return( inv_logit(logodds) )
}

p_raw <- sapply(1:4, function(i) p_link(i, actor=4, block_id = 1) )
p_mu_raw <- apply( p_raw, 2, mean )
p_ci_raw <- apply( p_raw, 2, PI )


## R code 13.36
p_link_abar <- function( treatment ) {
  logodds <- with( post , a_bar + b[,treatment] )
  return( inv_logit(logodds) )
}

## R code 13.37
post <- extract.samples(m13.4)
p_raw <- sapply( 1:4 , function(i) p_link_abar( i ) )
p_mu <- apply( p_raw , 2 , mean )
p_ci <- apply( p_raw , 2 , PI )

plot( NULL , xlab="treatment" , ylab="proportion pulled left" ,
      ylim=c(0,1) , xaxt="n" , xlim=c(1,4) )
axis( 1 , at=1:4 , labels=c("R/N","L/N","R/P","L/P") )
lines( 1:4 , p_mu )
shade( p_ci , 1:4 )

## R code 13.38
a_sim <- with( post , rnorm( length(post$a_bar) , a_bar , sigma_a ) )
p_link_asim <- function( treatment ) {
  logodds <- with( post , a_sim + b[,treatment] )
  return( inv_logit(logodds) )
}
p_raw_asim <- sapply( 1:4 , function(i) p_link_asim( i ) )

## R code 13.39
plot( NULL , xlab="treatment" , ylab="proportion pulled left" ,
      ylim=c(0,1) , xaxt="n" , xlim=c(1,4) )
axis( 1 , at=1:4 , labels=c("R/N","L/N","R/P","L/P") )
for ( i in 1:200 ) lines( 1:4 , p_raw_asim[i,] , col=grau(0.25) , lwd=2 )

## R code 13.40


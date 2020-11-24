library(rethinking)
library(rstan)
prob_drink <- 0.2
rate_work <- 1

N <- 365

set.seed(365)

drink <- rbinom( N, 1, prob_drink )

y <- (1-drink)*rpois(N, rate_work)

simplehist(y, xlab='manuscripts completed', lwd=4)
zeros_drink <- sum(drink)
zeros_work <-sum(y==0 & drink==0)
zeros_total <- sum(y==0)
lines(c(0,0), c(zeros_work, zeros_total), lwd=4, col=rangi2)

dat <- list(y=y)

m12.3 <- ulam(
  alist(
    y ~ dzipois( p, lambda ),
    logit( p ) <- ap,
    log( lambda ) <- al,
    ap ~ dnorm(0,1),
    al ~ dnorm(1, 0.5)
  ), data=dat, chains=4, cores=4)

precis(m12.3)

post <- extract.samples(m12.3)
mean(inv_logit(post$ap))
mean(exp(post$al))
library(rethinking)
data("Trolley")
d <- Trolley
head(Trolley)

simplehist(d$response, xlim=c(1,7), xlab='response')

pr_k <- table(d$response)/nrow(d)

cum_pr_k <- cumsum(pr_k)

plot(1:7, cum_pr_k, type='b', xlab='response')

odds <- function(x) x/(1-x)
cum_odds <- round(lodds<-odds(cum_pr_k),3)

logit <- function(x) log(x/(1-x))
log_cum_pr_k <- round(lco<-logit(cum_pr_k), 3 )

plot(1:7, log_cum_pr_k, type='b', xlab='response')


m12.4 <- ulam(
  alist(
    R ~ dordlogit(0, cutpoints),
    cutpoints ~ dnorm(0, 1.5)
  ), data=list(R=d$response), chains=4, cores=4
)

m12.4q <- quap(
  alist(
    response ~ dordlogit( 0 , c( a1 , a2 , a3 , a4 , a5 , a6 ) ),
    c( a1 , a2 , a3 , a4 , a5 , a6 ) ~ dnorm( 0 , 1.5 ),
  ), data = d , start = list( a1 = -2 , a2 = -1 , a3 = 0, a4 = 1 , a5 = 2 , a6 = 3 ))

precis(m12.4, 2)
precis(m12.4q, 2)

round(inv_logit(coef(m12.4)), 3)
round(inv_logit(coef(m12.4q)), 3)

round(pk <- dordlogit(1:7, 0, coef(m12.4)), 2 )
round(pk <- dordlogit(1:7, 0, coef(m12.4q)), 2 )

dat <- list(
  R = d$response,
  A = d$action,
  I = d$intention,
  C = d$contact
)
m12.5 <- ulam(
  alist(
    R ~ dordlogit( phi, cutpoints ),
    phi <- bA*A + bC*C + BI*I ,
    BI <- bI + bIA * A + bIC*C,
    c(bA,bI,bC,bIA,bIC) ~ dnorm(0, 0.5),
    cutpoints ~ dnorm(0, 1.5)
  ), data=dat, chains=4, cores=4
)
precis(m12.5)

plot(precis(m12.5), xlim=c(-1.4,0))

plot(NULL, type='n', xlab='intention', ylab='probability',
     xlim=c(0,1), ylim=c(0,1), xaxp=c(0,1,1), yaxp=c(0,1,2))
kA <-0
kC <- 0
kI <- 0:1
pdat <- data.frame(A=kA, C=kC, I=kI)
phi <- link(m12.5, data=pdat)$phi

post <- extract.samples(m12.5)
for (s in 1:50) {
  pk <- pordlogit(1:6, phi[s,], post$cutpoints[s,])
  for (i in 1:6 ) lines(kI, pk[,i], col=grau(0.1))
}


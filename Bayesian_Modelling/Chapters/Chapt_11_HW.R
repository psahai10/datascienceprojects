library(rethinking)
data(Kline)
d <- Kline

d$P <- scale(log(d$population))
d$contact_id <- ifelse(d$contact=='high',2,1)

curve(dlnorm(x,0,10), from=0,to=100,n=200)

a <- rnorm(1e4, 0, 10)
lambda <- exp(a)
mean(lambda)


curve(dlnorm(x,3,0.5), from=0,to=100,n=200)

N <- 100
a <- rnorm(N, 3, 0.5)
b <- rnorm(N, 0, 0.1)
plot(NULL, xlim=c(-2,2), ylim=c(0,100))
for (i in 1:N) curve( exp(a[i] + b[i]*x), add=TRUE, col=grau())


x_seq <- seq(from=log(100), to=log(2e5), length.out=100)
lambda <- sapply(x_seq, function(k) exp(a+b*k))
plot(NULL, xlim=range(x_seq), ylim=c(0,500), xlab='log population', ylab='total tools')
for (i in 1:N) lines(x_seq, lambda[i,], col=grau(), lwd=1.25)

plot(NULL, xlim=range(exp(x_seq)), ylim=c(0,500), xlab='population', ylab='total tools')
for (i in 1:N) lines(exp(x_seq), lambda[i,], col=grau(), lwd=1.25)


dat <- list(
  T=d$total_tools,
  log_P=d$P,
  cid=d$contact_id,
  P=d$population
  
)

m11.9 <- ulam(
  alist(
    T ~ dpois( lambda ),
    log(lambda) <- a,
    a ~ dnorm(3, 0.5)
  ), data=dat, chains=4, log_lik = TRUE
)

m11.10 <- ulam(
  alist(
    T ~ dpois( lambda ),
    #log(lambda) <- a[cid] + b[cid]*P,
    lambda <- exp(a[cid])*P^b[cid]/g,
    a[cid] ~ dnorm(1, 1),
    b[cid] ~ dexp(1),
    g ~ dexp(1)
  ), data=dat, chains=4, cores=4, log_lik = TRUE
)

m11.10b <- ulam(
  alist(
    T ~ dpois( lambda ),
    #log(lambda) <- a[cid] + b[cid]*P,
    lambda <- exp(a[cid])*P^b[cid]/(g-exp(a[cid])),
    a[cid] ~ dnorm(1, 1),
    b[cid] ~ dexp(1),
    g ~ dexp(1)
  ), data=dat, chains=4, cores=4, log_lik = TRUE
)

m11.10a <- ulam(
  alist(
    T ~ dpois( lambda ),
    #log(lambda) <- a[cid] + b[cid]*P,
    log(lambda) <- a[cid]+log_P*(b[cid]-g),
    a[cid] ~ dnorm(3, 0.5),
    b[cid] ~ dnorm(0, 0.2),
    g ~ dnorm(0, 0.2)
  ), data=dat, chains=4, cores=4, log_lik = TRUE
)

compare(m11.9, m11.10, func=PSIS)

k <- PSIS( m11.10, pointwise = TRUE)$k
plot( dat$P, dat$T, xlab='log Population (std)', ylab='total tools',
      col=rangi2, pch=ifelse(dat$cid==1,1,16), lwd=2,
      ylim=c(0,120), xlim=c(0, 275000),cex=1+normalize(k))

ns <- 100
P_seq <- seq( from=0, to=275000, length.out = ns)

lambda <- link(m11.10, data=data.frame(P=P_seq, cid=1 ) )
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines( P_seq, lmu, lty=2, lwd=1.5)
shade(lci, P_seq, xpd=TRUE)

lambda <- link(m11.10, data=data.frame(P=P_seq, cid=2 ) )
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines( P_seq, lmu, lty=1, lwd=1.5)
shade(lci, P_seq, xpd=TRUE)


plot( dat$P, dat$T, xlab='log Population (std)', ylab='total tools',
      col=rangi2, pch=ifelse(dat$cid==1,1,16), lwd=2,
      ylim=c(0,120), xlim=c(0, 275000),cex=1+normalize(k))

ns <- 100
P_seq <- seq( from=0, to=275000, length.out = ns)

lambda <- link(m11.10b, data=data.frame(P=P_seq, cid=1 ) )
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines( P_seq, lmu, lty=2, lwd=1.5)
shade(lci, P_seq, xpd=TRUE)

lambda <- link(m11.10b, data=data.frame(P=P_seq, cid=2 ) )
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines( P_seq, lmu, lty=1, lwd=1.5)
shade(lci, P_seq, xpd=TRUE)



plot( dat$log_P, dat$T, xlab='log Population (std)', ylab='total tools',
      col=rangi2, pch=ifelse(dat$cid==1,1,16), lwd=2,
      ylim=c(0,120), xlim=c(-1, 3),cex=1+normalize(k))

ns <- 100
P_seq <- seq( from=-1, to=3, length.out = ns)


lambda <- link(m11.10a, data=data.frame(log_P=P_seq, cid=1 ) )
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines( P_seq, lmu, lty=2, lwd=1.5)
shade(lci, P_seq, xpd=TRUE)

lambda <- link(m11.10a, data=data.frame(log_P=P_seq, cid=2 ) )
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines( P_seq, lmu, lty=1, lwd=1.5)
shade(lci, P_seq, xpd=TRUE)


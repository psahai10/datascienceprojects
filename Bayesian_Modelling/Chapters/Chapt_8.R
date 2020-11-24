library(rethinking)
data(rugged)
d <- rugged

head(d)
library(psych)
psych::describe(d)

d$log_gdp <- log(d$rgdppc_2000 )

dd <- d[complete.cases(d$log_gdp), ]

dd$log_gdp_std <- dd$log_gdp/mean(dd$log_gdp)
dd$rugged_std <- dd$rugged/max(dd$rugged)

library(MASS)
dd$boxcox_rugged <- boxcox(dd$rugged)

mean(dd$rugged_std)

m8.1 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a + b*(rugged_std-0.215),
    a ~ dnorm(1, 1),
    b ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=dd
)

set.seed(7)
prior <- extract.prior(m8.1)

plot(NULL, xlim=c(0,1), ylim=c(0.5, 1.5),
     xlab='ruggedness', ylab='log GDP' )
abline( h=c(min(dd$log_gdp_std), max(dd$log_gdp_std)), col=c('black', 'red'), lty=c(2, 2))

rugged_seq <- seq(from=-0.1, to=1.1, length.out =30)
mu <- link(m8.1, post=prior, data=data.frame(rugged_std=rugged_seq))
for (i in 1:50) lines(rugged_seq, mu[i,] , col=col.alpha('black', 0.3))

sum(abs(prior$b) > 0.6) / length(prior$b)



m8.1b <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a + b*(rugged_std-0.215),
    a ~ dnorm(1, 0.1),
    b ~ dnorm(0, 0.4),
    sigma ~ dexp(1)
  ), data=dd
)

set.seed(7)
prior <- extract.prior(m8.1b)

plot(NULL, xlim=c(0,1), ylim=c(0.5, 1.5),
     xlab='ruggedness', ylab='log GDP' )
abline( h=c(min(dd$log_gdp_std), max(dd$log_gdp_std)), col=c('black', 'red'), lty=c(2, 2))

rugged_seq <- seq(from=-0.1, to=1.1, length.out =30)
mu <- link(m8.1b, post=prior, data=data.frame(rugged_std=rugged_seq))
for (i in 1:50) lines(rugged_seq, mu[i,] , col=col.alpha('black', 0.3))

sum(abs(prior$b) > 0.6) / length(prior$b)

precis(m8.1b)

dd$cid <- ifelse(dd$cont_africa==1, 1, 2)

m8.2 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b*(rugged_std-0.215),
    a[cid] ~ dnorm(1, 0.1),
    b ~ dnorm(0, 0.4),
    sigma ~ dexp(1)
  ), data=dd
)
compare(m8.1b, m8.2)

precis(m8.2, depth=2)

post <- extract.samples(m8.2)
str(post)

diff_a1_a2 <- post$a[,1] - post$a[,2] 
PI( diff_a1_a2 )

rugged_seq <- seq(from=-0.1, to=1.1, length.out =30)
mu.NotAfrica <- link(m8.2, data=data.frame(cid=2, rugged_std=rugged_seq))
mu.Africa <- link(m8.2, data=data.frame(cid=1, rugged_std=rugged_seq))

mu.NotAfrica_mu <- apply(mu.NotAfrica, 2, mean)
mu.NotAfrica_ci <- apply(mu.NotAfrica, 2, PI, prob=0.97)

mu.Africa_mu <- apply(mu.Africa, 2, mean)
mu.Africa_ci <- apply(mu.Africa, 2, PI, prob=0.97)

plot(log_gdp_std ~ rugged_std, data= dd, col=col.alpha(rangi2, 0.8))


lines(rugged_seq, mu.NotAfrica_mu)
shade(mu.NotAfrica_ci, rugged_seq)

lines(rugged_seq, mu.Africa_mu)
shade(mu.Africa_ci, rugged_seq)

m8.3 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1, 0.1),
    b[cid] ~ dnorm(0, 0.4),
    sigma ~ dexp(1)
  ), data=dd
)

precis(m8.3, depth=2)

compare(m8.1b, m8.2, m8.3, func=PSIS)
plot(compare(m8.1b, m8.2, m8.3, func=PSIS))

waic_m8.3 <- WAIC(m8.3, pointwise = TRUE)$WAIC
waic_m8.2 <- WAIC(m8.2, pointwise = TRUE)$WAIC
n <- length(waic_m8.3)
diff_m6.7_m6.8 <- waic_m8.3 - waic_m8.2
sqrt(n*var(diff_m6.7_m6.8))

PSIS_m8.3 <- PSIS(m8.3, pointwise = TRUE)
WAIC_m8.3 <- WAIC(m8.3, pointwise = TRUE)
plot(PSIS_m8.3$k, WAIC_m8.3$penalty, xlab = 'PSIS Pareto k',
     ylab='WAIC penalty', col=rangi2, lwd=2)


7 + c(-1, 1)*6.75*2.6


plot( PSIS(m8.3, pointwise=TRUE)$k)

m8.1t <- quap(
  alist(
    log_gdp_std ~ dstudent(2, mu, sigma),
    mu <- a + b*(rugged_std-0.215),
    a ~ dnorm(1, 0.1),
    b ~ dnorm(0, 0.4),
    sigma ~ dexp(1)
  ), data=dd
)


m8.2t <- quap(
  alist(
    log_gdp_std ~ dstudent(2, mu, sigma),
    mu <- a[cid] + b*(rugged_std-0.215),
    a[cid] ~ dnorm(1, 0.1),
    b ~ dnorm(0, 0.4),
    sigma ~ dexp(1)
  ), data=dd
)

m8.3t <- quap(
  alist(
    log_gdp_std ~ dstudent(2, mu, sigma),
    mu <- a[cid] + b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1, 0.1),
    b[cid] ~ dnorm(0, 0.4),
    sigma ~ dexp(1)
  ), data=dd
)

precis(m8.3t, depth=2)

compare(m8.1t, m8.2t, m8.3t, func=PSIS)
plot(compare(m8.1t, m8.2t, m8.3t, func=PSIS))

plot( PSIS(m8.3, pointwise=TRUE)$k)
precis(m8.3t, depth=2)

compare(m8.3, m8.3t, func=PSIS)

plot(compare(m8.3, m8.3t, func=PSIS))

plot( PSIS(m8.3t, pointwise=TRUE)$k)

PSIS_m8.3t <- PSIS(m8.3t, pointwise = TRUE)
WAIC_m8.3t <- WAIC(m8.3t, pointwise = TRUE)
plot(PSIS_m8.3t$k, WAIC_m8.3t$penalty, xlab = 'PSIS Pareto k',
     ylab='WAIC penalty', col=rangi2, lwd=2)

waic_m8.3t <- WAIC(m8.3t, pointwise = TRUE)$WAIC
waic_m8.2t <- WAIC(m8.2t, pointwise = TRUE)$WAIC
n <- length(waic_m8.1t)
diff_m6.7_m6.8 <- waic_m8.3t - waic_m8.2t
sqrt(n*var(diff_m6.7_m6.8))

7 + c(-1, 1)*7.44*2.6


# plot Africa- cid=1
d.A1 <- dd[dd$cid==1, ]
plot(d.A1$rugged_std, d.A1$log_gdp_std, pch=16, col=rangi2, xlim=c(0,1),
     xlab='ruggedness (std)', ylab='log GDP' )
rugged_seq <- seq(from=-0.1, to=1.1, length.out =30)
mu <- link(m8.3, data=data.frame(cid=1, rugged_std=rugged_seq))
mu_mean <- apply(mu, 2, mean)
mu_ci <- apply(mu, 2, PI, prob=0.97)
lines(rugged_seq, mu_mean, lwd=2)
shade(mu_ci, rugged_seq, col=col.alpha(rangi2, 0.3))
mtext('African nations')

# plot Africa- cid=1
d.A0 <- dd[dd$cid==2, ]
plot(d.A0$rugged_std, d.A0$log_gdp_std, pch=1, col='black', xlim=c(0,1),
     xlab='ruggedness (std)', ylab='log GDP' )
rugged_seq <- seq(from=-0.1, to=1.1, length.out =30)
mu <- link(m8.3, data=data.frame(cid=2, rugged_std=rugged_seq))
mu_mean <- apply(mu, 2, mean)
mu_ci <- apply(mu, 2, PI, prob=0.97)
lines(rugged_seq, mu_mean, lwd=2)
shade(mu_ci, rugged_seq, col=col.alpha(rangi2, 0.3))
mtext('Non-African nations')

# plot Africa- cid=1
d.A1 <- dd[dd$cid==1, ]
plot(d.A1$rugged_std, d.A1$log_gdp_std, pch=16, col=rangi2, xlim=c(0,1),
     xlab='ruggedness (std)', ylab='log GDP' )
rugged_seq <- seq(from=-0.1, to=1.1, length.out =30)
mu <- link(m8.3t, data=data.frame(cid=1, rugged_std=rugged_seq))
mu_mean <- apply(mu, 2, mean)
mu_ci <- apply(mu, 2, PI, prob=0.97)
lines(rugged_seq, mu_mean, lwd=2)
shade(mu_ci, rugged_seq, col=col.alpha(rangi2, 0.3))
mtext('African nations')

# plot Africa- cid=1
d.A0 <- dd[dd$cid==2, ]
plot(d.A0$rugged_std, d.A0$log_gdp_std, pch=1, col='black', xlim=c(0,1),
     xlab='ruggedness (std)', ylab='log GDP' )
rugged_seq <- seq(from=-0.1, to=1.1, length.out =30)
mu <- link(m8.3t, data=data.frame(cid=2, rugged_std=rugged_seq))
mu_mean <- apply(mu, 2, mean)
mu_ci <- apply(mu, 2, PI, prob=0.97)
lines(rugged_seq, mu_mean, lwd=2)
shade(mu_ci, rugged_seq, col=col.alpha(rangi2, 0.3))
mtext('Non-African nations')


sim_dat <- data.frame(rugged_std=rugged_seq)

rugged_seq <- seq(from=0, to=1, length.out = 30)
muA <- link(m8.3, data=data.frame(cid=1, rugged_std=rugged_seq))
muA.mean <- apply(muA, 2, mean)
mu.PIA <- apply(muA, 2, PI)
plot( log_gdp_std ~ rugged_std , data=dd )
lines( rugged_seq , muA.mean , lwd=2 )
shade( mu.PIA , rugged_seq )


muN <- link(m8.3, data=data.frame(cid=2, rugged_std=rugged_seq))
muN.mean <- apply(muN, 2, mean)
mu.PIN <- apply(muN, 2, PI)
plot( log_gdp_std ~ rugged_std , data=dd )
lines( rugged_seq , muN.mean , lwd=2 )
shade( mu.PIN , rugged_seq )


mu.mean <- muA.mean - muN.mean
mu.PI <- mu.PIA - mu.PIN
plot( log_gdp_std ~ rugged_std , data=dd )
lines(rugged_seq, mu.mean , lwd=2 )
shade( mu.PI , rugged_seq )

data(tulips)
d <- tulips
str(d)

d$blooms_std <- d$blooms / max(d$blooms)
d$water_cent <- d$water - mean(d$water)
d$shade_cent <- d$shade - mean(d$shade)
d$bid <- ifelse(d$bed=='a', 1, ifelse(d$bed=='b', 2, 3))



a <- rnorm(1e4, 0.5, 1); sum(a < 0 | a > 1)/length(a)

a <- rnorm(1e4, 0.5, 0.25); sum(a < 0 | a > 1)/length(a)
library(rethinking)
m8.4 <- quap(
  alist(
    blooms_std ~ dnorm(mu, sigma),
    mu <- a + bw*water_cent + bs*shade_cent,
    a ~ dnorm(0.5, 0.25),
    bw ~ dnorm(0, 0.25),
    bs ~ dnorm(0, 0.25),
    sigma ~ dexp(1)
  ), data= d
)

set.seed(7)
prior <- extract.prior(m8.4)
plot(NULL, xlim=c(-1.5,1), ylim=c(-0.5, 0.5),
     xlab='Water and Shade', ylab='Bloom' )

bloom_seq <- seq(from=-1.5, to=1.1, length.out=30)
mu <- link(m8.4, post=prior, data=data.frame(blooms_std=bloom_seq))
for (i in 1:50) lines(bloom_seq, mu[i,] , col=col.alpha('black', 0.3))


m8.5 <- quap(
  alist(
    blooms_std ~ dnorm(mu, sigma),
    mu <- a + bw*water_cent + bs*shade_cent + bws*water_cent*shade_cent,
    a ~ dnorm(0.5, 0.25),
    bw ~ dnorm(0, 0.25),
    bs ~ dnorm(0, 0.25),
    bws ~ dnorm(0, 0.25),
    sigma ~ dexp(1)
  ), data= d
)
precis(m8.4)
precis(m8.5)

m8.5b <- quap(
  alist(
    blooms_std ~ dnorm(mu, sigma),
    mu <- a[bid] + bw[bid]*water_cent + bs[bid]*shade_cent + bws[bid]*water_cent*shade_cent,
    a[bid] ~ dnorm(0.5, 0.25),
    bw[bid] ~ dnorm(0, 0.25),
    bs[bid] ~ dnorm(0, 0.25),
    bws[bid] ~ dnorm(0, 0.25),
    sigma ~ dexp(1)
  ), data= d
)
precis(m8.4)
precis(m8.5)
precis(m8.5b, depth=2)

par(mfrow=c(1,3))

for (s in -1:1) {
  idx <- which(d$shade_cent==s)
  plot( d$water_cent[idx], d$blooms_std[idx], xlim=c(-1,1), ylim=c(0,1),
        xlab='water', ylab='blooms', pch=16, col=rangi2)
  mu <- link(m8.4, data=data.frame( shade_cent=s, water_cent=-1:1))
  
  for (i in 1:20) lines(-1:1, mu[i,], col=col.alpha('black', 0.3))
}


par(mfrow=c(1,3))

for (s in -1:1) {
  idx <- which(d$shade_cent==s)
  plot( d$water_cent[idx], d$blooms_std[idx], xlim=c(-1,1), ylim=c(0,1),
        xlab='water', ylab='blooms', pch=16, col=rangi2)
  mu <- link(m8.5, data=data.frame( shade_cent=s, water_cent=-1:1))
  
  for (i in 1:20) lines(-1:1, mu[i,], col=col.alpha('black', 0.3))
}



set.seed(7)
prior <- extract.prior(m8.5)

plot(NULL, xlim=c(0,1), ylim=c(0.5, 1.5),
     xlab='water', ylab='bloom' )
bloom_seq <- seq(from=0, to=1, length.out =30)
mu <- link(m8.1, post=prior, data=data.frame(bloom_std=bloom_seq))
for (i in 1:50) lines(bloom_seq, mu[i,] , col=col.alpha('black', 0.3))


# plot bid =1
d.1 <- d[d$bid==1, ]
plot(d.1$water_cent, d.1$blooms_std, pch=16, col=rangi2, xlim=c(-1,1),
     xlab='water (std)', ylab='Blooms' )
shade_seq <- seq(from=-1.1, to=1.1, length.out =30)
water_seq <- seq(from=-1.1, to=1.1, length.out =30)
mu <- link(m8.4, data=data.frame(bid=1, water_cent=water_seq, shade_cent=shade_seq))
mu_mean <- apply(mu, 2, mean)
mu_ci <- apply(mu, 2, PI, prob=0.97)
lines(water_seq, mu_mean, lwd=2)
shade(mu_ci, water_seq, col=col.alpha(rangi2, 0.3))
mtext('Bed A')

# plot bid=1
d.2 <- d[d$bid==2, ]
plot(d.2$blooms_std, d.2$water_cent, pch=16, col=rangi2, xlim=c(-1,1),
     xlab='water (std)', ylab='Blooms' )
shade_seq <- seq(from=-1.1, to=1.1, length.out =30)
water_seq <- seq(from=-1.1, to=1.1, length.out =30)
mu <- link(m8.4, data=data.frame(bid=2, water_cent=water_seq, shade_cent=shade_seq))
mu_mean <- apply(mu, 2, mean)
mu_ci <- apply(mu, 2, PI, prob=0.97)
lines(water_seq, mu_mean, lwd=2)
shade(mu_ci, water_seq, col=col.alpha(rangi2, 0.3))
mtext('Bed B')

# plot Africa- cid=1
d.3 <- d[d$bid==3, ]
plot(d.3$blooms_std, d.3$water_cent, pch=16, col=rangi2, xlim=c(-1,1),
     xlab='water (std)', ylab='Blooms' )
shade_seq <- seq(from=-1.1, to=1.1, length.out =30)
water_seq <- seq(from=-1.1, to=1.1, length.out =30)
mu <- link(m8.4, data=data.frame(bid=3, water_cent=water_seq, shade_cent=shade_seq))
mu_mean <- apply(mu, 2, mean)
mu_ci <- apply(mu, 2, PI, prob=0.97)
lines(water_seq, mu_mean, lwd=2)
shade(mu_ci, water_seq, col=col.alpha(rangi2, 0.3))
mtext('Bed C')

library(rethinking)
data("Primates301")
d <- Primates301
head(d)

d <- d[complete.cases(d$brain, d$body, d$longevity), ] 

d$L <- log(d$longevity)
d$M <- log(d$body)
d$B <- log(d$brain)

m



library(rethinking)
library(ggplot2)

p_grid <- seq(from=0, to=1, length.out=10)
prior <- rep(1,10)
prior <- ifelse(p_grid<0.5,0,1)
prior <- exp(-5*abs(p_grid-0.5))
likelihood <- dbinom(6, 9, prob=p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)
plot(p_grid, posterior, type='b',
     xlab='probability of water', ylab='posterior probability')
mtext('20 points')

### Quadratic Functions

globe.qa <- quap(
  alist(
    W ~ dbinom( W+L, p ) ,
    p ~ dunif(0,1)
  ),
    data=list(W=6, L=3))
precis(globe.qa)

W<-24
L<-12
curve(dbeta(x, W+1, L+1), from=0, to=1)
curve(dnorm(x, 0.67, 0.16), lty=2, add=TRUE)

## FIRST MCMC
n_samples <- 1000
p <- rep(NA, n_samples)
p[1] <-0.5
W <- 6
L <- 3
for (i in 2:n_samples) {
  p_new <-rnorm(1, p[i-1], 0.1)
  if (p_new < 0) p_new <- abs(p_new)
  if (p_new > 1) p_new <- 2-p_new
  q0 <- dbinom(W, W+L, p[i-1])
  q1 <- dbinom(W, W+L, p_new)
  p[i] <- ifelse(runif(1) < q1/q0, p_new, p[i-1])
}
dens(p, xlim=c(0,1))
curve(dbeta(x, W+1, L+1), lty=2, add=TRUE)

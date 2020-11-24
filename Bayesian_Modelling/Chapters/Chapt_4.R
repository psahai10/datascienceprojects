library(rethinking)
pos <- replicate(1000, sum( runif(1e5,-1, 1)))
hist(pos)
plot(density(pos))

var <- prod(1+runif(12,0,0.1))
growth <-replicate(1e5, var)
dens(growth, norm.comp = TRUE)

w <- 6; n <- 9
p_grid <- seq(0,1,length.out = 100)
posterior <- dbinom(w, n, p_grid)*dunif(p_grid,0,1)
posterior <- posterior/sum(posterior)
posterior

data(Howell1)
d <- Howell1
head(d)
str(d)
library(psych)
psych::describe(d)
hist(d)

d2 <- d[ d$age>=18, ]

dens(d2$height)
curve(dnorm(x, 178, 20), from=100, to=250)
curve(dunif(x, 0, 50), from=-10, to=60)



sample_mu <- rnorm(1e4, 178, 20)
sample_sigma <- runif(1e4, 0, 50)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)

####
# Step1: generate a sequence of mu and sigma
mu.list <- seq(150, 160, length.out = 100)
sigma.list <- seq(7, 9, length.out = 100)

# Step2: Combine this sequence in a grid that will equal length = 100X100
post <- expand.grid(mu=mu.list, sigma=sigma.list)

# Create a function that takes the height from the data, runs a dnorm
# simulation for all height with respective mu and sigma and sums it up
# and return it as sum total
func <- function(x) { sum(
    dnorm( d2$height, post$mu[x], post$sigma[x], log=TRUE))
}

# Uses sapply and runs this sigma, mu calculations for all combinations
# in mu and sigma in the post matrix
post$LL <- sapply(1:nrow(post), func)

# Sums post$LL with a simulation for pnorm, with all mu, and 178,20 paramter 
# and sums it for with a dunif simulation
post$prod <- post$LL + dnorm( post$mu, 178, 20, TRUE) +
  dunif(post$sigma, 0, 50, TRUE)

# runs an exponential for the diference between the product of the combination
# with the lowest value which will equal the probability 
post$prob <- exp(post$prod - max(post$prod))

contour_xyz(post$mu, post$sigma, post$prob)
image_xyz(post$mu, post$sigma, post$prob)



sample.rows <- sample(1:nrow(post), size=1e4, replace=TRUE, prob=post$prob)
sample.mu <- post$mu[ sample.rows]
sample.sigma <- post$sigma[ sample.rows]

plot(sample.mu, sample.sigma, cex=0.5, pch=16, col=col.alpha(rangi2,0.3))

dens(sample.mu)
dens(sample.sigma)

PI(sample.mu)
PI(sample.sigma)

###############
# Step1:
d3 <- sample(d2$height, size=20)
mu.list2 <- seq(150, 170, length.out = 200)
sigma.list2 <- seq(4, 20, length.out = 200)

# Step2:
post2 <- expand.grid(mu=mu.list2, sigma=sigma.list2)

func <- function(x) { sum(
  dnorm( d3, mean=post2$mu[x], sd=post2$sigma[x], log=TRUE))
}

# Step3:
post2$LL <- sapply(1:nrow(post2), func)


# Step4:
post2$prod <- post2$LL + dnorm( post2$mu, 178, 20, TRUE) +
  dunif(post2$sigma, 0, 50, TRUE)


# Step5:
post2$prob <- exp(post2$prod - max(post2$prod))

contour_xyz(post2$mu, post2$sigma, post2$prob)
image_xyz(post2$mu, post2$sigma, post2$prob)
sample2.rows <- sample(1:nrow(post2), size=1e4, replace=TRUE, prob=post2$prob)
sample2.mu <- post2$mu[ sample2.rows]
sample2.sigma <- post2$sigma[ sample2.rows]

plot(sample2.mu, sample2.sigma, cex=0.5, col=col.alpha(rangi2,0.3),
     xlab='mu', ylab='sigma', pch=16)

dens(sample2.mu)
dens(sample2.sigma, norm.comp = TRUE)

PI(sample2.mu)
PI(sample2.sigma)

data("Howell1")
d <- Howell1
d2 <- d[d$age>=18,]

flist <- alist(
  height ~ dnorm(mu, sigma),
  mu ~ dnorm(178, 20),
  sigma ~ dunif(0, 50)
)
m4.1 <- quap(flist, data=d2)
precis(m4.1)
precis(m4.1, prob=.99)

start <- list(
  mu=mean(d2$height),
  sigma=sd(d2$height)
)
m4.1 <- quap(flist, data=d2, start=start)

m4.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu ~ dnorm(178, 0.1),
    sigma ~ dunif(0, 50)
  ) , data=d2
)
precis(m4.2)
precis(m4.2, prob=0.99)

vcov(m4.1)
diag( vcov(m4.1))
cov2cor(vcov(m4.1))

library(rethinking)
post <- extract.samples(m4.1, n=1e4)
precis(post, hist=FALSE)

library(MASS)
post <- mvrnorm(1e4, mu=coef(m4.1), Sigma=vcov(m4.1))
plot(d2$height~d2$weight)

set.seed(2971)
N <- 100
a <- rnorm(N, 178, 20)
b <- rnorm(N, 1, 10)

plot( NULL, xlim=range(d2$weight), ylim=c(-100, 400),
     xlab='weight', ylab='height')
abline(h=0, lty=2)
abline(h=272, lty=1, lwd=0.5)
mtext('b ~ dnorm(0, 10)')
xbar <- mean(d2$weight)
for (i in 1:N) curve(a[i] + b[i]*(x-xbar),
    from=min(d2$weight), to=max(d2$weight), add=TRUE,
    col=col.alpha('black', 0.2) )

b <-rlnorm(1e4, 0, 1)
dens(b, xlim=c(0,10), adj=0.1)

xbar <- mean(d2$weight)
m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight - xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ) , data=d2)

xbar <- mean(d2$weight)

m4.3b <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + exp(log_b)*(weight - xbar),
    a ~ dnorm(178, 20),
    log_b ~ dnorm(0, 1),
    sigma ~ dunif(0, 50)
  ) , data=d2)

precis(m4.3)

round(vcov(m4.3), 3)

plot(height ~ weight, data=d2, col=rangi2)
post <- extract.samples(m4.3)
a_map <- mean(post$a)
b_map <- mean(post$b)
curve(a_map, b_map*(x-xbar), add=TRUE)

post <- extract.samples(m4.3)
post[1:5,]

N <- 100
dN <- d2[1:N, ]
mN <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight-mean(weight)),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0, 50)
  ), data=d2
)

post <- extract.samples(mN, n=100)
plot(dN$weight, dN$height,
     xlim=range(d2$weight), ylim=range(140, 170), 
     col=rangi2, xlab='weight', ylab='height')
mtext(concat("N= ", N))

for (i in 1:100) curve(post$a[i] + post$b[i]*(x-mean(dN$weight)), 
        col=col.alpha('black', 0.1), add=TRUE)

post <- extract.samples(m4.3)
mu_at_50 <- post$a + post$b * (50-xbar)

dens(mu_at_50, col=rangi2, lwd=2, xlab='mu|weight50')
PI(mu_at_50, prob=0.89)

mu<-link(m4.3)
str(mu)

weight.seq <- seq(from=20, to=75, by=1)
mu <- link(m4.3, data=data.frame(weight=weight.seq))
str(mu)
plot(height~weight, d2, type='n')
for (i in 1:100) points(weight.seq, mu[i,], pch=16, col=col.alpha(rangi2,0.1))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)

plot(height~weight, data=d2, col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)


post <- extract.samples(m4.3)
mu.link <- function(weight) post$a + post$b*(weight-xbar)
weight.seq <- seq(from=25, to=70, by=1)
mu <- sapply(weight.seq, mu.link)
mu.mean <- apply(mu,2,mean)
mu.CI <- apply(mu, 2, PI, prob=0.89)

sim.height <- sim(m4.3, data=list(weight=weight.seq))
sim.height <- sim(m4.3, data=list(weight=weight.seq), n=1e4)
str(sim.height)
height.PI <- apply(sim.height, 2, PI, prob=0.67)
height.PI

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)

plot(height~weight, data=d2, col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)


post <- extract.samples(m4.3)
weight.seq <- 25:70
func <- function(x) {
  rnorm(n=row(post), mean=post$a+post$b*(x-xbar), sd=post$sigma)
}
sim.height <- sapply(weight.seq, func)
height.PI <- apply(sim.height, 2, PI, prob=0.89)


d<-Howell1
plot(height~weight, d)

d$weight_s <- (d$weight-mean(d$weight))/sd(d$weight)
d$weight_s2 <- d$weight_s^2

m4.5 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2,
    a ~ dnorm(178, 20),
    b1 ~ dlnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dunif(0, 50)
    ) , data=d)
precis(m4.5)

weight.seq <- seq(from=-2.2, to=2, length.out=30)
pred_dat <- list(weight_s=weight.seq, weight_s2=weight.seq^2)
mu <- link(m4.5, data=pred_dat)
mu.PI <- apply(mu, 2, PI, prob=0.89)
sim.height <- sim(m4.5, data=pred_dat)
height.PI <- apply(sim.height, 2, PI, prob=0.89)

plot(height~weight_s, data=d, col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)

d$weight_s3 <- d$weight_s^3
m4.6 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2 + b3*weight_s3,
    a ~ dnorm(178, 20),
    b1 ~ dlnorm(0, 1),
    b2 ~ dnorm(0, 10),
    b3 ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ) , data=d)
precis(m4.6)

weight.seq <- seq(from=-2.2, to=2, length.out=30)
pred_dat <- list(weight_s=weight.seq, weight_s2=weight.seq^2, weight_s3=weight.seq^3)
mu <- link(m4.6, data=pred_dat)
mu.PI <- apply(mu, 2, PI, prob=0.89)
sim.height <- sim(m4.6, data=pred_dat)
height.PI <- apply(sim.height, 2, PI, prob=0.89)

mu.mean <- apply(mu, 2, mean)

plot(height~weight_s, data=d, col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)

plot(height~weight_s, data=d, col=col.alpha(rangi2, 0.5), xaxt='n')
at <- c(-2,-1,0,1,2)
labels <- at*sd(d$weight) + mean(d$weight)
axis(side=1, at=at, labels=round(labels,1))


data(cherry_blossoms)
d <- cherry_blossoms
precis(d, hist=FALSE)
names(d)
plot(doy~year, d, col=col.alpha(rangi2,0.5))


d2 <- d[ complete.cases(d$doy), ]
num_knots <- 15
knot_list <- quantile(d2$year, probs=seq(0, 1,length.out=15))
knot_list
library(splines)
B <- bs(d2$year, 
        knots=knot_list[-c(1, num_knots)],
        degree=3, intercept=TRUE)

plot(NULL, xlim=range(d2$year), ylim=c(0,1), xlab='year', ylab='basis')
for (i in 1:ncol(B) ) lines(d2$year, B[,i])

library(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce
names(d)
library(psych)
psych::describe(d)
plot(d)

d$D <-standardize(d$Divorce)
d$M <-standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

sd(d$MedianAgeMarriage)

m5.1 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data= d
)

set.seed(10)
prior <- extract.prior(m5.1)
mu <- link(m5.1, post=prior, data=list(A=c(-2,2)))
plot(NULL, xlim=c(-2,2), ylim=c(-2,2))
for (i in 1:50) lines(c(-2,2), mu[i,], col=col.alpha('black', 0.4)) 

round(vcov(m5.1), 3)

#Get the uncertainty for every Z-point in Age of Marriage with divorce rate
A_seq <- seq(from=-3, to=3.2, length.out = 30)

mu <- link(m5.1, data=list(A=A_seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(D~A, data=d, col=rangi2)
lines(A_seq, mu.mean, lwd=2)
shade(mu.PI, A_seq)
precis(m5.1)

######################## AGE of marriage and Marriage Rate ####################

m5.1b <- quap(
  alist(
    M ~ dnorm(mu, sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data= d
)
set.seed(10)
prior <- extract.prior(m5.1b)
mu <- link(m5.1b, post=prior, data=list(A=c(-2,2)))
plot(NULL, xlim=c(-2,2), ylim=c(-2,2))
for (i in 1:50) lines(c(-2,2), mu[i,], col=col.alpha('black', 0.4)) 

A_seq <- seq(from=-3, to=3.2, length.out = 30)
mu <- link(m5.1b, data=list(A=A_seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(M~A, data=d, col=rangi2)
lines(A_seq, mu.mean, lwd=2)
shade(mu.PI, A_seq)
precis(m5.1b)

##############################################################################


m5.2 <- quap(
    alist(
      D ~ dnorm(mu, sigma),
      mu <- a + bM * M,
      a ~ dnorm(0, 0.2),
      bM ~ dnorm(0, 0.5),
      sigma ~ dexp(1)
    ), data= d
  )

set.seed(10)
prior <- extract.prior(m5.2)
mu <- link(m5.2, post=prior, data=list(M=c(-2,2)))
plot(NULL, xlim=c(-2,2), ylim=c(-2,2))
for (i in 1:50) lines(c(-2,2), mu[i,], col=col.alpha('black', 0.4)) 
  
M_seq <- seq(from=-3, to=3.2, length.out = 30)
mu <- link(m5.2, data=list(M=A_seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
  
plot(D~M, data=d, col=rangi2)
lines(M_seq, mu.mean, lwd=2)
shade(mu.PI, A_seq)
precis(m5.2)


DMA_dag1 <- dagitty('dag{ D <- A -> M -> D}')
impliedConditionalIndependencies(DMA_dag1)

DMA_dag2 <- dagitty('dag{ D <- A -> M}')
impliedConditionalIndependencies(DMA_dag2)


m5.3 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M + bA * A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data= d
)
precis(m5.3)

compare(m5.1, m5.2, m5.3)

plot(coeftab(m5.1, m5.1b, m5.2, m5.3), par=c('bA', 'bM'))

m5.4 <- quap(
  alist(
    M ~ dnorm(mu, sigma),
    mu <- a + bAM * A,
    a ~ dnorm(0, 0.2),
    bAM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
    ), data = d)

mu <- link(m5.4)
mu_mean <- apply(mu, 2, mean)
mu_resid <- d$M - mu_mean 

mu <- link(m5.3)
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)

D_sim <- sim(m5.3, 10000)
D_PI <- apply(sim, 2, PI)

plot(mu_mean ~ d$D, col=rangi2, ylim=range(mu_PI),
     xlab='Observed divorce', ylab='Predictive divorce')
abline(a=0, b=1, lty=2)
for (i in 1:nrow(d)) lines(rep(d$D[i],2), mu_PI[,i], col=rangi2)

########## SIMULATION ################
N <- 200
x_real <- rnorm(N)
x_spur <- rnorm(N, x_real)
y <- rnorm(N, x_real)
d <- data.frame(y, x_real, x_spur)

mS <- quap(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a + bx_real*x_real + bx_spur * x_spur,
    a ~ dnorm(0, 0.2),
    bx_real ~ dnorm(0, 0.5),
    bx_spur ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)
pairs(d)
precis(mS)
plot(coeftab(mS), par=c('bx_real', 'bx_spur'))
########## SIMULATION END ################





data(WaffleDivorce)
d <- list()
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage)

m5.3A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm(mu, sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1),
    ## A -> M
    M ~ dnorm(mu_M, sigma_M),
    mu_M <- aM + bAM*A,
    aM ~ dnorm(0, 0.2),
    bAM ~ dnorm(0, 0.5),
    sigma_M ~ dexp(1)
  ), data=d )

precis(m5.3_A)

A_seq <- seq(from=-2, to=2, length.out = 30)

sim_dat <- data.frame(A=A_seq)

M_sim <- with(post, sapply(1:30, 
      function(x) rnorm(1e3, aM + bAM*A_seq[x], sigma_M)) )

D_sim <- with(post, sapply(1:30, 
      function(x) rnorm(1e3, a + bA*A_seq[x] + bM*M_sim[,x], sigma)) )



plot(sim_dat$A, colMeans(D_sim), ylim=c(-2,2), type='l',
     xlab='manipulated A', ylab='counterfactual D')
shade(apply(D_sim,2,PI), sim_dat$A)
mtext('Total counterfactual effect of A on D')


plot(sim_dat$A, colMeans(M_sim), ylim=c(-2,2), type='l',
     xlab='manipulated A', ylab='counterfactual M')
shade(apply(M_sim,2,PI), sim_dat$A)
mtext('Total counterfactual effect of A on M')



####################### WORKING PROPERLY???????????? #######################

sim_dat <- data.frame( M=seq(from=-2, to=2, length.out = 30), A=0)


s <- with(post, sapply(1:30, 
          function(x) rnorm(1e3, aM + bAM*sim_dat$M[x], sigma)) )

s <- with(post, sapply(1:30, 
          function(x) rnorm(1e3, a + bA*sim_dat$A[x] + bM*M_sim[,x], sigma)) )


plot(sim_dat$M, colMeans(s), ylim=c(-2,2), type='l',
     xlab='manipulated M', ylab='counterfactual D')
shade(apply(s,2,PI), sim_dat$M)
mtext('Total counterfactual effect of M on D')


#####################################################################


# plot(A_seq, mean(D_sim[,2]))

data(milk)
d <- milk
str(d)

d$K <- standardize(d$kcal.per.g)
d$N <- standardize(d$neocortex.perc)
d$M <- standardize(log(d$mass))

dcc <- d[complete.cases(d$K, d$N, d$M) , ]

m5.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N,
    a ~ dnorm(0,0.2),
    bN ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=dcc
)

prior <- extract.prior(m5.5)
xseq <- c(-2,2)
mu <-link(m5.5, post=prior, data=list(N=xseq))
plot(NULL, xlim=xseq, ylim=xseq)
for (i in 1:50) lines(xseq, mu[i,], col=col.alpha('black', 0.3))
precis(m5.5)

xseq <- seq(from=min(dcc$N)-0.15, to=max(dcc$N)+0.15, length.out = 30)
mu <- link(m5.5, data=list(N=xseq))
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu, 2, PI)
plot(K~N, data=dcc)
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)


m5.6 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bM*M,
    a ~ dnorm(0,0.2),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=dcc
)
precis(m5.6)

prior <- extract.prior(m5.6)
xseq <- c(-2,2)
mu <-link(m5.6, post=prior, data=list(M=xseq))
plot(NULL, xlim=xseq, ylim=xseq)
for (i in 1:50) lines(xseq, mu[i,], col=col.alpha('black', 0.3))
precis(m5.6)

xseq <- seq(from=min(dcc$M)-0.15, to=max(dcc$M)+0.15, length.out = 30)
mu <- link(m5.6, data=list(M=xseq))
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu, 2, PI)
plot(K~M, data=dcc)
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)

library(rethinking)
m5.7 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N + bM*M,
    a ~ dnorm(0,0.2),
    bN ~ dnorm(0, 0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=dcc
)
precis(m5.7)

plot( coeftab( m5.5, m5.6, m5.7) , pars=c("bM","bN"))
pairs(~K + M + N, dcc)


xseq <- seq(from=min(dcc$M)-0.15, to=max(dcc$M)+0.15, length.out = 30)
mu <- link(m5.7, data=list(M=xseq, N=0))
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu, 2, PI)
plot(NULL, xlim=range(dcc$M), ylim=range(dcc$K))
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)

# M -> K <- N
# M -> N
n <- 100
M <- rnorm(n)
N <- rnorm(n, M)
K <- rnorm(n, N-M)
d_sim <- data.frame(K=K, N=N, M=M)

# M -> K <- N
# N -> M
n <- 100
N <- rnorm(n)
M <- rnorm(n, N)
K <- rnorm(n, N-M)
dsim2 <- data.frame(K=K, N=N, M=M)

# M -> K <- N
# M <- U -> N
U <- rnorm(n)
N <- rnorm(n, U)
M <- rnorm(n, U)
K <- rnorm(n, N-M)
dsim3 <- data.frame(K=K, N=N, M=M)

dag5.7 <- dagitty('dag{
                  M -> K <- N
                  M -> N}')
coordinates(dag5.7) <- list(x=c(M=0, K=1, N=2), y=c(M=0.5, K=1, N=0.5))
MElist <- equivalentDAGs(dag5.7)
drawdag(MElist)

data("Howell1")
d <- Howell1
str(d)

mu_female <- rnorm(1e4,168, 20)
mu_male <- rnorm(1e4,168, 20) + rnorm(1e4, 10, 10)
precis(data.frame(mu_female, mu_male), hist=FALSE)
d$sex <- ifelse(d$male == 1, 2, 1)
str(d$sex)

m5.8 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a[sex],
    a[sex] ~ dnorm(178, 20),
    sigma ~ dunif(0, 50)
  ), data=d
)
precis(m5.8, depth=2)
library(rethinking)
data(milk)
d <- milk
str(d)
levels(d$clade)
d$clade_id <- as.integer(d$clade)

d$K <- standardize(d$kcal.per.g)
m5.9 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a[clade_id],
    a[clade_id] ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)

labels <- paste('a[', 1:4, ']:', levels(d$clade), sep='')

plot(precis(m5.9, depth = 2, pars='a'), labels=labels, xlab='exp kcal (std)')

set.seed(63)
d$house <- sample(rep(1:4, each=8), size=nrow(d))
d$house

m5.10 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a[clade_id] + h[house],
    a[clade_id] ~ dnorm(0, 0.5),
    h[house] ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data =d
)

precis(m5.10, depth=2)
labels <- paste('house[', 1:4, ']:', levels(c('Gryffindor',)), sep='')
plot(precis(m5.10, depth = 2, pars='a'), labels=labels, xlab='exp kcal (std)')

post <- extract.samples(m5.10)



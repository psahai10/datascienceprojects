set.seed(1914)
N <- 200
p <- 0.1

nw <- rnorm(N)
tw <- rnorm(N)

s <- nw + tw 

q <- quantile(s, 1-p)

selected <- ifelse( s>=q , TRUE , FALSE)

df <- data.frame(newsworthiness=nw, trustworthiness=tw, total_score=s, quantile=q, selected=selected)

cor(tw[selected], nw[selected])
plot(df$newsworthiness, df$trustworthiness)

N <- 100
set.seed(909)

height <- rnorm(N, 10, 2)

leg_prop <- runif(N, 0.4, 0.5)

leg_left <- leg_prop*height + rnorm(N, 0, 0.02)

leg_right <- leg_prop*height + rnorm(N, 0, 0.02)

d <- data.frame(height, leg_left, leg_right)

m6.1 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left + br*leg_right,
    a ~ dnorm(10, 100),
    bl ~ dnorm(2, 10),
    br ~ dnorm(2, 10),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.1)
plot(precis(m6.1))


post <- extract.samples(m6.1)
plot(bl~br, post, col=col.alpha(rangi2, 0.1), pch=16)

sum_blbr <- post$bl + post$br
dens(post$bl, col=rangi2, lwd=2, xlab='bl density')
dens(post$br, col=rangi2, lwd=2, xlab='br density')
dens(sum_blbr, col=rangi2, lwd=2, xlab='sum of bl and br')

m6.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left,
    a ~ dnorm(10, 100),
    bl ~ dnorm(2, 10),
    sigma ~ dexp(1)
  ), data=d
)

post <- extract.samples(m6.2)
dens(post$bl, col=rangi2, lwd=2, xlab='bl density')

m6.2b <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + br*leg_right,
    a ~ dnorm(10, 100),
    br ~ dnorm(2, 10),
    sigma ~ dexp(1)
  ), data=d
)

post <- extract.samples(m6.2b)
dens(post$br, col=rangi2, lwd=2, xlab='br density')

precis(m6.2)
precis(m6.2b)
plot(coeftab(m6.1, m6.2, m6.2b), par=c('bl', 'br'))


library(rethinking)
data(milk)
d <- milk
d$K <- standardize(d$kcal.per.g)
d$F <- standardize(d$perc.fat)
d$L <- standardize(d$perc.lactose)


m6.3 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bF*F,
    a ~ dnorm(0, 0.2),
    bF ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)

m6.4 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bL*L,
    a ~ dnorm(0, 0.2),
    bL ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)

m6.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bL*L + bF*F,
    a ~ dnorm(0, 0.2),
    bL ~ dnorm(0, 0.2),
    bF ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)

precis(m6.3)
precis(m6.4)
precis(m6.5)
plot(coeftab(m6.3, m6.4, m6.5), par=c('bL', 'bF'))

pairs(~ kcal.per.g + perc.fat + perc.lactose, data=d, col=rangi2)

# Sim perc.fat -> kcal.per.g
sim.coll <- function(r=0.9) {
  d$x <- rnorm(nrow(d), mean=r*d$perc.fat,
               sd=sqrt( (1-r^2)*var(d$perc.fat) ))
  m <- lm(kcal.per.g ~ perc.fat + x, data=d)
  sqrt( diag( vcov(m)))[2]
}


rep.sim.coll <- function(r=0.9, n=100) {
  stddev <- replicate(n, sim.coll(r))
  mean(stddev)
}

r.seq <- seq(from=0, to=0.99, by=0.01)
stddev <- sapply(r.seq, function(z) rep.sim.coll(r=z, n=100))
plot(stddev ~ r.seq, typle='l', col=rangi2, lwd=2, xlab='correlation')


# Sim kcal.per.g -> perc.fat 
sim.coll <- function(r=0.9) {
  d$x <- rnorm(nrow(d), mean=r*d$kcal.per.g,
               sd=sqrt( (1-r^2)*var(d$kcal.per.g) ))
  m <- lm(perc.fat ~ kcal.per.g + x, data=d)
  sqrt( diag( vcov(m)))[2]
}


rep.sim.coll <- function(r=0.9, n=100) {
  stddev <- replicate(n, sim.coll(r))
  mean(stddev)
}

r.seq <- seq(from=0, to=0.99, by=0.01)
stddev <- sapply(r.seq, function(z) rep.sim.coll(r=z, n=100))
plot(stddev ~ r.seq, typle='l', col=rangi2, lwd=2, xlab='correlation')


set.seed(71)
N <- 100
h0 <- rnorm(N, 10, 2)

treatment <- rep(0:1, each=N/2)
fungus <- rbinom(N, size=1, prob=0.5-treatment*0.4)
h1 <- h0 + rnorm(N, 5-3*fungus)

d <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)
precis(d, hist=FALSE)

sim_p <- rlnorm(1e4, 0, 0.25)
precis(data.frame(sim_p))

m6.6 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p ~ dlnorm(0, 0.25),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.6)

m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    bf ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.7)

m6.8 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p <- a + bt*treatment,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.8)

plot(coeftab(m6.7, m6.8), par=c('bt', 'bf'))

library(dagitty)
plant_dag <- dagitty('dag {
                     H_0 -> H_1
                     F -> H_1
                     T -> F}
                     ')

coordinates(plant_dag) <- list(x=c(H_0=0, T=2, F=1.5, H_1=1),
                               y=c(H_0=0, T=0, F=0, H_1=0))

drawdag(plant_dag)

impliedConditionalIndependencies(plant_dag)

set.seed(71)
N <- 1000
h0 <- rnorm(N, 10, 2)
treatment <- rep(0:1, each=N/2)

# M is the hiddent Variable: Moisture
M <- rbern(N)
fungus <- rbinom(N, size=1, prob=0.5 - treatment*0.4 + 0.4*M)
h1 <- h0 + rnorm(N, 0, 5+3*M)
d2 <-data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)

m6.6b <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p <- a + bf*fungus,
    a ~ dlnorm(0, 0.2),
    bf ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d2
)
precis(m6.6b)

m6.7b <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    bf ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d2
)
precis(m6.7b)

m6.8b <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p <- a + bt*treatment,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d2
)
precis(m6.8b)

plot(coeftab(m6.6b, m6.7b, m6.8b), par=c('bt', 'bf'))

library(dagitty)
plant_dag_moisture <- dagitty('dag {
                     H_0 -> H_1
                     M -> F
                     M -> H_1
                     T -> F}
                     ')

coordinates(plant_dag_moisture) <- list(x=c(H_0=0,H_1=0.75, M=1.125, F=1.5, T=2.25),
                               y=c(H_0=0.5, T=0.5, M=0.6, F=0.5, H_1=0.5))
drawdag(plant_dag_moisture)

impliedConditionalIndependencies(plant_dag_moisture)

library(rethinking)

d <- sim_happiness(seed=1977, N_years=1e3)
precis(d, hist=FALSE)

d2 <- d[d$age>=18, ]
d2$A <- (d2$age-18)/(65-18)
d2$mid <- d2$married + 1

m6.9 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a[mid] + bA*A,
    a[mid] ~ dnorm(0, 1),
    bA ~ dnorm(0, 2),
    sigma ~ dexp(1)
  ), data = d2
)
precis(m6.9, depth = 2)


m6.10 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a + bA*A,
    a ~ dnorm(0, 1),
    bA ~ dnorm(0, 2),
    sigma ~ dexp(1)
  ), data = d2
)

precis(m6.10, depth = 2)

N <- 200
b_GP <- 1
b_GC <- 0.25
b_PC <- 1
b_U <- 2

set.seed(1)
U <- 2*rbern(N, 0.5) -1 
G <- (2*rbern(N, 0.5) -1)*0.5
G <- rnorm(N)
P <- rnorm(N, b_GP*G + b_U*U)
C <- rnorm(N, b_PC*P + b_GC*G + b_U*U)
d <- data.frame(C=C, P=P, G=G, U=U)


m6.11b <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P,
    a ~ dnorm(0, 1),
    b_PC ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=d 
)
precis(m6.11)


m6.11c <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_GC*G,
    a ~ dnorm(0, 1),
    b_GC ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=d 
)
precis(m6.11)

m6.11 <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G,
    a ~ dnorm(0, 1),
    c(b_PC, b_GC) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=d 
)
precis(m6.11)

m6.12b <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_U*U,
    a ~ dnorm(0, 1),
    b_U ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=d 
)

m6.12c <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G + b_U*U,
    a ~ dnorm(0, 1),
    c(b_GC, b_U) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=d 
)

m6.12 <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G + b_U*U,
    a ~ dnorm(0, 1),
    c(b_PC, b_GC, b_U) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=d 
)
precis(m6.12)


plot(coeftab(m6.11, m6.12), par=c('b_PC', 'b_GC', 'b_U'))


plot(coeftab(m6.11, m6.11c, m6.12, m6.12c), par=c('b_PC', 'b_GC', 'b_U'))


plot(coeftab(m6.12, m6.12b, m6.12c), par=c('b_PC', 'b_GC', 'b_U'))

library(rethinking)
library(dagitty)
dag_6.1 <- dagitty('dag {
      U [unobserved]
      X -> Y
      X <- U <- A -> C -> Y
      U -> B <- C
  }')
adjustmentSets(dag_6.1, exposure='X', outcome='Y')

dag_6.1b <- dagitty('dag {
      U [unobserved]
      V [unobserved]
      X -> Y
      X <- U <- A -> C -> Y
      U -> B <- C
      C <- V -> Y
  }')
adjustmentSets(dag_6.1b, exposure='X', outcome='Y')

dag_6M2<- dagitty('dag {
      X -> Z -> Y
  }')
adjustmentSets(dag_6M2, exposure='X', outcome='Y')

# X -> Z -> Y
n <- 100
X <- rnorm(n)
Z <- rnorm(n, -X)
Y <- rnorm(n, -Z)
#Z <- rnorm(n)
#Z <- rnorm(Y)
d_sim6M2 <- data.frame(X=X, Y=Y, Z=Z)
precis(d_sim6M2, hist=FALSE)

m6M.2a <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + bt*X + bf*Z,
    a ~ dnorm(0, 0.5),
    bt ~ dnorm(0, 0.25),
    bf ~ dnorm(0, 0.25),
    sigma ~ dexp(1)
  ), data = d_sim6M2
)

m6M.2b <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + bt*X,
    a ~ dnorm(0, 0.5),
    bt ~ dnorm(0, 0.25),
    sigma ~ dexp(1)
  ), data = d_sim6M2
)

m6M.2c <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + bf*Z,
    a ~ dnorm(0, 0.5),
    bf ~ dnorm(0, 0.25),
    sigma ~ dexp(1)
  ), data = d_sim6M2
)

precis(m6M.2a)
precis(m6M.2b)
precis(m6M.2c)

compare(m6M.2a, m6M.2b, m6M.2c)
plot(compare(m6M.2a, m6M.2b, m6M.2c))


n <- 100
A <- rnorm(n)
U <- rnorm(n, A)
C <- rnorm(n, A)
B <- rnorm(n, U-C)
X <- rnorm(n, U)
Y <- rnorm(n, X-C)
d <- data.frame(A=A, U=U, C=C, B=B, X=X, Y=Y)
##
V <- rnorm(n)
C <- rnorm(n, A-V)
Y <- rnorm(n, X-C-V)
##
d <- data.frame(A=A, U=U, C=C, B=B, X=X, Y=Y)

m6.4.2a <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + bX*X + bC*C + bA*A + bB*B,
    a ~ dnorm(0, 1),
    c(bX, bC, bA, bB) ~ dnorm(0,0.25),
    sigma ~ dexp(1)
  ), data= d
)

m6.4.2b <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + bX*X,
    a ~ dnorm(0, 1),
    bX ~ dnorm(0,0.25),
    sigma ~ dexp(1)
  ), data= d
)

m6.4.2c <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + bB*B,
    a ~ dnorm(0, 1),
    bB ~ dnorm(0,0.25),
    sigma ~ dexp(1)
  ), data= d
)

m6.4.2d <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + bA*A,
    a ~ dnorm(0, 1),
    bA ~ dnorm(0,0.25),
    sigma ~ dexp(1)
  ), data= d
)

m6.4.2e <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + bX*X + bC*C,
    a ~ dnorm(0, 1),
    c(bX, bC) ~ dnorm(0,0.25),
    sigma ~ dexp(1)
  ), data= d
)


m6.4.2f <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + bX*X + bB*B,
    a ~ dnorm(0, 1),
    c(bX, bB) ~ dnorm(0,0.25),
    sigma ~ dexp(1)
  ), data= d
)

m6.4.2g <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + bX*X + bA*A,
    a ~ dnorm(0, 1),
    c(bX, bA) ~ dnorm(0,0.25),
    sigma ~ dexp(1)
  ), data= d
)

plot(coeftab(m6.4.2a, m6.4.2e, m6.4.2f, m6.4.2g), pars=c('bX', 'bB', 'bA', 'bC'))

plot(coeftab(m6.4.2a, m6.4.2b, m6.4.2e, m6.4.2f, m6.4.2g), pars=c('bX', 'bB', 'bA', 'bC'))

plot(coeftab(m6.4.2a, m6.4.2d, m6.4.2g), pars=c('bX', 'bB', 'bA', 'bC'))

plot(coeftab(m6.4.2a, m6.4.2c, m6.4.2f), pars=c('bX', 'bB', 'bA', 'bC'))

compare(m6.4.2a, m6.4.2c, m6.4.2b, m6.4.2d,m6.4.2e, m6.4.2f, m6.4.2g, func=WAIC)


compare(m6.4.2a, m6.4.2c, m6.4.2b, m6.4.2d,m6.4.2e, m6.4.2f, m6.4.2g, func=PSIS)

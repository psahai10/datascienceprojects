library(rethinking)
library(rstan)

p <- list()
p[[1]] <- c(1/4,1/4,1/4,1/4)
p[[2]] <- c(2/6,1/6,1/6,2/6)
p[[3]] <- c(1/6,2/6,2/6,1/6)
p[[4]] <- c(1/8,4/8,2/8,1/8)


sapply(p, function(p) sum(p*c(0,1,1,2)) )

sapply(p, function(p) -sum(p*log(p) ) )

p <- 0.7
(A <- c( (1-p)^2, p*(1-p), (1-p)*p, p^2 ) )

-sum(A*log(A))

sim.p <- function(G=1.4) {
  x123 <- runif(3)
  x4 <- ( (G)*sum(x123)-x123[2]-x123[3] ) / (2-G)
  z <- sum( c(x123, x4) )
  p <- c( x123, x4 )/z
  list(H = -sum(p*log(p) ), p=p)
}

H <- replicate(1e5, sim.p(1.4))

dens( as.numeric(H[1,]), adj=0.1)

entropies <- as.numeric(H[1,])
dist <- H[2,]

max(entropies)

dist[ which.max(entropies)]


library(rethinking)
data("chimpanzees")
d <- chimpanzees

d$treatment <- 1 + d$prosoc_left + 2 * d$condition

xtabs( ~ treatment + prosoc_left + condition, d)

m11.1 <- quap(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a,
    a ~ dnorm(0, 1.5)
  ), data=d
)
set.seed(1999)
prior <- extract.prior(m11.1, n=1e4)

p <- inv_logit(prior$a)

dens(p, adj=0.1)

m11.2 <- quap(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a + b[treatment] ,
    a ~ dnorm(0, 1.5),
    b[treatment] ~ dnorm(0, 10)
  ), data=d
)
set.seed(1999)
prior <- extract.prior(m11.2, n=1e4)
p <- sapply(1:4, function(k) inv_logit(prior$a + prior$b[,k]))

dens( abs( p[,1] - p[,2]), adj=0.1)
dens( abs( p[,3] - p[,4]), adj=0.1)

m11.3 <- quap(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a + b[treatment] ,
    a ~ dnorm(0, 1.5),
    b[treatment] ~ dnorm(0, 0.5)
  ), data=d
)
set.seed(1999)
prior <- extract.prior(m11.3, n=1e4)
p <- sapply(1:4, function(k) inv_logit(prior$a + prior$b[,k]))

dens( abs( p[,1] - p[,2]), adj=0.1)
mean( abs( p[,1] - p[,2]), adj=0.1)

dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = as.integer(d$treatment)
)

m11.4 <- ulam(
  alist(
    pulled_left ~ dbinom(1,p),
    logit(p) <- a[actor] + b[treatment],
    a[actor] ~ dnorm(0, 1.5),
    b[treatment] ~ dnorm(0, 0.5)
  ), data=dat_list, chains = 4, log_lik = TRUE
  
)

precis(m11.4, depth = 2)

m11.4_stan_code <- stancode(m11.4)

post <- extract.samples(m11.4)
p_left <- inv_logit(post$a)
plot(precis(as.data.frame(p_left)), xlim=c(0,1))


labs <- c("R/N", "L/N", "R/P", "L/P")
plot(precis(m11.4, depth=2, pars='b'), labels=labs)

diffs <- list(
  db12 = post$b[,1] - post$b[,3],
  db24 = post$b[,2] - post$b[,4])
plot(precis(diffs))

pl <- by(d$pulled_left, list(d$actor, d$treatment), mean)
pl[5,]


plot(NULL, xlim=c(1,28), ylim=c(0,1), xlab='',
     ylab='proportion left lever', xaxt='n', yaxt='n')
axis( 2, at=c(0, 0.5, 1), labels=c(0,0.5,1))
abline(h=0.5, lty=2)
for (j in 1:7) abline( v=(j-1)*4+4.5, lwd=0.5)
for (j in 1:7) text( (j-1)*4+2.5, 1.07, concat('actor ', j), xpd=TRUE)
for (j in (1:7)[-2]) {
  lines( (j-1)*4+c(1,3), pl[j,c(1,3)] , lwd=2, col=rangi2)
  lines( (j-1)*4+c(2,4), pl[j,c(2,4)] , lwd=2, col=rangi2)
}
yoff <- 0.01
points(1:28, t(pl), pch=16, col='white', cex=1.7)
points(1:28, t(pl), pch=c(1,1,16,16), col=rangi2, lwd=2)


###

dat <- list(actor=rep(1:7, each=4), treatment=rep(1:4, times=7))
p_post <- link(m11.4, data=dat)
p_mu <- apply(p_post, 2, mean)
p_ci <- apply(p_post, 2, PI)
dat$mu <- p_mu
dat$ci <- p_ci 



plot(NULL, xlim=c(1,28), ylim=c(0,1), xlab='',
     ylab='proportion left lever', xaxt='n', yaxt='n')
axis( 2, at=c(0, 0.5, 1), labels=c(0,0.5,1))
abline(h=0.5, lty=2)
for (j in 1:7) abline( v=(j-1)*4+4.5, lwd=0.5)
for (j in 1:7) text( (j-1)*4+2.5, 1.07, concat('actor ', j), xpd=TRUE)
x=0
for (j in (1:7)[-2]) {
  lines( (j-1)*4+c(1,3), c(dat[[3]][j+x],dat[[3]][j+2+x]) , lwd=2, col=rangi2)
  lines( (j-1)*4+c(2,4), c(dat[[3]][j+1+x],dat[[3]][j+3+x]) , lwd=2, col=rangi2)
  x+3
}

yoff <- 0.01
points(dat$mu, pch=16, col='white', cex=1.7)
points(dat$mu, pch=c(1,1,16,16), col=rangi2, lwd=2)

data("chimpanzees")
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2 * d$condition


d$side <- d$prosoc_left+1
d$cond <- d$condition+1

d_aggregated <- aggregate(
  d$pulled_left,
  list(treatment=d$treatment, actor=d$actor,
       side=d$side, cond=d$cond), sum)


colnames(d_aggregated)[5] <- 'left_pulls'

d_aggregated

dat <- with( d_aggregated , list(
  left_pulls = left_pulls,
  treatment = treatment,
  actor = actor,
  side = side,
  cond = cond ) )

m11.6 <- ulam(
  alist(
    left_pulls ~ dbinom( 18 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm( 0 , 1.5 ) ,
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat , chains=4 , log_lik=TRUE )
compare(m11.6, m11.4, func=PSIS)

-2*dbinom(6,9,0.2, log=TRUE)

-2*sum(dbern(c(1,1,1,1,1,1,0,0,0), 0.2, log=TRUE))


library(rethinking)
data(UCBadmit)
d <- UCBadmit
head(UCBadmit)

d$gid <- ifelse(d$applicant.gender=='male', 1, 2)

dat_list <- list(
  admit=d$admit,
  applications = d$applications,
  gid = d$gid
)

m11.7 <- ulam(
  alist(
    admit ~ dbinom(applications, p),
    logit(p) <- a[gid],
    a[gid] ~ dnorm(0, 1.5)
  ) , data=dat_list, chains=4
)
precis(m11.7, depth = 2)


post <- extract.samples(m11.7)
diff_a <- post$a[,1] - post$a[,2]
diff_b <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
precis(list(diff_a=diff_a,diff_p=diff_b), hist=FALSE)

postcheck(m11.7)

for ( i in 1:6 ) {
  x <- 1 + 2*(i-1)
  y1 <- d$admit[x]/d$applications[x]
  y2 <- d$admit[x+1]/d$applications[x+1]
  lines( c(x,x+1) , c(y1,y2) , col=rangi2 , lwd=2 )
  text( x+0.5 , (y1+y2)/2 + 0.05 , d$dept[x] , cex=0.8 , col=rangi2 )
}

dat_list_B <- list(
  admit=d$admit,
  applications = d$applications,
  gid = d$gid,
  dept_id = rep(1:6, each=2)
)

m11.8 <- ulam(
  alist(
    admit ~ dbinom(applications, p),
    logit(p) <- a[gid] + delta[dept_id],
    a[gid] ~ dnorm(0, 1.5),
    delta[dept_id] ~ dnorm(0, 1.5)
  ) , data=dat_list_B, chains=4, iter=4000
)
precis(m11.8, depth = 2)

post <- extract.samples(m11.8)
diff_a <- post$a[,1] - post$a[,2]
diff_b <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
precis(list(diff_a=diff_a,diff_p=diff_b), hist=FALSE)


postcheck(m11.8)

for ( i in 1:6 ) {
  x <- 1 + 2*(i-1)
  y1 <- d$admit[x]/d$applications[x]
  y2 <- d$admit[x+1]/d$applications[x+1]
  lines( c(x,x+1) , c(y1,y2) , col=rangi2 , lwd=2 )
  text( x+0.5 , (y1+y2)/2 + 0.05 , d$dept[x] , cex=0.8 , col=rangi2 )
}

pg <- with( dat_list_B , sapply( 1:6 , function(k)
  applications[dept_id==k]/sum(applications[dept_id==k]) ) )
rownames(pg) <- c("male","female")
colnames(pg) <- unique(d$dept)
round( pg , 2 )

m11.8b <- ulam(
  alist(
    admit ~ dbinom(applications, p),
    logit(p) <- a + delta[dept_id],
    a ~ dnorm(0, 1.5),
    delta[dept_id] ~ dnorm(0, 1.5)
  ) , data=dat_list_B, chains=4, iter=4000
)
precis(m11.8b, depth = 2)

post <- extract.samples(m11.8b)
diff_a <- post$a[,1] - post$a[,2]
diff_b <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
precis(list(diff_a=diff_a,diff_p=diff_b), hist=FALSE)


postcheck(m11.8b)

for ( i in 1:6 ) {
  x <- 1 + 2*(i-1)
  y1 <- d$admit[x]/d$applications[x]
  y2 <- d$admit[x+1]/d$applications[x+1]
  lines( c(x,x+1) , c(y1,y2) , col=rangi2 , lwd=2 )
  text( x+0.5 , (y1+y2)/2 + 0.05 , d$dept[x] , cex=0.8 , col=rangi2 )
}

y <- rbinom(10000, 1000, 1/1000)

c(mean(y), var(y))

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
  P=d$P,
  cid=d$contact_id
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
    log(lambda) <- a[cid] + b[cid]*P,
    a[cid] ~ dnorm(3, 0.5),
    b[cid] ~ dnorm(0, 0.2)
  ), data=dat, chains=4, log_lik = TRUE
)

compare(m11.9, m11.10, func=PSIS)

k <- PSIS( m11.10, pointwise = TRUE)$k
plot( dat$P, dat$T, xlab='log Population (std)', ylab='total tools',
      col=rangi2, pch=ifelse(dat$cid==1,1,16), lwd=2,
      ylim=c(0,75), cex=1+normalize(k))

ns <- 100
P_seq <- seq( from=-1.4, to=3, length.out = ns)

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


plot( d$population, d$total_tools, xlab='Population', ylab='total tools',
      col=rangi2, pch=ifelse(dat$cid==1,1,16), lwd=2,
      ylim=c(0,75), cex=1+normalize(k))

ns <- 100
P_seq <- seq( from=-5, to=3, length.out = ns)
pop_seq <- exp(P_seq*1.53+9)
P_seq2 <- log(pop_seq)/1.53-9
lambda <- link(m11.10, data=data.frame(P=P_seq, cid=1 ) )
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines( pop_seq, lmu, lty=2, lwd=1.5)
shade(lci, pop_seq, xpd=TRUE)

lambda <- link(m11.10, data=data.frame(P=P_seq, cid=2 ) )
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines( pop_seq, lmu, lty=1, lwd=1.5)
shade(lci, pop_seq, xpd=TRUE)


N <- 100
a <- rnorm(N, 1, 1)
b <- rexp(N, 1)
g <- rexp(N, 1)
plot(NULL, xlim=c(-2,2), ylim=c(0,100))

x_seq <- seq(from=log(100), to=log(2e5), length.out=100)
lambda <- sapply(x_seq, function(x) exp(a)*x^b/g)

plot(NULL, xlim=range(exp(x_seq)), ylim=c(0,500), xlab='population', ylab='total tools')
for (i in 1:N) lines(exp(x_seq), lambda[i,], col=grau(), lwd=1.25)


dat2 <- list(T=d$total_tools, P=d$population, cid=d$contact_id)

m11.11 <- ulam(
  alist(
    T ~ dpois( lambda ),
    lambda <- exp(a[cid])*P^b[cid]/g,
    a[cid] ~ dnorm(1, 1),
    b[cid] ~ dexp(1),
    g ~ dexp(1)
  ), data=dat2, chains=4, log_lik = TRUE
)

k <- PSIS( m11.11, pointwise = TRUE)$k

plot( d$population, d$total_tools, xlab='Population', ylab='total tools',
      col=rangi2, pch=ifelse(dat$cid==1,1,16), lwd=2,
      ylim=c(0,75), cex=1+normalize(k))

ns <- 100
P_seq <- seq( from=0, to=6, length.out = ns)
pop_seq <- exp(P_seq*1.53+9)

lambda <- link(m11.11, data=data.frame(P=P_seq, cid=1 ) )
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines( pop_seq, lmu, lty=2, lwd=1.5)
shade(lci, pop_seq, xpd=TRUE)

lambda <- link(m11.11, data=data.frame(P=P_seq, cid=2 ) )
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines( pop_seq, lmu, lty=1, lwd=1.5)
shade(lci, pop_seq, xpd=TRUE)

library(rethinking)
num_days <- 30
y <- rpois( num_days, 1.5)
mean(y)

num_weeks <- 4
y_new <- rpois( num_weeks, 0.5*7)
mean(y_new)

y_all <- c(y, y_new)
exposure <- c(rep(1, 30), rep(7, 4))
monastery <- c(rep(0,30), rep(1,4))
d <- data.frame(y=y_all, days=exposure, monastery=monastery)

d$log_days <- log(d$days)

m11.12 <- quap(
  alist(
    y ~ dpois( lambda ),
    log(lambda) <- log_days + a + b*monastery,
    a ~ dnorm(0, 1),
    b ~ dnorm(0, 1)
  ), data=d
)

m11.12b <- ulam(
  alist(
    y ~ dpois( lambda ),
    log(lambda) <- log_days + a + b*monastery,
    a ~ dnorm(0, 1),
    b ~ dnorm(0, 1)
  ), data=d
)

precis(m11.12)
precis(m11.12b)


post <- extract.samples(m11.12)
lambda_old <- exp(post$a)
lambda_new <- exp(post$a + post$b)
precis(data.frame(lambda_old, lambda_new), hist=FALSE)

# simulate career choices among 500 individuals
N <- 500             # number of individuals
income <- c(1,2,5)   # expected income of each career
score <- 0.5*income  # scores for each career, based on income
# next line converts scores to probabilities
p <- softmax( score[1], score[2], score[3] )

# now simulate choice
# outcome career holds event type values, not counts
career <- rep( NA, N )  # empty vector of choices for each individual
# sample chosen career for each individual
set.seed(34302)
for ( i in 1:N ) career[i] <- sample( 1:3 , size=1 , prob=p )


m11.13 <- quap(
  alist(
    career ~ dcategorical(softmax(s1,s2,s3)),
    s1 <- a1 + b*1, # 1 = career income [1]   # linear model for event type 2
    s2 <- a2 + b*2, # 2 = career income [2]   # linear model for event type 3
    s3 <- 0,
    c(a1, a2) ~ dnorm(0,1),
    b ~ dnorm(0, 0.5)
  ) ,
  data=list(career=career) )
precis(m11.13, 2)

## R code 11.56
code_m11.13 <- "
data{
    int N; // number of individuals
    int K; // number of possible careers
    int career[N]; // outcome
    vector[K] career_income;
}
parameters{
    vector[K-1] a; // intercepts
    real<lower=0> b; // association of income with choice
}
model{
    vector[K] p;
    vector[K] s;
    a ~ normal( 0 , 1 );
    b ~ normal( 0 , 0.5 );
    s[1] = a[1] + b*career_income[1];
    s[2] = a[2] + b*career_income[2];
    s[3] = 0; // pivot
    p = softmax( s );
    career ~ categorical( p );
}
"

## R code 11.57
dat_list <- list( N=N , K=3 , career=career , career_income=income )
m11.13 <- stan( model_code=code_m11.13 , data=dat_list , chains=4)
precis( m11.13 , 2 )

post <- extract.samples(m11.13)
s1 <- with( post, a[,1] + b*income[1] )
s2_orig <- with( post, a[,2] + b*income[2] )
s2_new <- with( post, a[2] + b*income[2]*2)

p_orig <- sapply( 1:length(post$b), function(i)
  softmax( c(s1[i], s2_orig[i], 0) ) )

p_new <- sapply( 1:length(post$b), function(i)
  softmax( c(s1[i], s2_new[i], 0) ) )

p_diff <- p_new[2,] - p_orig[2,]
precis(p_diff, hist=FALSE)


N <- 500
# simulate family incomes for each individual
family_income <- runif(N)
# assign a unique coefficient for each type of event
b <- c( -2, 0 ,2 )
career <- rep( NA , N )  # empty vector of choices for each individual
for ( i in 1:N ) {
  score <- 0.5*(1:3) + b*family_income[i]
  p <- softmax( score[1] , score[2] , score[3] )
  career[i] <- sample( 1:3 , size=1 , prob=p )
}

m11.14 <- quap(
  alist(
    career ~ dcategorical(softmax(s1,s2,s3)),
    s1 <- a1 + b1*family_income,
    s2 <- a2 + b2*family_income,
    s3 <- 0,
    c(a1,a2) ~ dnorm(0,1.5),
    c(b1,b2) ~ dnorm(0,1)),
  data=list(career=career,family_income=family_income))

precis(m11.14)


code_m11.14 <- "
data{
    int N; // number of observations
    int K; // number of outcome values
    int career[N]; // outcome
    real family_income[N];
}
parameters{
    vector[K-1] a; // intercepts
    vector[K-1] b; // coefficients on family income
}
model{
    vector[K] p;
    vector[K] s;
    a ~ normal(0,1.5);
    b ~ normal(0,1);
    for ( i in 1:N ) {
        for ( j in 1:(K-1) ) s[j] = a[j] + b[j]*family_income[i];
        s[K] = 0; // the pivot
        p = softmax( s );
        career[i] ~ categorical( p );
    }
}
"

dat_list <- list( N=N , K=3 , career=career , family_income=family_income )
m11.14 <- stan( model_code=code_m11.14 , data=dat_list , chains=4 )
precis( m11.14 , 2 )

## R code 11.60
library(rethinking)
data(UCBadmit)
d <- UCBadmit

## R code 11.61
# binomial model of overall admission probability
m_binom <- quap(
  alist(
    admit ~ dbinom(applications,p),
    logit(p) <- a,
    a ~ dnorm( 0 , 1.5 )
  ), data=d )

library(rethinking)
data("UCBadmit")
d <- UCBadmit

m_binom <- quap(
  alist(
    admit ~ dbinom(applications, p),
    logit(p) <- a,
    a ~ dnorm(0, 1.5)
  ), data=d)

dat <- list( admit=d$admit, rej=d$reject)
m_pois <- ulam(
  alist(
    admit ~ dpois(lambda1),
    rej ~ dpois(lambda2),
    log(lambda1) <- a1,
    log(lambda2) <- a2,
    a(a1, a2) ~ dnorm(0, 1.5)
  ), data=dat, chains=4, cores=4
)


inv_logit(coef(m_binom))

k <- coef(m_pois)
a1 <- k['a1']
a2 <- k['a2']

exp(a1)/(exp(a1)+exp(a2))

logit(0.35)

inv_logit(3.2)

library(rethinking)
data("UCBadmit")
d <- UCBadmit

m_binom <- map2stan(
  alist(
      admit ~ dbinom(applications, p),
      logit(p) <- a,
      a ~ dnorm(0, 100)
  ), data=d
)

d$rej <- d$reject

m_pois <- map2stan(
  alist(
    admit ~ dpois(lambda1),
    rej ~ dpois(lambda2),
    log(lambda1) <- a1,
    log(lambda2) <- a2,
    c(a1, a2) ~ dnorm(0,100)
  ), data=d, chains=4, cores=4
)

N <- 2
x <- replicate(1e5, min(runif(2,1,100)))
y <- replicate(1e5, min(runif(5,1,100)))
d1 <- data.frame(dens=x, val='a')
d2 <- data.frame(dens=y, val='b')
d <- rbind(d1, d2)

library(ggplot2)

ggplot(d, aes(x=dens, fill=val)) + geom_density(alpha=0.5)

x <- replicate(1e5, sort(runif(10,1,100))[2])
y <- replicate(1e5, sort(runif(10,1,100))[5])
d1 <- data.frame(dens=x, val='a')
d2 <- data.frame(dens=y, val='b')
d <- rbind(d1, d2)

ggplot(d, aes(x=dens, fill=val)) + geom_density(alpha=0.5)

library(rethinking)
data("AustinCats")
d <- AustinCats

d$adopt <- ifelse( d$out_event=="Adoption", 1L, 0L)

dat <- list(
  days_to_event = as.numeric(d$days_to_event),
  color_id = ifelse(d$color=='Black', 1L, 2L),
  adopted = d$adopt
)

m11.14 <- ulam(
  alist(
    days_to_event|adopted==1 ~ exponential(lambda),
    days_to_event|adopted==0 ~ custom(exponential_lccdf(!Y|lambda)),
    lambda <- 1.0/mu,
    log(mu) <- a[color_id],
    a[color_id] ~ normal(0,1)
  ), data= dat, chains=4, cores=4
)

precis(m11.14, depth=2, hist=FALSE)

post <- extract.samples(m11.14)
post$D <- exp(post$a)
precis(post, 2, hist=FALSE)

library(rethinking)
ns <- 200
P_seq <- seq(from=0, 200, length.out = 200)
# Predictions for black cats
library(rethinking)
N <- 500
income <- c(1,2,5)
score <- 0.5*income

pp <- softmax(score[1], score[2], score[3])
career <- rep(NA,N)
for (i in 1:N) career[i] <- sample(1:3, size=1, prob=p)

d <- data.frame(career=career, income=0, scores=0, p=0)
library(dplyr)

d <- d %>%
  mutate(income = ifelse(career==3, 5, ifelse(career==2, 2, 1)))

d <- d %>%
  mutate(scores = ifelse(career==3, score[3], ifelse(career==2, score[2], score[1])))

d <- d %>%
  mutate(p = ifelse(career==3, pp[3], ifelse(career==2, pp[2], pp[1])))

dat_list <- list( N=N , K=3 , career=career , income=income )

m11.13a <- ulam(
  alist(
    career ~ categorical( p ),
    p <- softmax(score[career]),
    score[career] ~ a[career] + b*income[career],
    a[career] ~ dnorm(0,1),
    b ~ dnorm(0, 0.5)
  ), data=d, chains=4, cores=4
)

m11.13a <- ulam(
  alist(
    career ~ dmultinom(x=3, N=500, p),
    p <- softmax(score[career]),
    score[career] ~ a[career] + b*income[career],
    a[career] ~ dnorm(0,1),
    b ~ dnorm(0, 0.5)
  ), data=d, chains=4, cores=4
)

m11.13b <- ulam(
  alist(
    income ~ dpois(lambda[career]),
    log(lambda[career]) <- a[career],
    a[career] ~ dnorm(0 , 1.5),
  ), data=d, chains=4, cores=4
)


m11.4_stan_code <- stancode(m11.4)



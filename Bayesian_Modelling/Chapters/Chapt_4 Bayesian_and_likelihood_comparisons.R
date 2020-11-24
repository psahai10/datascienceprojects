library(BAS)
data(bodyfat)
summary(bodyfat)

# lm Fit linear Models
# Frequentist OLS linear regression
bodyfat.lm = lm(Bodyfat ~ Abdomen, data = d)
summary(bodyfat.lm)

# Extract coefficients
beta = coef(bodyfat.lm)

# Visualize regression line on the scatter plot
library(ggplot2)
ggplot(data = bodyfat, aes(x = Abdomen, y = Bodyfat)) +
  geom_point(size=3, color = "blue", alpha=.3) +
  geom_abline(intercept = beta[1], slope = beta[2], size = 1) +
  xlab("abdomen circumference (cm)") 

library(rethinking)
data(Howell1)
d <- Howell1
d <- d[d$age>=18,]
heights.lm = lm(height~weight, data=d)
summary(heights.lm)
beta = coef(heights.lm)

ggplot(data = d, aes(x = weight, y = height)) +
  geom_point(size=3, color = "blue", alpha=.3) +
  geom_abline(intercept = beta[1], slope = beta[2], size = 1) +
  xlab("weight (kg)") 

resid = residuals(heights.lm)
n = length(resid)

# Calculate MSE
MSE = 1/ (n - 2) * sum((resid ^ 2))
MSE

result = data.frame(fitted_values = fitted.values(heights.lm),
                    residuals = residuals(heights.lm))

# Load library and plot residuals versus fitted values
ggplot(data = result, aes(x = fitted_values, y = residuals)) +
  geom_point(pch = 1, size = 2) + 
  geom_abline(intercept = 0, slope = 0) + 
  xlab(expression(paste("fitted value ", widehat(heights)))) + 
  ylab("residuals")

# Find the observation with the largest fitted value
which.max(as.vector(fitted.values(heights.lm)))

# Shows this observation has the largest Abdomen
which.max(d$weight)

plot(heights.lm, which = 2)

summary(bodyfat)
d <- bodyfat
psych::describe(d)
hist(d)
names(d)
dF <- data.frame(Bodyfat = d$Bodyfat, Abdomen = d$Abdomen)
psych::describe(dF)
dens(dF$Bodyfat)
curve(dnorm(x, 19, 6), from=-10, to= 50)
curve(dunif(x, 0, 20), from=-10, to=30)

sample_mu <- rnorm(1e4, 19, 6)
sample_sigma <- runif(1e4, 0, 20)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)

# Step1: generate a sequence of mu and sigma
mu.list <- seq(15, 24, length.out = 100)
sigma.list <- seq(7, 10, length.out = 100)

# Step2: Combine this sequence in a grid that will equal length = 100X100
post <- expand.grid(mu=mu.list, sigma=sigma.list)

# Create a function that takes the bodyfat from the data, runs a dnorm
# simulation for all bodyfat with respective mu and sigma and sums it up
# and return it as sum total
func <- function(x) { sum(
  dnorm(dF$Bodyfat, post$mu[x], post$sigma[x], log=TRUE))
}

post$LL <- sapply(1:nrow(post), func)

post$prod <- post$LL + dnorm( post$mu, 19, 6, TRUE) +
  dunif(post$sigma, 0, 20, TRUE)

post$prob <- exp(post$prod - max(post$prod))

contour_xyz(post$mu, post$sigma, post$prob)
image_xyz(post$mu, post$sigma, post$prob)


sample.rows <- sample(1:nrow(post), size=1e4, replace=TRUE, prob=post$prob)
sample.mu <- post$mu[ sample.rows]
sample.sigma <- post$sigma[ sample.rows]
sample.rows


plot(sample.mu, sample.sigma, cex=0.5, pch=16, col=col.alpha(rangi2,0.4))

dens(sample.mu)
dens(sample.sigma)

PI(sample.mu)
PI(sample.sigma)
HPDI(sample.mu)
HPDI(sample.sigma)




mBF <- quap(
  alist(
    Bodyfat ~ dnorm(mu, sigma),
    mu ~ dnorm(19, 6),
    sigma ~ dunif(0, 20)
  ) , data=dF
)

precis(mBF)
precis(mBF, prob=0.99)

vcov(mBF)
diag( vcov(mBF))
cov2cor(vcov(mBF))

post <- extract.samples(mBF, n=1e4)
precis(post, hist=FALSE)
plot(post)

library(MASS)
post <- mvrnorm(1e4, mu=coef(mBF), Sigma=vcov(mBF))
plot(dF$Bodyfat~dF$Abdomen)


set.seed(2971)
N <- 100
a <- rnorm(N, 19, 6)
b <- rlnorm(N, 0, 1)

plot( NULL, xlim=range(dF$Abdomen), ylim=c(-20, 80),
      xlab='Abdomen', ylab='Bodyfat')
abline(h=0, lty=2)
abline(h=40, lty=1, lwd=0.5)
mtext('b ~ dnorm(0, 10)')
xbar <- mean(dF$Abdomen)
for (i in 1:N) curve(a[i] + b[i]*(x-xbar),
                     from=min(dF$Abdomen), to=max(dF$Abdomen), add=TRUE,
                     col=col.alpha('black', 0.2) )

b <-rlnorm(1e4, 0, 1)
dens(b, xlim=c(0,10), adj=0.1)


xbar <- mean(dF$Abdomen)
mBF.2 <- quap(
  alist(
    Bodyfat ~ dnorm(mu, sigma),
    mu <- a + b*(Abdomen - xbar),
    a ~ dnorm(19, 3),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 20)
  ) , data=dF)

precis(mBF.2)
round(vcov(mBF.2), 3)
pairs(mBF.2)

xbar <- mean(dF$Abdomen)

plot(Bodyfat~Abdomen, data=dF, col=rangi2)
post<-extract.samples(mBF.2)
a_map <-mean(post$a)
b_map <-mean(post$b)
curve(a_map, + b_map*(x-xbar), add=TRUE)


N <- 10
dN <- dF[1:N, ]
mN <- quap(
  alist(
    Bodyfat ~ dnorm(mu, sigma),
    mu <- a + b*(Abdomen-mean(Abdomen)),
    a ~ dnorm(19, 6),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0, 20)
  ), data=dN
)

post <- extract.samples(mN, n=100)
plot(dN$Abdomen, dN$Bodyfat,
     xlim=range(dN$Abdomen), ylim=range(dN$Bodyfat), 
     col=rangi2, xlab='Abdomen', ylab='Bodyfat')
mtext(concat("N= ", N))

for (i in 1:20) curve(post$a[i] + post$b[i]*(x-mean(dN$Abdomen)), 
                       col=col.alpha('black', 0.3), add=TRUE)

post <- extract.samples(mBF.2)
mu_at_88 <- post$a + post$b *(88-xbar)
dens(mu_at_88, col=rangi2, lwd=2, xlab='mu|Abdomen=88')

mu<-link(mBF.2)
str(mu)
######################################
Abdomen.seq <- seq(from=60, to=140, by=1)
mu <- link(mBF.2, data=data.frame(Abdomen=Abdomen.seq))
str(mu)
plot(Bodyfat~Abdomen, dF, type='n')
for (i in 1:100) points(Abdomen.seq, mu[i,], pch=16, col=col.alpha(rangi2,0.1))

######################################

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)
plot(Bodyfat~Abdomen, data=dF, col=col.alpha(rangi2, 0.5))
lines(Abdomen.seq, mu.mean)
shade(mu.PI, Abdomen.seq)

##########################################

#sim.Bodyfat <- sim( mBF.2, data=list(Abdomen=Abdomen.seq) )
str(sim.Bodyfat)

sim.Bodyfat <- sim(mBF.2, data=list(Abdomen=Abdomen.seq), n=1e4)
#str(sim.Bodyfat)

Bodyfat.PI <- apply(sim.Bodyfat, 2, PI, prob=0.89)

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)

plot(Bodyfat~Abdomen, data=dF, col=col.alpha(rangi2, 0.5))
lines(Abdomen.seq, mu.mean)
shade(mu.PI, Abdomen.seq)
shade(Bodyfat.PI, Abdomen.seq)


######################################
dF$Abdomen_s <- (dF$Abdomen-mean(dF$Abdomen))/sd(dF$Abdomen)
dF$Abdomen_s2 <- dF$Abdomen_s^2

mBF4.5 <- quap(
  alist(
    Bodyfat ~ dnorm(mu, sigma),
    mu <- a + b1*Abdomen_s + b2*Abdomen_s2,
    a ~ dnorm(19, 3),
    b1 ~ dlnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ) , data=dF)
precis(mBF4.5)

Abdomen.seq <- seq(from=-2.2, to=4, length.out=40)
pred_dat <- list(Abdomen_s=Abdomen.seq, Abdomen_s2=Abdomen.seq^2)
mu <- link(mBF4.5, data=pred_dat)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)
sim.Bodyfat <- sim(mBF4.5, data=pred_dat)
Bodyfat.PI <- apply(sim.Bodyfat, 2, PI, prob=0.89)

plot(Bodyfat~Abdomen_s, data=dF, col=col.alpha(rangi2, 0.5))
lines(Abdomen.seq, mu.mean)
shade(mu.PI, Abdomen.seq)
shade(Bodyfat.PI, Abdomen.seq)

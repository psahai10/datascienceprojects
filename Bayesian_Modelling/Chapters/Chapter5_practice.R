library(rethinking)
getwd()
setwd("C:/Users/psahai/Documents/R/R_Bayesian_Statistical_Rethinking/practice")
getwd()
library("readxl")
data <- read_excel("Concrete_Data.xls")
head(data)
names(data)

d <- data.frame(Strenght = data$'Concrete compressive strength(MPa, megapascals)', 
                Cement = data$'Cement (component 1)(kg in a m^3 mixture)', 
                Blast = data$'Blast Furnace Slag (component 2)(kg in a m^3 mixture)',
                Fly = data$'Fly Ash (component 3)(kg in a m^3 mixture)',
                Water = data$'Water  (component 4)(kg in a m^3 mixture)',
                Plasticizer = data$'Superplasticizer (component 5)(kg in a m^3 mixture)',
                Aggregate = data$'Fine Aggregate (component 7)(kg in a m^3 mixture)',
                Age = data$'Age (day)')


d$S <-standardize(d$Strenght)
d$C <-standardize(d$Cement)
d$B <- standardize(d$Blast)
d$F <-standardize(d$Fly)
d$W <-standardize(d$Water )
d$P <- standardize(d$Plasticizer )
d$Age <-standardize(d$Age)
d$A <- standardize(d$Aggregate )


mC1.1 <- quap(
  alist(
    S ~ dnorm(mu, sigma),
    mu <- a + bC*C + bB*B + bF*F + bW*W + bP*P + bAge*Age + bA*A,
    a ~ dnorm(0, 0.2),
    bC ~ dnorm(0, 0.5),
    bB ~ dnorm(0, 0.5),
    bF ~ dnorm(0, 0.5),
    bW ~ dnorm(0, 0.5),
    bP ~ dnorm(0, 0.5),
    bAge ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)

mC2.1 <- quap(
  alist(
    S ~ dnorm(mu, sigma),
    mu <- a + bC*C,
    a ~ dnorm(0, 0.2),
    bC ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)

mC3.1 <- quap(
  alist(
    S ~ dnorm(mu, sigma),
    mu <- a + bA*A,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)

mC4.1 <- quap(
  alist(
    S ~ dnorm(mu, sigma),
    mu <- a + bA*A +  bC*C,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    bC ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)

plot(coeftab(mC2.1, mC3.1, mC4.1), par=c('bC','bA'))

plot(coeftab(mC1.1), par=c('bC','bA','bB', 'bW', 'bAge', 'bF', 'bP'))
set.seed(10)
prior <- extract.prior(mC1.1)
mu <- link(mC1.1, post=prior, data=list(C=c(-2,2)))
plot(NULL, xlim=c(-2,2), ylim=c(-2,2))
for (i in 1:50) lines(c(-2,2), mu[i,], col=col.alpha('black', 0.4)) 

C_seq <- seq(from=-3, to=3.2, length.out = 30)
mu <- link(mC1.1, data=list(C=C_seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(D~M, data=d, col=rangi2)
lines(M_seq, mu.mean, lwd=2)
shade(mu.PI, A_seq)
precis(m5.2)





mC2.1 <- quap(
  alist(
    ## C -> S <-B ; S <- F
    S ~ dnorm(mu, sigma),
    mu <- a + bC*C + bB*B + bF*F,
    a ~ dnorm(0, 0.2),
    bC ~ dnorm(0, 0.5),
    bB ~ dnorm(0, 0.5),
    bF ~ dnorm(0, 0.5),
    sigma ~ dexp(1),
    ## B --> C <-- F; 
    C ~ dnorm(mu_C, simga_C),
    mu_C <- aC + bBC*B + bFC*F,
    aC ~ dnorm(0, 0.2),
    bBC ~ dnorm(0, 0.5),
    bFC ~ dnorm(0, 0.5),
    simga_C ~ dexp(1),
    ## C --> F <-- B; 
    F ~ dnorm(mu_F, simga_F),
    mu_F <- aF + bCF*C + bBF*B,
    aF ~ dnorm(0, 0.2),
    bCF ~ dnorm(0, 0.5),
    bBF ~ dnorm(0, 0.5),
    simga_F ~ dexp(1)
    ## C --> B <-- F; 
    B ~ dnorm(mu_B, simga_B),
    mu_B <- aB + bCB*C + bFB*F,
    aF ~ dnorm(0, 0.2),
    bCB ~ dnorm(0, 0.5),
    bFB ~ dnorm(0, 0.5),
    simga_B ~ dexp(1)
  ), data=d
)

precis(mC1.1)
precis(mC2.1)



plot(coeftab(mC2.1), par=c('bC', 'bB', 'bF', 'bBC', 'bFC', 'bCF' , 'bBF'))

     
kidiq <- haven::read_dta("http://www.stat.columbia.edu/~gelman/arm/examples/child.iq/kidiq.dta")
library(rethinking)
library(dplyr)
head(kidiq)

idiq100 <- kidiq %>%
  mutate(mom_iq = mom_iq/100,
         kid_score = kid_score/100,
         mom_iq_c = mom_iq - 1,
         mom_hs = factor(mom_hs, labels= c('no', 'yes')))

d <- sample_n(idiq100, 7)

# R code 7.7
m7.2 <- quap(
  alist(
    kid_score ~ dnorm( mu , exp(log_sigma) ),
    mu <- a + b[1]*mom_iq + b[2]*mom_iq^2,
    a ~ dnorm( 90 , 20 ),
    b ~ dnorm( 0 , 10 ),
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,2)) )

## R code 7.8
m7.3 <- quap(
  alist(
    kid_score ~ dnorm( mu , exp(log_sigma) ),
    mu <- a + b[1]*mom_iq + b[2]*mom_iq^2 +
      b[3]*mom_iq^3,
    a ~ dnorm( 90 , 20 ),
    b ~ dnorm( 0 , 10 ),
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,3)) )

m7.4 <- quap(
  alist(
    kid_score ~ dnorm( mu , exp(log_sigma) ),
    mu <- a + b[1]*mom_iq + b[2]*mom_iq^2 +
      b[3]*mom_iq^3 + b[4]*mom_iq^4,
    a ~ dnorm( 90 , 20 ),
    b ~ dnorm( 0 , 10 ),
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,4)) )

m7.5 <- quap(
  alist(
    kid_score ~ dnorm( mu , exp(log_sigma) ),
    mu <- a + b[1]*mom_iq + b[2]*mom_iq^2 +
      b[3]*mom_iq^3 + b[4]*mom_iq^4 +
      b[5]*mom_iq^5,
    a ~ dnorm( 90 , 20 ),
    b ~ dnorm( 0 , 10 ),
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,5)) )

## R code 7.9
m7.6 <- quap(
  alist(
    kid_score ~ dnorm( mu , 0.001 ),
    mu <- a + b[1]*mom_iq + b[2]*mom_iq^2 +
      b[3]*mom_iq^3 + b[4]*mom_iq^4 +
      b[5]*mom_iq^5 + b[6]*mom_iq^6,
    a ~ dnorm( 90 , 20 ),
    b ~ dnorm( 0 , 10 )
  ), data=d , start=list(b=rep(0,6)) )

## R code 7.10
post <- extract.samples(m7.2)
mass_seq <- seq( from=min(d$mom_iq) , to=max(d$mom_iq) , length.out=100 )
l <- link( m7.2 , data=list( mom_iq=mass_seq ) )
mu <- apply( l , 2 , mean )
ci <- apply( l , 2 , PI )
plot( kid_score ~ mom_iq , data=d )
lines( mass_seq , mu )
shade( ci , mass_seq )

post <- extract.samples(m7.3)
mass_seq <- seq( from=min(d$mom_iq) , to=max(d$mom_iq) , length.out=100 )
l <- link( m7.3 , data=list( mom_iq=mass_seq ) )
mu <- apply( l , 2 , mean )
ci <- apply( l , 2 , PI )
plot( kid_score ~ mom_iq , data=d )
lines( mass_seq , mu )
shade( ci , mass_seq )

post <- extract.samples(m7.4)
mass_seq <- seq( from=min(d$mom_iq) , to=max(d$mom_iq) , length.out=100 )
l <- link( m7.4 , data=list( mom_iq=mass_seq ) )
mu <- apply( l , 2 , mean )
ci <- apply( l , 2 , PI )
plot( kid_score ~ mom_iq , data=d )
lines( mass_seq , mu )
shade( ci , mass_seq )


post <- extract.samples(m7.5)
mass_seq <- seq( from=min(d$mom_iq) , to=max(d$mom_iq) , length.out=100 )
l <- link( m7.5 , data=list( mom_iq=mass_seq ) )
mu <- apply( l , 2 , mean )
ci <- apply( l , 2 , PI )
plot( kid_score ~ mom_iq , data=d )
lines( mass_seq , mu )
shade( ci , mass_seq )


post <- extract.samples(m7.6)
mass_seq <- seq( from=min(d$mom_iq) , to=max(d$mom_iq) , length.out=100 )
l <- link( m7.6 , data=list( mom_iq=mass_seq ) )
mu <- apply( l , 2 , mean )
ci <- apply( l , 2 , PI )
plot( kid_score ~ mom_iq , data=d )
lines( mass_seq , mu )
shade( ci , mass_seq )


# Function to compute DIC
dic_brmsfit <- function(object) {
  Dbar <- -2 * mean(rowSums(log_lik(object)))
  coef_pmean <- unname(fixef(m1)[ , "Estimate"])
  X <- model.matrix(as.formula(object$formula), object$data)
  res <- res <- residuals(m1)[ , "Estimate"]
  N <- length(res)
  sigma <- posterior_summary(m1, pars = "sigma")[ , "Estimate"]
  Dhat <- -2 * sum(dnorm(res, sd = sigma, log = TRUE))
  p <- Dbar - Dhat
  elpd <- Dhat / -2 - p
  data.frame(elpd_dic = elpd, p_dic = p, dic = Dhat + 2 * p, 
             row.names = "Estimate")
}
dic_brmsfit(m1)

set.seed(1)
sum(lppd(m7.2, n=1e4))

sapply(list(m7.2, m7.3, m7.4, m7.5, m7.6), function(m) sum(lppd(m)) )

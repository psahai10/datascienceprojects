#https://towardsdatascience.com/estimating-probabilities-with-bayesian-modeling-in-python-7144be007815

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


library('pscl')
library('dplyr')

d <- UKHouseOfCommons

uk92_data <- within(list(), {
  y <- as.matrix(dplyr::select(UKHouseOfCommons, y1, y2))
  X <- model.matrix(~ 0 + y1lag + y2lag + coninc + labinc + libinc, data = UKHouseOfCommons) %>% scale()
  N <- nrow(y)
  K <- ncol(y)
  P <- ncol(X)
  alpha_loc <- rep(0, K)
  alpha_scale <- rep(10, K)
  beta_loc <- matrix(0, K, P)
  beta_scale <- matrix(2.5, K, P)
  Sigma_corr_shape <- 2
  Sigma_scale_scale <- 5
})

glimpse(uk92_data)

install.packages(c('VGAM', 'AER', 'geepack', 'MCMCpack', 'Amelia', 'MatchIt', 'maxLik'))
install.packages("https://cran.r-project.org/src/contrib/Archive/Zelig/Zelig_5.0-12.tar.gz", repos=NULL, type="source")
library('Zelig')
data("turnout")
d <- turnout
d$racegroup <- ifelse(d$race=='white', 1, 2)
data_list <- list(
  voted = data$vote,
  age = (d$age-min(d$age))/mean(d$age),
  income = (d$income-min(d$income))/mean(d$income),
  race_group = d$racegroup,
  education = (d$educate-min(d$educate))/mean(d$educate))
)

library(dplyr)

data <- d %>% mutate(agegroup = case_when( age >= 50  ~ '3',
                                             age >= 30  & age < 50 ~ '2',
                                             age > 0  & age < 30 ~ '1'))

data <- data %>% mutate(incomegroup = case_when( income >= 8 ~ '3',
                                              income >=4 & income < 8 ~ '2',
                                              income < 4 ~ '1'))

data$racegroup <- ifelse(d$race=='white', 1, 2)


data <- data %>% mutate(educategroup = case_when( educate >= 14 ~ '3',
                                              educate >=10 & educate < 14 ~ '2',
                                              educate < 10 ~ '1'))



dat_list <- list(
  voted = data$vote,
  age_group = as.integer(data$agegroup),
  income_group = as.integer(data$incomegroup),
  race_group = data$racegroup,
  education = as.integer(data$educategroup)
)

mvot <- ulam(
  alist(
    voted ~ dbinom(1, p),
    logit(p) ~ a + b*income + c[race_group] + d*education + e*age,
    a ~ dnorm(0, 1.5),
    b ~ dnorm(0.5, 0.5),
    c[race_group] ~ dnorm(0, 0.5),
    d ~ dnorm(0.5, 0.5),
    e ~ dnorm(0.5, 0.5)
  ), data=data_list, chains=4, cores=4
)

mvot <- quap(
  alist(
    voted ~ dbinom(1, p),
    logit(p) ~ a + b*income + c[race_group] + d*education + e*age,
    a ~ dnorm(0, 1.5),
    b ~ dnorm(0.5, 0.5),
    c[race_group] ~ dnorm(0, 1.5),
    d ~ dnorm(0.5, 0.5),
    e ~ dnorm(0.5, 0.5)
  ), data=data_list
)

precis(mvot, depth=2)
post <- extract.samples(mvot)
vote <- inv_logit(post$a)
plot(precis(as.data.frame(vote)))


mvot <- quap(
  alist(
    voted ~ dbinom(1, p),
    logit(p) ~ a + b[income] + c[race_group] + d[education] + e[age],
    a ~ dnorm(1, 1.5),
    b[income] ~ dnorm(1, 1.5),
    c[race_group] ~ dnorm(1, 1.5),
    d[education] ~ dnorm(1.0, 1.5),
    e[age] ~ dnorm(1, 1.5)
  ), data=data_list
)

mvot <- quap(
  alist(
    voted ~ dbinom(1, p),
    logit(p) ~ a[age_group],
    a[age_group] ~ dnorm(1, 1.5)
  ), dat=dat_list)

precis(mvot, depth=2)
post <- extract.samples(mvot)
vote <- inv_logit(post$a)
plot(precis(as.data.frame(vote)))


library('foreign')
library('nnet')
library('ggplot2')
library('reshape2')
library('dplyr')
library('rethinking')
library('rstan')

d <- read.dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta")

d$prog2 <- relevel(d$prog, ref = "academic")


d <- d %>% mutate(session = case_when( ses == 'low'  ~ 1,
                                       ses == 'middle'  ~ 2,
                                       ses == 'high'  ~ 3))

d <- d %>% mutate(program = case_when( prog == 'vocation'  ~ 1,
                                       prog == 'general'  ~ 2,
                                       prog == 'academic'  ~ 3))

d$treatment <- d$session + 3*d$program - 3

dat_list <- list(N=200, K=9, treatment=d$treatment, write=d$write/mean(d$write))

dat_list <- list(N=200, K=3, treatment=d$program, writing_score=d$write/mean(d$write))

library('mlogit')



slim_d <- data.frame(key=d$id, ses=as.factor(d$ses), write=as.numeric(d$write), 
                     chid=d$cid, prog=as.factor(d$prog), read=as.numeric(d$read))
slim_d$id <- ifelse(slim_d$school=='vocation', 1, ifelse(slim_d$school=='academic', 2, 3))

slim_d$response <- with(slim_d, cbind(write, read))

brms_fit <- brm(response | trials(200) ~ (1|ses) + (1|prog), data = slim_d, family = multinomial())


idx <- rep(1:nrow(slim_d), 3)

df <- data.frame(school=list1,list2,list3)
dupdf <- slim_d[idx,]



list1 <- data.frame(schoolz = rep('vocation', 200))
list2 <- data.frame(schoolz= rep('academic', 200))
list3 <- data.frame(schoolz = rep('general', 200))
df <- rbind(list1, list2, list3)

dupdf <- cbind(dupdf, df)


dupdf$prog <- ifelse(dupdf$school==dupdf$schoolz, TRUE, FALSE)

d <- subset(dupdf, select=-schoolz)


programLong <- mlogit.data(d, alt.levels=c('vocation', 'academic', 'general'))


slim <- data.frame(id=as.integer(d$key), prog=as.numeric(d$prog), ses=as.factor(d$ses), 
                   write=as.numeric(d$write), chid=d$chid, alt=d$school)

slim <- slim[order(slim$id),]

programLong <- mlogit.data(slim, alt.levels=c('vocation', 'academic', 'general'))


multi_mod = mlogit(prog ~ 1|write + ses, data=programLong)
summary(multi_mod)
##
multinomregML = function(par, X, y) {
  levs = levels(y)
  ref = levs[1]               # reference level (category label 1)
  
  y0 = y==ref
  y1 = y==levs[2]             # category 2
  y2 = y==levs[3]             # category 3
  
  beta = matrix(par, ncol=2)
  
  # more like mnlogit function
  # V1 = X %*% beta[,1]
  # V2 = X %*% beta[,2]
  # ll = sum(-log(1 + exp(V1)+exp(V2))) + sum(V1[y1],V2[y2])
  
  V = X %*% beta               #  a vectorized approach
  baseProbVec = 1/(1 + rowSums(exp(V)))  # reference group probabilities
  
  loglik = sum(log(baseProbVec))  + crossprod(c(V), c(y1, y2))
  loglik
}

out = optim(runif(8, -.1, .1), multinomregML, X=model.matrix(prog ~ ses + write, data = slim_d),
            y=slim_d$prog, control=list(maxit=1000, reltol=1e-12, ndeps=rep(1e-8, 8),
                                         trace=T, fnscale=-1, type=3),
            method='BFGS')
out$par

cbind(out$par, coef(multi_mod)[c(1,5,7,3,2,6,8,4)])

X=model.matrix(prog ~ ses + write, data = slim_d)
y=slim_d$prog
y1 = y=='general'
y2 = y=='vocation'
pars = matrix(out$par, ncol=2)
V = X %*% pars
acadprob = 1/(1+rowSums(exp(V)))
fitnonacad = exp(V) * matrix(rep(acadprob, 2), ncol = 2)
fits = cbind(acadprob, fitnonacad)
yind = model.matrix(~-1+prog, data=slim_d)


# because dmultinom can't take matrix for prob
ll = 0
for (i in 1:200){
  ll = ll+dmultinom(yind[i,], size=1, prob=fits[i,], log=T)
}
ll

out$value

library(rethinking)
sppnames <- c('afarensis', 'africanus', 'habilis', 'boisei', 'rudolfensis', 'ergaster', 'sapiens')
brainvolcc <- c(438, 452, 612, 521, 752, 871, 1350)
masskg <- c(37.0, 35.5, 34.5, 41.5, 55.5, 61.0, 53.5)
d <- data.frame(species= sppnames, brain=brainvolcc, mass=masskg)

d$mass_std <- standardize(d$mass)

d$brain_std <- d$brain/max(d$brain)

m7.1 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b*mass_std,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), data=d
)
precis(m7.1)

mu <- link(m7.1)
str(mu)
mass.seq <- seq(from=-2, to=2, by=0.01)
mu <- link(m7.1, data=data.frame(mass_std=mass.seq))

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)
plot(brain_std ~ mass_std, data=d, col=col.alpha(rangi2, 0.1))
lines(mass.seq, mu.mean)
shade(mu.PI, mass.seq)

m7.1_OLS <- lm(brain_std ~ mass_std, data=d)
post <- extract.samples(m7.1_OLS)

mass.seq <- seq(from=-2, to=2, by=0.01)
mu <- link(m7.1_OLS, data=data.frame(mass_std=mass.seq))

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)
plot(brain_std ~ mass_std, data=d, col=col.alpha(rangi2, 0.1))
lines(mass.seq, mu.mean)
shade(mu.PI, mass.seq)


set.seed(12)
s <- sim(m7.1)
str(s)

mu_mean <- apply(s, 2, mean)
r <-  mu_mean - d$brain_std
resid_var <- var2(r)
outcome_var <- var2(d$brain_std)
1 - resid_var/outcome_var


s <- sim(m7.1_OLS)
str(s)
mu_mean <- apply(s, 2, mean)
r <-  mu_mean - d$brain_std
resid_var <- var2(r)
outcome_var <- var2(d$brain_std)
1 - resid_var/outcome_var

R2_is_bad <- function(quap_fit) {
  s <- sim(quap_fit, refresh=0)
  r <- apply(s,2,mean) - d$brain_std
  1 - var2(r)/var2(d$brain_std)
}

p <- c(0.3, 0.7)
-sum(p*log(p))

q <- c(0.25, 0.75)

sum(p*log(p/q))

set.seed(1)
lppd(m7.1, n=1e4)

logprob <- sim(m7.1, ll=TRUE, n=1e4)
n <- ncol(logprob)
ns <- nrow(logprob)
f <- function(x) log_sum_exp( logprob[,x]) - log(ns)
(lppd <- sapply(1:n, f))


library(rethinking)

N <- 20
kseq <- 1:5
dev <- sapply(kseq, function(k) {
  print(k);
  r <- replicate(1e4, sim_train_test(N=N, k=k));
  c( mean(r[1, ]), mean(r[2,]), sd(r[1,]), sd(r[2,]))
})

plot( 1:5, dev[1,], ylim=c( min(dev[1:2,])-5, max(dev[1:2,])+10 ),
      xlim=x(1,5.1), xlab='number of parameters', ylab='deviance',
      pch=16, col=rangi2)
mtext(concat('N= ',N))
points( (1:5)+0.1, dev[2,])

for (i in kseq) {
    pts_in <- dev[1, i] + c(-1, +1)*dev[3,i]
    pts_out <- dev[2, i] + c(-1,+1)*dev[4,i]
    lines( c(i,i), pts_in, col=rangi2)
    lines( c(i,i) +0.1, pts_out)
 }

data(cars)
m <- quap(
  alist(
    dist ~ dnorm(mu, sigma),
    mu <- a + b*speed,
    a ~ dnorm(0,100),
    b ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ), data=cars
)
set.seed(94)
post <- extract.samples(m, n=1e3)

data(cars)
m <- quap(
  alist(
    dist ~ dnorm(mu, sigma),
    mu <- a + b*speed,
    a ~ dnorm(0, 100),
    b ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ), data=cars
)
set.seed(94)

post <- extract.samples(m, n=1000)


n_samples <- 1000

logprob <- sapply(1:n_samples,
      function(s) {
        mu <- post$a[s] + post$b[s]*cars$speed
        dnorm( cars$dist, mu , post$sigma[s], log=TRUE )
      } )

str(logprob)

n_cases <- nrow(cars)

lppd <- sapply( 1:n_cases, function(i) log_sum_exp(logprob[i,]) - log(n_samples) )


pWAIC <- sapply(1:n_cases, function(i) var(logprob[i,]))

-2*(sum(lppd) - sum(pWAIC))

waic_vec <- -2*(lppd - pWAIC)
sqrt(n_cases*var(waic_vec))


set.seed(77)
WAIC(m6.7)
compare(m6.6, m6.7, m6.8, func=WAIC)
compare(m6.6, m6.7, m6.8, func=PSIS)

########### Calcularte WAIC from scratch #################

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
post <- extract.samples(m6.8, n=1000)

n_samples <- 300

f <- function(s) {
  mu <- post$a[s] + post$bt[s]*d$treatment
  dnorm(d$h1, mu, post$sigma[s], log=TRUE)
}

log_prob <- sapply(1:n_samples, f) 

log_prob <- sapply(1:n_samples, 
      function(s) {
          mu <- post$a[s] + post$bt[s]*d$treatment
          dnorm(d$h1, mu, post$sigma[s], log=TRUE)
     })

str(log_prob)

n_cases <- nrow(d)
f <- function(i) log_sum_exp(logprob[,i]) - log(n_samples)
lppd <- sapply(1:n_cases, f)

f <- function(i) var(logprob[,i])

pWAIC <- sapply(1:n_cases, f)

-2*(sum(lppd) - sum(pWAIC))

set.seed(91)
waic_m6.6 <- WAIC(m6.6, pointwise = TRUE)$WAIC
waic_m6.7 <- WAIC(m6.7, pointwise = TRUE)$WAIC
waic_m6.8 <- WAIC(m6.8, pointwise = TRUE)$WAIC

n<- length(waic_m6.7)
diff_m6.7_m6.8 <- waic_m6.7 - waic_m6.8

sqrt(n*var(diff_m6.7_m6.8))
40.0 + c(-1,1)*10.4*2.6

plot(compare(m6.6,m6.7,m6.8))

plot(compare(m6.6,m6.7,m6.8, func=PSIS))

set.seed(92)
waic_m6.6 <- WAIC(m6.6, pointwise = TRUE)$WAIC
waic_m6.8 <- WAIC(m6.8, pointwise = TRUE)$WAIC

n<- length(waic_m6.6)
diff_m6.6_m6.8 <- waic_m6.6 - waic_m6.8
sqrt(n*var(diff_m6.6_m6.8))

compare(m6.6,m6.7,m6.8)@dSE

library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce
d$A <- standardize( d$MedianAgeMarriage )
d$D <- standardize( d$Divorce )
d$M <- standardize( d$Marriage )

m5.1 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bA * A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )

m5.2 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM * M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )

m5.3 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )

set.seed(24071847)
compare( m5.1 , m5.2 , m5.3 , func=PSIS )
plot(compare( m5.1 , m5.2 , m5.3 , func=PSIS ))

PSIS_m5.3 <- PSIS(m5.3, pointwise = TRUE)
WAIC_m5.3 <- WAIC(m5.3, pointwise = TRUE)

plot(PSIS_m5.3$k, WAIC_m5.3$penalty, xlab='PSIS Pareto k',
     ylab='WAIC penalty', col=rangi2, lwd=2)

WAIC_m5.1 <- WAIC(m5.1, pointwise = TRUE)

x <- d$M - d$A
x <- x - min(x)
x <- x / max(x)
# draw the plot
dif_waic <- WAIC_m5.1 - WAIC_m5.3



plot( dif_waic$WAIC , d$D ,
      xlab="pointwise difference in WAIC" , ylab="Divorce Rate(std)" , pch=21 ,
      col=col.alpha("black",0.8) , cex=1+x , lwd=2 , bg=col.alpha(rangi2,0.4) )

abline( v=0 , lty=2 )
abline( h=0 , lty=2 )


df <- data.frame(dif_waic= dif_waic$WAIC, D=d$D, x = x, name=d$Location)
library(ggplot2)
library(ggrepel)

joint_score_plot <- ggplot(df,
                           aes(x = dif_waic, y = D, label=name)) +
  geom_point(aes(color = 'red'), alpha = 0.6, size = x*8) + xlim(-0.3, 0.3) + ylim(-1.5,1.5) +
  labs(x = "pointwise difference in WAIC", y = "Divorce Rate", title = "Figure 7.1") +
  #color = "TBD",
  geom_text_repel(data=subset(df,D>0.7 | D<(-0.9) | dif_waic>0.15| dif_waic<(-0.1))) +
  #geom_text(aes(label=ifelse(df$name %in% species.to.label, as.character(name),'')),hjust=0,vjust=1) +
  #geom_text(aes(label=ifelse(log_L>1.2,as.character(name),'')),hjust=0,vjust=1) +   
  #geom_text(aes(label=ifelse(dif_waic.WAIC>0.2,as.character(name),'')),hjust=0,vjust=1) +  
  #geom_text(aes(label=ifelse(log_L<(-1.25),as.character(name),'')),hjust=0,vjust=1) +   
  #geom_text(aes(label=ifelse(dif_waic.WAIC<(-0.2),as.character(name),'')),hjust=0,vjust=1) +  
  #geom_text(aes(label=ifelse(log_L<-1.2,as.character(name),'')),hjust=0,vjust=1) + 
  # Make the color scale go from darkorange (bad) to midpoint gray (usually 
  # associate white with missing) to darkblue (good):
  #scale_color_gradient2(low = "darkorange", mid = "gray", high = "darkblue",
  #midpoint = 0) +
  # Add dashed red lines through the origin:
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  # Add some text for the different quadrants for ease of reading:
  annotate("text", label = "Large Body +\n Long Life", x = .25, y = .25,
           size = 5, color = "darkred") +
  annotate("text", label = "Large Body +\n Short Life", x = .25, y = -.20,
           size = 5, color = "darkred") +
  annotate("text", label = "Large Body +\n Long Life", x = -.25, y = .25,
           size = 5, color = "darkred") +
  annotate("text", label = "Large Body +\n Short Life", x = -.25, y = -.20,
           size = 5, color = "darkred") +
  annotate("text", label = "<-- m5.1 better", x = -.1, y = -1.50,
           size = 5, color = "black") +
  annotate("text", label = "m5.3 better -->", x = .1, y = -1.50,
           size = 5, color = "black") +
  
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16))


joint_score_plot



m5.1t <- quap(
  alist(
    D ~ dstudent( 2, mu , sigma ) ,
    mu <- a + bA * A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )

m5.2t <- quap(
  alist(
    D ~ dstudent(2,  mu , sigma ) ,
    mu <- a + bM * M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )


m5.3t <- quap(
  alist(
    D ~ dstudent(2,  mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )


set.seed(24071847)

PSIS_m5.3t <- PSIS(m5.3t, pointwise = TRUE)
WAIC_m5.3t <- WAIC(m5.3t, pointwise = TRUE)
plot(PSIS_m5.3t$k, WAIC_m5.3t$penalty, xlab='PSIS Pareto k',
     ylab='WAIC penalty', col=rangi2, lwd=2)

plot(compare(m5.1,  m5.2 , m5.3, func=WAIC ))
plot(compare(m5.1t , m5.2t , m5.3t , func=WAIC ))

precis(m5.3)
precis(m5.3t)

WAIC(m5.3)
WAIC(m5.3t)

compare(m5.3,  m5.3t)
plot(compare(m5.3,  m5.3t))


library(rethinking)
data("Primates301")
d <- Primates301
head(d)

d$log_L <- scale( log(d$longevity) )
d$log_B <- scale( log(d$brain) )
d$log_M <- scale( log(d$body) )


sapply( d[,c("log_L","log_B","log_M")] , function(x) sum(is.na(x)) )

d2 <- d[ complete.cases( d$log_L , d$log_M , d$log_B ) , ]
nrow(d2)

m7.8 <- quap(
  alist(
    log_L ~ dnorm( mu , sigma ),
    mu <- a + bM*log_M + bB*log_B,
    a ~ dnorm(0,0.1),
    bM ~ dnorm(0,0.5),
    bB ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ) , data=d2 )

m7.9 <- quap(
  alist(
    log_L ~ dnorm( mu , sigma ),
    mu <- a + bB*log_B,
    a ~ dnorm(0,0.1),
    bB ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ) , data=d2 )
m7.10 <- quap(
  alist(
    log_L ~ dnorm( mu , sigma ),
    mu <- a + bM*log_M,
    a ~ dnorm(0,0.1),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ) , data=d2 )


set.seed(301)
compare( m7.8 , m7.9 , m7.10 )

plot( compare( m7.8 , m7.9 , m7.10 ) )

plot( coeftab( m7.8 , m7.9 , m7.10 ) , pars=c("bM","bB") )


cor( d2$log_B , d2$log_M )

plot(d2$log_B ~ d2$log_M)

#post <- extract.samples(m7.8)
#post 
#library(dplyr)

#dd <- sample_n(post, 800)
#plot(dd$bB ~ dd$bM)

waic_m7.8 <- WAIC( m7.8 , pointwise=TRUE )
waic_m7.9 <- WAIC( m7.9 , pointwise=TRUE )

# compute point scaling
x <- d2$log_B - d2$log_M
x <- x - min(x)
x <- x / max(x)
# draw the plot
dif_waic <- waic_m7.8 - waic_m7.9

df <- data.frame(dif_waic= dif_waic, log_L=d2$log_L, x = x, name=d2$name)


plot( dif_waic$WAIC , d2$log_L ,
      xlab="pointwise difference in WAIC" , ylab="log longevity (std)" , pch=21 ,
      col=col.alpha("black",0.8) , cex=1+x , lwd=2 , bg=col.alpha(rangi2,0.4) )

abline( v=0 , lty=2 )
abline( h=0 , lty=2 )




m7.11 <- quap(
  alist( 
    log_B ~ dnorm( mu , sigma ), 
    mu <- a + bM*log_M + bL*log_L, 
    a ~ dnorm(0,0.1), 
    bM ~ dnorm(0,0.5), 
    bL ~ dnorm(0,0.5), 
    sigma ~ dexp(1) 
    ) , data=d2
  ) 

precis( m7.11 )

m7.12 <- quap(
  alist( 
    log_B ~ dnorm( mu , sigma ), 
    mu <- a + bM*log_M , 
    a ~ dnorm(0,0.1), 
    bM ~ dnorm(0,0.5), 
    sigma ~ dexp(1) 
  ) , data=d2
) 
precis( m7.12 )

m7.13 <- quap(
  alist( 
    log_B ~ dnorm( mu , sigma ), 
    mu <- a + bL*log_L, 
    a ~ dnorm(0,0.1), 
    bL ~ dnorm(0,0.5), 
    sigma ~ dexp(1) 
  ) , data=d2
)
precis( m7.13 )

plot( coeftab(m7.11, m7.12, m7.13) , pars=c("bM", 'bL') )





library(ggplot2)
library(ggrepel)

# First generate the scattterplot with the points
joint_score_plot <- ggplot(df,
                           aes(x = dif_waic.WAIC, y = log_L, label=name)) +
  geom_point(aes(color = 'red'), alpha = 0.6, size = x*8) + xlim(-0.3, 0.3) + ylim(-1.5,1.5) +
  labs(x = "pointwise difference in WAIC", y = "log longevity", title = "Figure 7.1") +
       #color = "TBD",
  #geom_text_repel(data=subset(df,log_L>1.2 | log_L<(-1.25) | dif_waic.WAIC>0.2 | dif_waic.WAIC<(-0.2))) +
  geom_text(aes(label=ifelse(df$name %in% species.to.label, as.character(name),'')),hjust=0,vjust=1) +
  #geom_text(aes(label=ifelse(log_L>1.2,as.character(name),'')),hjust=0,vjust=1) +   
  #geom_text(aes(label=ifelse(dif_waic.WAIC>0.2,as.character(name),'')),hjust=0,vjust=1) +  
  #geom_text(aes(label=ifelse(log_L<(-1.25),as.character(name),'')),hjust=0,vjust=1) +   
  #geom_text(aes(label=ifelse(dif_waic.WAIC<(-0.2),as.character(name),'')),hjust=0,vjust=1) +  
  #geom_text(aes(label=ifelse(log_L<-1.2,as.character(name),'')),hjust=0,vjust=1) + 
  # Make the color scale go from darkorange (bad) to midpoint gray (usually 
  # associate white with missing) to darkblue (good):
  #scale_color_gradient2(low = "darkorange", mid = "gray", high = "darkblue",
                        #midpoint = 0) +
  # Add dashed red lines through the origin:
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  # Add some text for the different quadrants for ease of reading:
  annotate("text", label = "Large Body +\n Long Life", x = .25, y = .25,
           size = 5, color = "darkred") +
  annotate("text", label = "Large Body +\n Short Life", x = .25, y = -.20,
           size = 5, color = "darkred") +
  annotate("text", label = "Large Body +\n Long Life", x = -.25, y = .25,
           size = 5, color = "darkred") +
  annotate("text", label = "Large Body +\n Short Life", x = -.25, y = -.20,
           size = 5, color = "darkred") +
  annotate("text", label = "<-- m7.8 better", x = -.1, y = -1.50,
           size = 5, color = "black") +
  annotate("text", label = "m7.9 black -->", x = .1, y = -1.50,
           size = 5, color = "black") +
  
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16))

species.to.label <- c('Cebus_capucinus', 'Cebus_apella', 'Cebus_olivaceus', 'Gorilla_gorilla_gorilla', 'Eulemur_fulvus',
                      'Lepilemur_leucopus', 'Cercopithecus_lhoesti', 'Cacajao_melanocephalus')
joint_score_plot
# Make passing marginal along x axis:
pass_dens <- axis_canvas(joint_score_plot, axis = "x") +
  geom_density(data = team_summary_df, aes(x = pass_epa_diff), 
               fill = "darkblue",
               alpha = 0.5, size = .2)

# Same thing for rushing but along y and with coord_flip = TRUE to make it vertical:
rush_dens <- axis_canvas(joint_score_plot, axis = "y", coord_flip = TRUE) +
  geom_density(data = team_summary_df, aes(x = rush_epa_diff), 
               fill = "darkblue",
               alpha = 0.5, size = .2) +
  coord_flip()

# Now generate by adding these objects to the main plot:
joint_score_plot_1 <- insert_xaxis_grob(joint_score_plot, pass_dens, 
                                        grid::unit(.2, "null"), position = "top")
joint_score_plot_2 <- insert_yaxis_grob(joint_score_plot_1, rush_dens, 
                                        grid::unit(.2, "null"), position = "right")
ggdraw(joint_score_plot_2)


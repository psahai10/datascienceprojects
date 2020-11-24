#link website = http://www.mas.ncl.ac.uk/~nmf16/teaching/mas3301/

d1 <- data.frame(chol = c(233, 291, 312, 250, 246, 197, 268, 224, 239, 239, 
                          254, 276, 234, 181, 248, 252, 202, 218, 212, 325),
                 category = 'A')
                 
d2 <- data.frame(chol = c(344, 185, 263, 246, 224, 212, 188, 250, 148, 169, 
                          226, 175, 242, 252, 153, 183, 137, 202, 194, 213),
                 category = 'B')

d <- rbind(d1, d2)
levels(d$category)

d$cid <- ifelse(d$category=='A', 1, 2)
d$chol_std <- d$chol/max(d$chol)

library(rethinking)

chol <- quap(
  alist(
    chol ~ dnorm( mu[cid], sigma),
    mu[cid] ~ dnorm( 230, 20),
    sigma ~ dunif(0, 50)
  ), data=d
)


precis(chol, depth=2)
vcov(chol)

d$cid <- as.integer(d$cid)

labels <- paste('mu[', 1:2, ']:', levels(d$cid), sep='')

plot(precis(chol, depth=2, pars='mu'), labels=labels)

post <- extract.samples(chol)
post$diff_cat <- post$mu[,1] - post$mu[,2]
head(post)

precis(post, depth=2, hist=FALSE)

mean(post$diff_cat )

quantile(post$mu[,1], c(0.1, 0.9))
quantile(post$mu[,2], c(0.1, 0.9))
quantile(post$diff_cat, c(0.1, 0.9))


PI(post$mu[,1], prob=0.9)
PI(post$mu[,2], prob=0.9)
PI(post$diff_cat, prob=0.9)









d <- rbind(d1, d2)
levels(d$category)

d$cid <- ifelse(d$category=='A', 1, 2)
d$chol_std <- d$chol/max(d$chol)

m8.1s <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a + b*(rugged_std-0.215),
    a ~ dnorm(1, 1),
    b ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=dd
)

sample_mu <- rnorm(1e4, )



set.seed(7)
prior <- extract.prior(m8.1)

plot(NULL, xlim=c(0,1), ylim=c(0.5, 1.5),
     xlab='ruggedness', ylab='log GDP' )
abline( h=c(min(dd$log_gdp_std), max(dd$log_gdp_std)), col=c('black', 'red'), lty=c(2, 2))

rugged_seq <- seq(from=-0.1, to=1.1, length.out =30)
mu <- link(m8.1s, post=prior, data=data.frame(rugged_std=rugged_seq))
for (i in 1:50) lines(rugged_seq, mu[i,] , col=col.alpha('black', 0.3))

sum(abs(prior$b) > 0.6) / length(prior$b)


t1<-1/2500
t2<-1/(20^2)
starta<-c(0.1,5)
tau<-c(t1,t2)
p<-c(0.5,0.9)


ab<-findab(starta,tau,p,20)
ab

cholest<-scan("http://www.mas.ncl.ac.uk/~nmf16/teaching/mas3301/cholesterol.txt")


group<-rep(1:2,c(20,20))
d0<-2*ab[1]
v0<-2*ab[2]/d0
C0<-matrix(c(2,-1,-1,2),ncol=2)*2/3
M0<-c(200,200)

post<-oneway(d0,v0,C0,M0,cholest,group)
post

w1<-solve(post$C1)
stepA <- 1
muA <- seq(180,280,stepA)
seA <- sqrt(post$v1*w1[1,1])
tA <- (muA-post$M1[1])/seA[,1]
densA<-dt(tA,post$d1)/seA[,1]
plot(muA,densA,type="l",xlab=expression(mu[A]),ylab="Density")

stepB<-1
muB<-seq(180,280,stepB)
seB<-sqrt(post$v1*w1[2,2])
tB<-(muB-post$M1[2])/seB[,1]
densB<-dt(tB,post$d1)/seB[,1]
plot(muB,densB,type="l",xlab=expression(mu[B]),ylab="Density")

mA<-matrix(muA,nrow=length(muA),ncol=length(muB))
mB<-matrix(muB,nrow=length(muA),ncol=length(muB),byrow=T)
dA<-mA-post$M1[1]
dB<-mB-post$M1[2]
B1<-post$d1*post$v1
q<-(dA^2)*post$C1[1,1] + (dB^2)*post$C1[2,2] + 2*dA*dB*post$C1[1,2]

dens<-(B1+q)^(-(post$d1+2)/2)
dens<-dens/(sum(dens)*stepA*stepB)
contour(muA,muB,dens)
contour(muA,muB,dens,xlab=expression(mu[A]),ylab=expression(mu[B]))

abline(0,1,lty=2)








findab<-function(starta,tau,p,niter)
{aa<-starta
a<-numeric(niter)
b<-numeric(niter)
diff<-numeric(niter)
c1<-findb(aa[1],tau,p)
c2<-findb(aa[2],tau,p)
while((c1[2]*c2[2])>0)
{aa[1]<-aa[1]/2
aa[2]<-aa[2]*2
c1<-findb(aa[1],tau,p)
c2<-findb(aa[2],tau,p)
}
for (i in 1:niter)
{a[i]<-(aa[1]+aa[2])/2
newc<-findb(a[i],tau,p)
diff[i]<-newc[2]
b[i]<-newc[1]
if ((diff[i]*c1[2])>0)
{aa[1]<-a[i]
c1<-newc
}
else
{aa[2]<-a[i]
c2<-newc
}
}
a<-signif(a,4)
b<-signif(b,4)
diff<-signif(diff,4)
table<-data.frame(a,b,diff)
write.table(table,file="")
c(a[niter],b[niter])
}


findb<-function(a,tau,p)
{q<-qgamma(p[2],a,1)
b<-q/tau[2]
diff<-pgamma(tau[1],a,b)-p[1]
c(b,diff)
}

oneway<-function(d0,v0,C0,M0,y,group)
{J<-max(group)
N<-numeric(J)
ybar<-numeric(J)
Sd<-0
C1<-C0
M1<-C0%*%M0
for (j in 1:J)
{N[j]<-sum(group==j)
yj<-y[group==j]
ybar[j]<-mean(yj)
Sd<-Sd+sum((yj-ybar[j])^2)
C1[j,j]<-C1[j,j]+N[j]
}
M1<-M1+N*ybar
M1<-solve(C1,M1)
R<-t(M0)%*%C0%*%M0+sum(N*ybar*ybar)-t(M1)%*%C1%*%M1
N<-sum(N)
Nvd<-Sd+R
d1<-d0+N
v1<-(d0*v0+Nvd)/d1
list(d1=d1,v1=v1,C1=C1,M1=M1)
}

oneway<-function(d0,v0,C0,M0,y,group)
{J<-max(group)
N<-numeric(J)
ybar<-numeric(J)
Sd<-0
C1<-C0
M1<-C0%*%M0
for (j in 1:J)
{N[j]<-sum(group==j)
yj<-y[group==j]
ybar[j]<-mean(yj)
Sd<-Sd+sum((yj-ybar[j])^2)
C1[j,j]<-C1[j,j]+N[j]
}
M1<-M1+N*ybar
M1<-solve(C1,M1)
R<-t(M0)%*%C0%*%M0+sum(N*ybar*ybar)-t(M1)%*%C1%*%M1
N<-sum(N)
Nvd<-Sd+R
d1<-d0+N
v1<-(d0*v0+Nvd)/d1
list(d1=d1,v1=v1,C1=C1,M1=M1)
}

library(rethinking)
a <- 3.5
b <- (-1)
sigma_a <- 1
sigma_b <- 0.5
rho <- (-0.7)

Mu <- c(a,b)

cov_ab <- sigma_a*sigma_b*rho

Sigma <- matrix(c(sigma_a^2, cov_ab, cov_ab, sigma_b^2), ncol=2)

sigmas <- c(sigma_a, sigma_b)

Rho <- matrix( c(1,rho,rho,1), nrow=2)

Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

N_cafes <- 20

library(MASS)
set.seed(5)
vary_effects <- mvrnorm(N_cafes, Mu, Sigma)

a_cafe <- vary_effects[,1]
b_cafe <- vary_effects[,2]

plot(a_cafe, b_cafe, col=rangi2)

library(ellipse)
for (l in c(0.1,0.3,0.5,0.8,0.99))
  lines(ellipse(Sigma, centre=Mu, level=l), col=col.alpha('black', 0.4))

set.seed(22)
N_visits <- 10
afternoon <- rep(0:1, N_visits*N_cafes/2)
cafe_id <- rep(1:N_cafes, each=N_visits)
mu <- a_cafe[cafe_id] + b_cafe[cafe_id]*afternoon
sigma <- 0.5
wait <- rnorm(N_visits*N_cafes, mu, sigma)
d <- data.frame(cafe=cafe_id, afternoon=afternoon, wait=wait)







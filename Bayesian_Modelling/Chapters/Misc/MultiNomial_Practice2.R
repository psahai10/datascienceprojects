# Poissin_Multinomial
library(Ecdat)
data(Yogurt)
d <- Yogurt

data_wide <- data.frame(obs=c(1,2,3,4), X1=c(0,0,1,1), X2=c(0,1,0,1), Y1= c(3,5,7,1), 
                        Y2=c(5,5,2,3), Y3=c(2,0,1,6), N=c(10,10,10,10))

data_wide$response <- with(data_wide, cbind(Y1, Y2, Y3))

brms_fit <- brm(response | trials(N) ~ (1|X1) + (1|X2), data = data_wide, family = multinomial())


dummy_data <- data.frame(obs=c(1,2,3,4), X1=c(0,0,1,1), X2=c(0,1,0,1), Y1= c(3,5,7,1), Y2=c(5,5,2,3), Y3=c(2,0,1,6))

library(tidyr)
data_long <- gather(dummy_data, C, Y, Y1:Y3, factor_key=TRUE)
data_long <- data_long[order(data_long$obs), ]
rownames(data_long) <- 1:nrow(data_long)
data_long$C <- ifelse(data_long$C=='Y1', 1, ifelse(data_long$C=='Y2', 2,3))
data_long$C <- as.factor(data_long$C )
data_long$obs <- as.factor(data_long$obs)

yogurt_wide <- d

library(dplyr)

names(d)
yogurt_wide <- yogurt_wide %>% 
  rename(
    feat.yoplait = yoplait,
    feat.dannon = dannon,
    feat.hiland = hiland,
    feat.weight = weight,
    price.yoplait = py,
    price.dannon = pd,
    price.hiland = ph,
    price.weight = pw
  )


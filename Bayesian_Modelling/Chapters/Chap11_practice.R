getwd()
setwd("C:/Users/psahai/Documents/R/R_Bayesian_Statistical_Rethinking/practice")
getwd()
 
d <- read.csv('all_data.csv')
head(d, 20)
tail(d)
names(d)
str(d)
levels(d$hormone)
levels(d$right.left)
levels(d$partner)

d$side_coding <- ifelse(d$right.left=='L', 1, 2)

No_hormone_data_set <- d[is.na(d$hormone), ]
hormone_data_set <- d[!is.na(d$hormone), ]
levels(hormone_data_set$hormone)

hormone_data_set$hormone_coding <- ifelse(hormone_data_set$hormone=='A',)
library(plyr)

hormone_data_set$hormone_coding <- revalue(hormone_data_set$hormone, c('A'=1, 'B'=2, 'C'=3, 'SAL'=4))
hormone_data_set

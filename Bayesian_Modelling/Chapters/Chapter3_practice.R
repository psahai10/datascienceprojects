# Practice Salinity vs water temp

getwd()
setwd("C:/Users/psahai/Documents/R/R_Bayesian_Statistical_Rethinking/practice")
getwd()
data <- read.csv("bottle.csv")
head(data)
library(psych)
psych::describe(data)
names(data)
df <- data.frame(data)
df2 <- data.frame(df$T_degC, df$Salnty)
df2
library(dplyr)

convert_temp <- function(x) {
  temp_F <- (x*9/5) + 32
  return(tempF)
}

df<- df2 %>% as_tibble() %>% mutate(
  df.T_degF = df.T_degC*9/5 + 32
)

df <- data.frame(Temp=df$df.T_degF, Salinity=df$df.Salnty)
df <- na.omit(df)

psych::describe(df)





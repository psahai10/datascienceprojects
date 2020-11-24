# https://www.rensvandeschoot.com/tutorials/brms-started/

library(brms) # for the analysis
library(haven) # to load the SPSS .sav file
library(tidyverse) # needed for data manipulation.
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(ggmcmc)
library(ggthemes)
library(ggridges)

popular2data <- read_sav(file = "https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book/blob/master/chapter%202/popularity/SPSS/popular2.sav?raw=true")

popular2data <- select(popular2data, pupil, class, extrav, sex, texp, popular) # we select just the variables we will use
head(popular2data) # we have a look at the first 6 observations

ggplot(data  = popular2data,
       aes(x = extrav,
           y = popular))+
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter")+# to add some random noise for plotting purposes
  theme_minimal()+
  labs(title = "Popularity vs. Extraversion")

ggplot(data  = popular2data,
       aes(x = extrav,
           y = popular))+
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  geom_smooth(method = lm,
              se     = FALSE, 
              col    = "black",
              size   = .5, 
              alpha  = .8)+ # to add regression line
  theme_minimal()+
  labs(title    = "Popularity vs. Extraversion",
       subtitle = "add regression line")

ggplot(data    = popular2data,
       aes(x   = extrav,
           y   = popular,
           col = class))+ #to add the colours for different classes
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(100))+
  labs(title    = "Popularity vs. Extraversion",
       subtitle = "add colours for different classes")

ggplot(data      = popular2data,
      aes(x     = extrav,
          y     = popular,
          col   = class,
          group = class))+ #to add the colours for different classes
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(100))+
  geom_smooth(method = lm,
              se     = FALSE,
              size   = .5, 
              alpha  = .8)+ # to add regression line
  labs(title    = "Popularity vs. Extraversion",
       subtitle = "add colours for different classes and regression lines")

# To colour code the extremes, we need to write a small function that calculates the regression lines and adds a collumn indicating which clusters have the most extreme.
f1 <- function(data, x, y, grouping, n.highest = 3, n.lowest = 3){
  groupinglevel <- data[,grouping]
  res           <- data.frame(coef = rep(NA, length(unique(groupinglevel))), group = unique(groupinglevel))
  names(res)    <- c("coef", grouping)
  for(i in 1:length(unique(groupinglevel))){
    data2    <- as.data.frame(data[data[,grouping] == i,])
    res[i,1] <- as.numeric(lm(data2[, y] ~ data2[, x])$coefficients[2])
  }
  top    <- res %>% top_n(n.highest, coef)
  bottom <- res %>% top_n(-n.lowest, coef)
  res    <- res %>% mutate(high_and_low = ifelse(coef %in% top$coef, "top",  ifelse(coef %in% bottom$coef, "bottom", "none")))
  data3  <- left_join(data, res)
  return(data3)
}

f1(data = as.data.frame(popular2data), 
  x    = "extrav",
  y    = "popular",
  grouping = "class",
  n.highest = 3, 
  n.lowest = 3) %>%
  ggplot()+
  geom_point(aes(x     = extrav,
                 y     = popular, 
                 fill  = class, 
                 group = class),
             size     =  1, 
             alpha    = .5, 
             position = "jitter", 
             shape    = 21, 
             col      = "white")+
  geom_smooth(aes(x     = extrav,
                  y     = popular,
                  col   = high_and_low,
                  group = class,
                  size  = as.factor(high_and_low),
                  alpha = as.factor(high_and_low)),
              method = lm,
              se     = FALSE)+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradientn(colours = rainbow(100))+
  scale_color_manual(values=c("top"      = "blue",
                              "bottom"   = "red",
                              "none"     = "grey40"))+
  scale_size_manual(values=c("top"       = 1.2,
                             "bottom"   = 1.2,
                             "none"     = .5))+
  scale_alpha_manual(values=c("top"      = 1,
                              "bottom"    = 1,
                              "none"      =.3))+
  labs(title="Linear Relationship Between Popularity and Extraversion for 100 Classes",
       subtitle="The 6 with the most extreme relationship have been highlighted red and blue")

interceptonlymodeltest <- brm(popular ~ 1 + (1 | class), 
                              data   = popular2data, 
                              warmup = 100, 
                              iter   = 200, 
                              chains = 2, 
                              inits  = "random",
                              cores  = 2)  #the cores function tells STAN to make use of 2 CPU cores simultaneously instead of just 1.

summary(interceptonlymodeltest)

interceptonlymodel <- brm(popular ~ 1 + (1|class),  
                          data = popular2data, 
                          warmup = 1000, iter = 3000, 
                          cores = 2, chains = 2, 
                          seed = 123) #to run the model

summary(interceptonlymodel)

hyp <- "sd_class__Intercept^2 / (sd_class__Intercept^2 + sigma^2) = 0"
hypothesis(interceptonlymodel, hyp, class = NULL)

model1 <- brm(popular ~ 1 + sex + extrav + (1|class),  
              data = popular2data, 
              warmup = 1000, iter = 3000, 
              cores = 2, chains = 2, 
              seed = 123) #to run the model

summary(model1)

model1tranformed <- ggs(model1) # the ggs function transforms the brms output into a longformat tibble, that we can use to make different types of plots.

ggplot(filter(model1tranformed, Parameter %in% c("b_Intercept", "b_extrav", "b_sex")),
       aes(x   = Iteration,
           y   = value, 
           col = as.factor(Chain)))+
  geom_line() +
  geom_vline(xintercept = 1000)+
  facet_grid(Parameter ~ . ,
             scale  = 'free_y',
             switch = 'y')+
  labs(title = "Caterpillar Plots", 
       col   = "Chains")

ggplot(filter(model1tranformed,
              Parameter == "b_Intercept", 
              Iteration > 1000),
       aes(x = value))+
  geom_density(fill  = "yellow", 
               alpha = .5)+
  geom_vline(xintercept = 0, 
             col  = "red",
             size = 1)+
  scale_x_continuous(name   = "Value",
                     limits = c(-1, 3)) + 
  geom_vline(xintercept = summary(model1)$fixed[1,3:4],
             col = "blue",
             linetype = 2) +
  theme_light() +
  labs(title = "Posterior Density of Intercept")

ggplot(filter(model1tranformed, Parameter == "b_extrav", Iteration > 1000), aes(x = value))+
  geom_density(fill = "orange", alpha = .5)+
  geom_vline(xintercept = 0, col = "red", size = 1)+
  scale_x_continuous(name = "Value", limits = c(-.2, .6))+ 
  geom_vline(xintercept = summary(model1)$fixed[3,3:4], col = "blue", linetype = 2)+
  theme_light()+
  labs(title = "Posterior Density of Regression Coefficient for Extraversion")

ggplot(filter(model1tranformed, Parameter == "b_sex", Iteration > 1000), aes(x = value))+
  geom_density(fill = "red", alpha = .5)+
  geom_vline(xintercept = 0, col = "red", size = 1)+
  scale_x_continuous(name = "Value", limits = c(-.2, 1.5))+ 
  geom_vline(xintercept = summary(model1)$fixed[2,3:4], col = "blue", linetype = 2)+
  theme_light()+
  labs(title = "Posterior Density of Regression Coefficient for Sex")

model2 <- brm(popular ~ 1 + sex + extrav + texp + (1|class),  
              data = popular2data, 
              warmup = 1000, iter = 3000, 
              cores = 2, chains = 2, 
              seed = 123)

summary(model2)

model3 <- brm(popular ~ 1 + sex + extrav + (1 + sex + extrav | class),  
              data = popular2data, 
              warmup = 1000, iter = 3000, 
              cores = 2, chains = 2, 
              seed = 123) #to run the model

summary(model3)

model4 <- brm(popular ~ 1 + sex + extrav + texp + (1 + extrav | class),  
              data = popular2data, 
              warmup = 1000, iter = 3000, 
              cores = 2, chains = 2, 
              seed = 123) #to run the model

summary(model4)

model5 <- brm(popular ~ 1 + sex + extrav + texp + extrav:texp + (1 + extrav|class), 
              data  = popular2data, warmup = 1000,
              iter  = 3000, chains = 2, 
              seed  = 123, control = list(adapt_delta = 0.97),
              cores = 2) # to reach a usuable number effective samples in the posterior 
              # distribution of the interaction effect, we need many more iteration. 
              # This sampler will take quite some time and you might want to run it with a 
              # few less iterations.

summary(model5)$fixed
summary(model5)$random

ggplot(data = popular2data, 
       aes(x   = extrav,
           y   = popular,
           col = as.factor(texp)))+
  viridis::scale_color_viridis(discrete = TRUE)+
  geom_point(size     = .7,
             alpha    = .8,
             position = "jitter")+
  geom_smooth(method = lm,
              se     = FALSE, 
              size   = 2,
              alpha  = .8)+
  theme_minimal()+
  labs(title    = "Linear Relationship for Different Years of Teacher Experience as Observed", 
       subtitle = "The linear relationship between the two is not the same for all classes", 
       col      = "Years of\nTeacher\nExperience")

simplemodel1 <- brm(popular ~ 1 + extrav + (1 + extrav | class), 
                    data = popular2data,
                    warmup = 1000, iter = 5000, chains = 2,  
                    seed = 123, control = list(adapt_delta = 0.96), 
                    save_all_pars = TRUE, cores= 2)


posteriorsimpelmodel1 <- as_tibble(t(posterior_samples(simplemodel1, pars = "extrav")[,-c(1:3)]))

teacherexperience <- popular2data %>%
  group_by(class) %>%
  summarise("teacherexperience" = mean(texp))

posteriorsimpelmodellong <- bind_cols(teacherexperience, posteriorsimpelmodel1) %>%
  gather(key = "key", value = "value", -teacherexperience, -class)%>%
  group_by(class) %>%
  mutate(meanperclass = mean(value))%>%
  ungroup()

ggplot()+
  ggridges::geom_density_ridges(data  = posteriorsimpelmodellong, 
                                aes(x      = value,
                                    y      = reorder(as.factor(class), meanperclass),
                                    height = ..density.., 
                                    fill   = as.factor(teacherexperience)),
                                scale = 3, 
                                alpha = .6) +
  scale_x_continuous(limits = c(-.5,.5))+
  geom_point(data = summarise(group_by(posteriorsimpelmodellong, class), mean = mean(meanperclass)),
             aes(x = mean, 
                 y = as.factor(class)),
             size = 1, 
             col  = "red")+
  viridis::scale_fill_viridis(discrete = TRUE)+
  geom_vline(xintercept = 0, 
             col        = "red")+
  labs(fill     = "Years of\nTeacher\nExperience",
       y        = "classes", 
       title    = expression(paste("Class Level Error of Regression Coefficient of Extraversion on Popularity (", u["2j"],")")),
       subtitle = expression(paste("posterior distribution of class level error of regression coefficient (", u["2j"],") per class with the means in red")), 
       caption  = expression(paste("Regression formula: popularity = ", gamma["00"], "+", gamma["20"],"*", extrav["ij"], "+", u["2j"], "*", extrav["ij"],"+", e["ij"] )))+
  annotate(geom     = "text", 
           x        = 0, 
           y        = 1.5, 
           label    = "Grand mean", 
           col      = "red", 
           family   = theme_get()$text[["family"]], 
           size     = theme_get()$text[["size"]]/3, 
           fontface = "italic")+
  theme_tufte()

distance95 <- posteriorsimpelmodellong %>%
  group_by(class) %>%
  summarise(lower95      = quantile(value, probs = .025),
            upper95      = quantile(value, probs = .975),
            distance     = upper95-lower95, 
            Meanestimate = mean(value)) %>%
  bind_cols(teacherexperience)%>%
  group_by(teacherexperience)%>%
  summarise(mean         = mean(distance), 
            Meanestimate = mean(Meanestimate),
            lower        = mean(lower95), 
            upper        = mean(upper95),
            meanCCI      = paste("[",sprintf("%.4f",round(lower,4)), ":", sprintf("%.4f",round(upper,4)), "]")) 

distance95 <- mutate(distance95, Quadratic = teacherexperience^2)

model <- lm(mean ~ teacherexperience + Quadratic, data = distance95)

summary(model)

dat <- data.frame(teacherexperience = c(2:25),
                  Quadratic         = c(2:25)^2)

dat$yhat <- predict(model, dat)

ggplot()+
  geom_line(data  = dat, 
            aes(x = teacherexperience,
                y = yhat), 
            linetype = "dotted",
            size     = 1)+
  geom_point(data = distance95, 
             aes(x   = teacherexperience, 
                 y   = mean, 
                 col = Meanestimate))+
  geom_text(data  = distance95, 
            aes(x     = teacherexperience, 
                y     = mean, 
                label = meanCCI,
                col   = Meanestimate),
            hjust    = .5, 
            vjust    = -.15, 
            family   = theme_get()$text[["family"]], 
            size     = 3, 
            fontface = "italic")+
  annotate(geom = "text", 
           x    = 15, 
           y    = 0.44,
           label    = "CCI distance=0.45-0.013*texp+0.00036*texp^2\nR^2=0.79", 
           colour   = "black",  
           family   = theme_get()$text[["family"]], 
           size     = theme_get()$text[["size"]]/2, 
           fontface = "italic")+
  viridis :: scale_color_viridis(discrete = F, direction = -1)+
  labs(y        = "95% CCI distance",
       title    = expression(paste("Mean CCI Distance of Posterior of (", u["2j"], ") for Different Years of Texp")),
       subtitle = "In brackets the actual CCIs and in colour the parameter estimate",
       col      = expression(paste("estimate ", u["2j"])))+
  scale_x_continuous(breaks = 2:25)+
  theme_tufte()


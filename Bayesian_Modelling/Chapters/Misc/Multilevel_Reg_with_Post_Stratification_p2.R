# https://cran.r-project.org/web/packages/rstanarm/vignettes/mrp.html
# https://github.com/stan-dev/rstanarm/blob/master/vignettes/mrp.Rmd


# not evaluated to avoid tidyverse dependency 
income_popn <- poststrat %>%
  group_by(income) %>%
  summarize(Num=sum(N)) %>%
  mutate(PROP=Num/sum(Num),TYPE='Popn',VAR='Income',CAT=income) %>%
  ungroup()
income_data <- sample %>%
  group_by(income) %>%
  summarise(Num=n()) %>%
  mutate(PROP=Num/sum(Num),TYPE='Sample',VAR='Income',CAT=income) %>%
  ungroup()
income<-rbind(income_data[,2:6],income_popn[,2:6])
age_popn <- poststrat%>%
  group_by(age)%>%
  summarize(Num=sum(N))%>%
  mutate(PROP=Num/sum(Num),TYPE='Popn',VAR='Age',CAT=age)%>%
  ungroup()
age_data <- sample%>%
  group_by(age)%>%
  summarise(Num=n())%>%
  mutate(PROP=Num/sum(Num),TYPE='Sample',VAR='Age',CAT=age)%>%
  ungroup()
age <- rbind(age_data[,2:6],age_popn[,2:6] )
eth_popn <- poststrat%>%
  group_by(eth)%>%
  summarize(Num=sum(N))%>%
  mutate(PROP=Num/sum(Num),TYPE='Popn',VAR='Ethnicity',CAT=eth)%>%
  ungroup()
eth_data <- sample%>%
  group_by(eth)%>%
  summarise(Num=n())%>%
  mutate(PROP=Num/sum(Num),TYPE='Sample',VAR='Ethnicity',CAT=eth)%>%
  ungroup()
eth<-rbind(eth_data[,2:6],eth_popn[,2:6])
male_popn <- poststrat%>%
  group_by(male)%>%
  summarize(Num=sum(N))%>%
  mutate(PROP=Num/sum(Num),TYPE='Popn',VAR='Male',CAT=male)%>%
  ungroup()
male_data <- sample%>%
  group_by(male)%>%
  summarise(Num=n())%>%
  mutate(PROP=Num/sum(Num),TYPE='Sample',VAR='Male',CAT=male)%>%
  ungroup()
male <- rbind(male_data[,2:6],male_popn[,2:6])
state_popn <- poststrat%>%
  group_by(state)%>%
  summarize(Num=sum(N))%>%
  mutate(PROP=Num/sum(poststrat$N),TYPE='Popn',VAR='State',CAT=state)%>%
  ungroup()
state_plot_data <- sample%>%
  group_by(state)%>%
  summarise(Num=n())%>%
  mutate(PROP=Num/nrow(sample),TYPE='Sample',VAR='State',CAT=state)%>%
  ungroup()
state_plot_data <- rbind(state_plot_data[,2:6],state_popn[,2:6])
state_plot_data$TYPE <- factor(state_plot_data$TYPE, levels = c("Sample","Popn"))
plot_data <- rbind(male,eth,age,income)
plot_data$TYPE <- factor(plot_data$TYPE, levels = c("Sample","Popn"))
save(state_plot_data, file = "state_plot_data.rda", version = 2)
save(plot_data, file = "plot_data.rda", version = 2)
getwd()

load("plot_data.rda") # created in previous chunk
ggplot(data=plot_data, aes(x=as.factor(CAT), y=PROP, group=as.factor(TYPE), linetype=as.factor(TYPE))) +
  geom_point(stat="identity",colour='black')+
  geom_line()+
  facet_wrap( ~ VAR, scales = "free",nrow=1,ncol=5)+
  theme_bw()+
  scale_fill_manual(values=c('#1f78b4','#33a02c',
                             '#e31a1c','#ff7f00','#8856a7'),guide=FALSE)+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1), labels=c('0%','25%',"50%","75%","100%"))+
  scale_alpha_manual(values=c(1, .3))+
  ylab('Proportion')+
  labs(alpha='')+
  theme(legend.position="bottom",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        axis.text=element_text(size=10),
        strip.text=element_text(size=10),
        strip.background = element_rect(fill='grey92'))

load("state_plot_data.rda") # created in previous chunk
ggplot(data=state_plot_data, aes(x=as.factor(CAT), y=PROP, group=as.factor(TYPE),    linetype=as.factor(TYPE))) +
  geom_point(stat="identity",colour='black')+
  geom_line()+
  facet_wrap( ~ VAR)+
  theme_bw()+
  scale_fill_manual(values=c('#1f78b4','#33a02c',
                             '#e31a1c','#ff7f00','#8856a7'),guide=FALSE)+
  scale_y_continuous(breaks=c(0,.025,.05,1), labels=c('0%','2.5%',"5%","100%"),expand=c(0,0),limits=c(0,.06))+
  scale_alpha_manual(values=c(1, .3))+
  ylab('Proportion')+
  labs(alpha='')+
  theme(legend.position="bottom",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=8,angle=90),
        strip.text=element_text(size=10),
        strip.background = element_rect(fill='grey92'))

# not evaluated to avoid dependency on tidyverse
#Summarise
summary_by_poststrat_var <- sample %>%
  gather(variable,category,c("income","eth","age","male")) %>%
  group_by(variable,category) %>%
  #Wald confidence interval
  summarise(y_mean=mean(cat_pref),y_sd=sqrt(mean(cat_pref)*(1-mean(cat_pref))/n())) %>%
  ungroup()
summary_by_poststrat_var$variable <- as.factor(summary_by_poststrat_var$variable)
levels(summary_by_poststrat_var$variable) <- list('Age'='age','Ethnicity'='eth','Income'='income','Male'='male')
save(summary_by_poststrat_var, file = "summary_by_poststrat_var.rda", 
     version = 2)

load("summary_by_poststrat_var.rda") # created in previous chunk
ggplot(data=summary_by_poststrat_var, aes(x=as.factor(category), y=y_mean,group=1)) +
  geom_errorbar(aes(ymin=y_mean-y_sd, ymax=y_mean+y_sd), width=0)+
  geom_line()+
  geom_point()+
  scale_colour_manual(values=c('#1f78b4','#33a02c','#e31a1c','#ff7f00',
                               '#8856a7'))+theme_bw()+
  facet_wrap(~variable,scales = "free_x",nrow=1,ncol=5)+
  scale_y_continuous(breaks=c(.5,.75,1), labels=c("50%","75%",
                                                  "100%"), limits=c(0.4-.4*.05,.9),expand = c(0,0))+
  labs(x="",y="Cat preference")+
  theme(legend.position="none",
        axis.title.y=element_text(size=10),
        axis.title.x=element_blank(),
        axis.text=element_text(size=10),
        strip.text=element_text(size=10),
        strip.background = element_rect(fill='grey92'))


#Summarise
interaction <- sample %>%
  gather(variable, category, c("age", "eth")) %>%
  group_by(variable, category, male) %>%
  summarise(y_mean = mean(cat_pref), 
            y_sd = sqrt(mean(cat_pref) * (1 - mean(cat_pref)) / n())) %>%
  ungroup()
#Tidy for nice facet labels
interaction$variable <- as.factor(interaction$variable)
levels(interaction$variable) <- list('Ethnicity' = 'eth', 'Age' = 'age')
save(interaction, file = "interaction.rda", version = 2)

load("interaction.rda") # created in previous chunk
ggplot(data=interaction, aes(x=as.factor(category), y=y_mean, colour=as.factor(male),group=as.factor(male))) +
  geom_errorbar(aes(ymin=y_mean-y_sd, ymax=y_mean+y_sd),width=0 )+
  geom_line(aes(x=as.factor(category), y=y_mean,colour=as.factor(male)))+
  geom_point()+
  facet_wrap(~variable,scales = "free_x",nrow=1,ncol=2)+
  labs(x="",y="Cat preference",colour='Gender')+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1), labels=c("0%",'25%',"50%","75%",
                                                        "100%"), limits=c(0,1),expand=c(0,0))+
  scale_colour_manual(values=c('#4575b4','#d73027'))+theme_bw()+
  theme(axis.title=element_text(size=10),
        axis.text=element_text(size=10),
        legend.position='none',
        strip.text=element_text(size=10),
        strip.background = element_rect(fill='grey92'))

preference_by_state <- sample %>%
  group_by(state) %>%
  summarise(y_mean = mean(cat_pref), 
            y_sd = sqrt(mean(cat_pref) * (1 - mean(cat_pref)) / n())) %>%
  ungroup()
save(preference_by_state, file = "preference_by_state.rda", version = 2)



load("preference_by_state.rda")
compare <- ggplot(data=preference_by_state, aes(x=state, y=y_mean,group=1)) +
  geom_ribbon(aes(ymin=y_mean-y_sd,ymax=y_mean+y_sd,x=state),fill='lightgrey',alpha=.7)+
  geom_line(aes(x=state, y=y_mean))+
  geom_point()+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1), 
                     labels=c("0%","25%","50%","75%","100%"), 
                     limits=c(0,1), expand=c(0,0))+
  scale_x_discrete(drop=FALSE)+
  scale_colour_manual(values=c('#1f78b4','#33a02c','#e31a1c','#ff7f00',
                               '#8856a7'))+
  theme_bw()+
  labs(x="States",y="Cat preference")+
  theme(legend.position="none",
        axis.title=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=90,size=8),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10))
compare2 <- ggplot()+
  geom_hline(yintercept = mean(sample$cat_pref),size=.8)+
  geom_text(aes(x = 5.2, y = mean(sample$cat_pref)+.025, label = "Sample"))+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1), 
                     labels=c("0%","25%","50%","75%","100%"),
                     limits=c(-0.25,1.25),expand=c(0,0))+
  theme_bw()+
  labs(x="Popn",y="")+
  theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=10),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10))
bayesplot_grid(compare,compare2, 
               grid_args = list(nrow=1, widths = c(8,1)))


fit <- stan_glmer(
  cat_pref ~ factor(male) + factor(male) * factor(age) + 
    (1 | state) + (1 | age) + (1 | eth) + (1 | income),
  family = binomial(link = "logit"),
  data = sample
)

posterior_prob <- posterior_linpred(fit, transform = TRUE, newdata = poststrat)
poststrat_prob <- posterior_prob %*% poststrat$N / sum(poststrat$N)
model_popn_pref <- c(mean = mean(poststrat_prob), sd = sd(poststrat_prob))
round(model_popn_pref, 3)

sample_popn_pref <- mean(sample$cat_pref)
round(sample_popn_pref, 3)

compare2 <- compare2 +
  geom_hline(yintercept = model_popn_pref[1], colour = '#2ca25f', size = 1) +
  geom_text(aes(x = 5.2, y = model_popn_pref[1] + .025), label = "MRP", colour = '#2ca25f')
bayesplot_grid(compare, compare2, 
               grid_args = list(nrow = 1, widths = c(8, 1)))

true_popn_pref <- sum(true_popn$cat_pref * poststrat$N) / sum(poststrat$N)
round(true_popn_pref, 3)

compare2 <- compare2 +
  geom_hline(yintercept = mean(true_popn_pref), linetype = 'dashed', size = .8) +
  geom_text(aes(x = 5.2, y = mean(true_popn_pref) - .025), label = "True")
bayesplot_grid(compare, compare2, 
               grid_args = list(nrow = 1, widths = c(8, 1)))

tate_df <- data.frame(
  State = 1:50,
  model_state_sd = rep(-1, 50),
  model_state_pref = rep(-1, 50),
  sample_state_pref = rep(-1, 50),
  true_state_pref = rep(-1, 50),
  N = rep(-1, 50)
)
for(i in 1:length(levels(as.factor(poststrat$state)))) {
  poststrat_state <- poststrat[poststrat$state == i, ]
  posterior_prob_state <- posterior_linpred(
    fit,
    transform = TRUE,
    draws = 1000,
    newdata = as.data.frame(poststrat_state)
  )
  poststrat_prob_state <- (posterior_prob_state %*% poststrat_state$N) / sum(poststrat_state$N)
  #This is the estimate for popn in state:
  state_df$model_state_pref[i] <- round(mean(poststrat_prob_state), 4)
  state_df$model_state_sd[i] <- round(sd(poststrat_prob_state), 4)
  #This is the estimate for sample
  state_df$sample_state_pref[i] <- round(mean(sample$cat_pref[sample$state == i]), 4)
  #And what is the actual popn?
  state_df$true_state_pref[i] <-
    round(sum(true_popn$cat_pref[true_popn$state == i] * poststrat_state$N) /
            sum(poststrat_state$N), digits = 4)
  state_df$N[i] <- length(sample$cat_pref[sample$state == i])
}
state_df[c(1,3:6)]
state_df$State <- factor(state_df$State, levels = levels(sample$state))

round(100 * c(
  mean = mean(abs(state_df$sample_state_pref-state_df$true_state_pref), na.rm = TRUE),
  max = max(abs(state_df$sample_state_pref-state_df$true_state_pref), na.rm = TRUE)
))

round(100 * c(
  mean = mean(abs(state_df$model_state_pref-state_df$true_state_pref)),
  max = max(abs(state_df$model_state_pref-state_df$true_state_pref))
))

compare <- compare +
  geom_point(data=state_df, mapping=aes(x=State, y=model_state_pref),
             inherit.aes=TRUE,colour='#238b45')+
  geom_line(data=state_df, mapping=aes(x=State, y=model_state_pref,group=1),
            inherit.aes=TRUE,colour='#238b45')+
  geom_ribbon(data=state_df,mapping=aes(x=State,ymin=model_state_pref-model_state_sd,
                                        ymax=model_state_pref+model_state_sd,group=1), 
              inherit.aes=FALSE,fill='#2ca25f',alpha=.3)+
  geom_point(data=state_df, mapping=aes(x=State, y=true_state_pref),
             alpha=.5,inherit.aes=TRUE)+
  geom_line(data=state_df, mapping=aes(x=State, y=true_state_pref),
            inherit.aes = TRUE,linetype='dashed')
bayesplot_grid(compare, compare2, 
               grid_args = list(nrow = 1, widths = c(8, 1)))

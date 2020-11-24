library('gapminder')
library('tidyr')
library('rstan')
library('ggplot2')
library('rethinking')
library('dplyr')
library('purrr')
library('broom')

gapminder %>% 
  ggplot(aes(year, lifeExp, group = country)) +
  geom_line(alpha = 5/9)

gapminder <- gapminder %>%
        mutate( year1950 = year - 1950 )

by_country <- gapminder %>% 
  group_by(country, continent) %>% 
  nest()

country_model <- function(df) {
  lm(lifeExp ~ year, data = df)
}

models <- by_country %>% 
  mutate(
    model = data %>% map(country_model))


by_country <- by_country %>% 
  mutate(model = map(data, country_model))


library(brms)

fit <-
  brm(data = gapminder,
      family = gaussian,
      lifeExp ~ 0 + Intercept + year1950,
      prior = c(prior(normal(0, 2), class = b),
                prior(student_t(3, 1, 2), class = sigma)),
      seed = 1)

country_model_bayes <- function(df) {
  update( fit, newdata = df )
}

bayes_mod <- by_country %>% 
  mutate(
    model = data %>% map(country_model_bayes))

by_country <- by_country %>% 
  mutate(bayes_mod = map(data, country_model_bayes))

by_country <- gapminder %>% 
  group_by(country, continent) %>% 
  nest %>%
  mutate(data2=map(data, ~select(., c(lifeExp, year1950))))

mll <- ulam(
  alist(
    lifeExp ~ dnorm(mu, sigma),
    mu <- a + b*year1950,
    a ~ dnorm(60, 10),
    b ~ dnorm(30, 10),
    sigma ~ dunif(0, 10)
  ), data=list(lifeExp=gapminder$lifeExp, year1950=gapminder$year1950),
               chains=4, cores=4
)

precis(mll, depth=2)

country_model_stan <- function(df) {
  stan(fit=mll@stanfit, data=df, chains=4, cores=4)
}


models <- by_country %>%
  mutate(stan_mod = map(data2, country_model_stan))

by_country$model[[1]]
by_country$bayes_mod[[1]]

lm_models <- models %>%
  mutate(
    tidy = map(model, broom::tidy),
    glance =  map(model, broom::glance),
    augment = map(model, broom::augment)
  )

bayes_models <- bayes_mod %>%
  mutate(
    tidy = map(bayes_mod, broom::tidy),
    glance =  map(bayes_mod, broom::glance),
    augment = map(bayes_mod, broom::augment)
  )



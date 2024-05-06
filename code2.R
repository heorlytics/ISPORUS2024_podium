library(tidyverse)
library(flexsurvcure)
library(pammtools)
library(dplyr)
library(broom)
library(survival)
library(flexsurv)
library(flexsurvcure)


setwd("C:\\Users\\ShubhramPandey\\OneDrive - PHARMACOEVIDENCE PRIVATE LIMITED\\Documents\\ISPOR US 24\\Podium")
data = read.csv("data.csv") %>% select(event,survtime,age,sex)
## Create age at diagnosis in days - used later for matching to expected rates
data$agedays <- floor(data$age * 365.25)
## Survival time in years
data$survyrs <- data$survtime / 365.25
data$sex <- as.factor(data$sex)

## Obtain attained age and attained calendar year in (whole) years
data <- data %>% mutate(attained.age.yr = floor(age + survtime/365.25))


# 2. lifetables ----
# We will use the US lifetables that come with the survival package
# First, let's reshape US lifetable to be a tidy data.frame and convert rates to
# per person-year as our survival analysis time scale will be in years
survexp.us.df <- as.data.frame.table(survexp.us[,,"2014"], responseName = "exprate") %>%
  mutate(exprate = 365.25 * exprate)
survexp.us.df$age <- as.numeric(as.character(survexp.us.df$age))
# Now we merge in (left join) the US rates at the event times in 
data <- data %>% left_join(survexp.us.df, by = c("attained.age.yr"="age", 
                                                    "sex"="sex")) 
surv = Surv(survyrs, event)~1

km <- survfit(surv, data = data)
p = summary(km, times = seq(0,2,0.01))

kap = broom::tidy(km)

(ggplot(kap, aes(time, estimate, ymin = conf.low, ymax = conf.high)) +  
    geom_stepribbon(fill = "black", alpha = .2) +  
    geom_step(color = "black")+
    annotate("text", x=1.5, y=1, label= "Median: 0.793 years") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format())+
    labs(x = "Time (years)", y = "Survival probability (in %)") +
    theme_minimal())

hazard <- tidy(muhaz::muhaz(data$survyrs, data$event))
(ggplot(hazard, aes(time, estimate)) +  
    geom_line(color = "black")+
    labs(x = "Time (years)", y = "Hazard") +
    theme_minimal())

# 3. Model fitting to ----
# Fit a standard parametric, excess hazard and mixture-cure model to the Good 
# prognosis group only.
models <- list()
# standard parametric
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  models[[dist]] <- flexsurvreg(surv, data=data, dist=dist)
}
# excess hazard models
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  method <- ifelse(dist=="exp", "BFGS", "Nelder-Mead")
  models[[paste0(dist,".excesshazard")]] <- 
    flexsurvreg(surv, 
                data=data, 
                dist=dist,
                bhazard=exprate, 
                method=method)
}

# excess hazard mixture-cure models
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  models[[paste0(dist,".excesshazardcure")]] <- 
    flexsurvcure(surv, 
                 data=data, 
                 dist=dist,
                 bhazard=exprate,
                 method="Nelder-Mead")
}

models <- tibble(Distribution = c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
                                  "gengamma")) %>%
  mutate(model = map(Distribution,
                     ~ flexsurvcure(surv, data = data, dist = .x, bhazard = exprate)),
         glanced = map(model, glance),
         tidied = map(model, tidy))

df3 = models %>%
  unnest(tidied) %>%
  filter(term == "theta") %>%
  transmute(Distribution, `Cure fraction` = estimate)

# 4. AIC statisics ----
# For the excess hazard models the reported AIC is from a partial likelihood
# Therefore the AIC from an excess hazard model cannot be directly compared to
# an AIC from a standard parametric model
aic.tbl <- tibble()
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  aic.tbl <- aic.tbl %>% bind_rows(tibble(dist=dist, Method="param", 
                                          aic=models[[dist]]$AIC))
  # Excess hazard models
  dist.eh <- paste0(dist,".excesshazard")
  aic <- models[[dist.eh]]$AIC
  aic.tbl <- aic.tbl %>% bind_rows(tibble(dist=dist, Method="eh", aic=aic))
  
  # Excess hazard cure models
  dist.cure <- paste0(dist,".excesshazardcure")
  aic <- models[[dist.cure]]$AIC
  aic.tbl <- aic.tbl %>% bind_rows(tibble(dist=dist, Method="cure", aic=aic))
}
aic.tbl.wide <- aic.tbl %>% pivot_wider(names_from = Method, values_from = aic) 
aic.tbl.wide %>% as.data.frame()

bic.tbl <- tibble()
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  bic.tbl <- bic.tbl %>% bind_rows(tibble(dist=dist, Method="param", 
                                          aic=BIC(models[[dist]])))
  # Excess hazard models
  dist.eh <- paste0(dist,".excesshazard")
  aic <- BIC(models[[dist.eh]])
  bic.tbl <- bic.tbl %>% bind_rows(tibble(dist=dist, Method="eh", aic=aic))
  
  # Excess hazard cure models
  dist.cure <- paste0(dist,".excesshazardcure")
  aic <- models[[dist.cure]] %>% BIC()
  bic.tbl <- bic.tbl %>% bind_rows(tibble(dist=dist, Method="cure", aic=aic))
}
bic.tbl.wide <- bic.tbl %>% pivot_wider(names_from = Method, values_from = aic) 
bic.tbl.wide %>% as.data.frame()


# 5. Predicted and extrapolated all-cause survival and hazard ----
ss.surv <- tibble()
# We use the standsurv command in the flexsurv package to make the predictions
# Predicting all-cause survival over next 30-years in increments of 0.2 years
# Standard parametric model
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  ss.surv.new <- standsurv(models[[dist]],
                           type="survival",
                           t=seq(0,20, by=0.01),
                           ci=F
  ) %>% 
    bind_cols(dist = dist) %>%
    bind_cols("Relative Survival" = FALSE) %>% 
    bind_cols(Cure = FALSE)
  ss.surv <- ss.surv %>% 
    bind_rows(ss.surv.new)
}

# Excess hazard model
# To get predictions of marginal all-cause survival, standsurv multiplies the 
# predicted relative survival function with the expected survival function for 
# each individual and then averages.
# We therefore must supply the expected ratetable for these calculations
# We also must scale from the ratetable time scale (days) to the regression 
# model time scale (years) using the scale.ratetable argument.
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  ss.surv.new <- standsurv(models[[paste0(dist,".excesshazard")]],
                           type="survival",
                           t=seq(0,20, by=0.01),
                           ci=F,
                           rmap=list(sex = sex,
                                     age = agedays
                           ),
                           ratetable = survexp.us[,,"2014"],
                           scale.ratetable = 365.25,
                           newdata = data
  ) %>% 
    bind_cols(dist = dist) %>%
    bind_cols("Relative Survival" = TRUE) %>% 
    bind_cols(Cure = FALSE)
  ss.surv <- ss.surv %>% 
    bind_rows(ss.surv.new)
}

## Excess hazard cure model
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  ss.surv.new <- standsurv(models[[paste0(dist,".excesshazardcure")]],
                           type="survival",
                           t=seq(0,20, by=0.01),
                           ci=F,
                           rmap=list(sex = sex,
                                     age = agedays
                           ),
                           ratetable = survexp.us[,,"2014"],
                           scale.ratetable = 365.25,
                           newdata = data
  ) %>% 
    bind_cols(dist = dist) %>%
    bind_cols("Relative Survival" = TRUE) %>% 
    bind_cols(Cure = TRUE)
  ss.surv <- ss.surv %>% 
    bind_rows(ss.surv.new)
  ss.surv
}
# 6. Plot of extrapolated survival (faceted)----
####################################################################

km1 = data.frame(time = p$time, 
                 at1 = p$surv, 
                 dist = "km",
                 `Relative Survival` = FALSE, 
                 Cure = FALSE, 
                 Method = "a) Standard parametric models",
                 check.names = F
                 ) %>% as.tibble()
km2 = data.frame(time = p$time, 
                 at1 = p$surv, 
                 dist = "km",
                 `Relative Survival` = TRUE, 
                 Cure = FALSE, 
                 Method = "b) Excess hazard (no cure) models",
                 check.names = F
) %>% as.tibble()

km3 = data.frame(time = p$time, 
                 at1 = p$surv, 
                 dist = "km",
                 `Relative Survival` = TRUE, 
                 Cure = TRUE, 
                 Method = "c) Excess hazard (cure) models",
                 check.names = F
) %>% as.tibble()

ss.surv  = ss.surv %>% 
          bind_rows(km1,km2,km3) %>% 
          mutate(dist=factor(dist)) %>% 
          mutate(dist=fct_relevel(dist,c("km","exp","weibull","gamma","gompertz","llogis","lnorm","gengamma")))


colnames(ss.surv)[3] = "Distribution"

ss.surv <- ss.surv %>% 
  mutate(Method = ifelse(`Relative Survival`==FALSE & Cure==FALSE,
                         "a) Standard parametric models",
                         ifelse(Cure==FALSE, 
                                "b) Excess hazard (no cure) models", "c) Excess hazard (cure) models")))
ss.surv$Method <- factor(ss.surv$Method, 
                         levels = c("a) Standard parametric models", 
                                    "b) Excess hazard (no cure) models", 
                                    "c) Excess hazard (cure) models"))

col  = RColorBrewer::brewer.pal(7,"Set3")
plot0 = ggplot(ss.surv) + geom_line(aes(x=time,y=at1,color=Distribution), lwd=0.5) + 
          facet_wrap(~ Method, ncol=1) + 
          theme_bw() +
          ylab("Survival") +
          xlab("Time (years)")+
          scale_color_manual(labels = c("Kaplan-Meier","Exponential","Weibull","Gamma","Gompertz","Log-logistic","Log-normal","Generalized Gamma"),
            values=c("black", col[1:7]))+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent'),
    strip.background = element_rect(
      color="black", fill="transparent", size=0, linetype="solid"
    ),
    strip.text = element_text(hjust = 0, size=10),
    legend.position="right",
    legend.title=element_blank()
  )
ggsave('myplot3.png', plot0, bg='transparent',width = 32,height = 15,units = "cm")


# Hazard plot -------------------------------------------------------------

# 5. Predicted and extrapolated all-cause survival and hazard ----
ss.surv <- tibble()
# We use the standsurv command in the flexsurv package to make the predictions
# Predicting all-cause survival over next 30-years in increments of 0.2 years
# Standard parametric model
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  ss.surv.new <- standsurv(models[[dist]],
                           type="hazard",
                           t=seq(0,20, by=0.01),
                           ci=F
  ) %>% 
    bind_cols(dist = dist) %>%
    bind_cols("Relative Survival" = FALSE) %>% 
    bind_cols(Cure = FALSE)
  ss.surv <- ss.surv %>% 
    bind_rows(ss.surv.new)
}

# Excess hazard model
# To get predictions of marginal all-cause survival, standsurv multiplies the 
# predicted relative survival function with the expected survival function for 
# each individual and then averages.
# We therefore must supply the expected ratetable for these calculations
# We also must scale from the ratetable time scale (days) to the regression 
# model time scale (years) using the scale.ratetable argument.
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  ss.surv.new <- standsurv(models[[paste0(dist,".excesshazard")]],
                           type="hazard",
                           t=seq(0,20, by=0.01),
                           ci=F,
                           rmap=list(sex = sex,
                                     age = agedays
                           ),
                           ratetable = survexp.us[,,"2014"],
                           scale.ratetable = 365.25,
                           newdata = data
  ) %>% 
    bind_cols(dist = dist) %>%
    bind_cols("Relative Survival" = TRUE) %>% 
    bind_cols(Cure = FALSE)
  ss.surv <- ss.surv %>% 
    bind_rows(ss.surv.new)
}

## Excess hazard cure model
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  ss.surv.new <- standsurv(models[[paste0(dist,".excesshazardcure")]],
                           type="hazard",
                           t=seq(0,20, by=0.01),
                           ci=F,
                           rmap=list(sex = sex,
                                     age = agedays
                           ),
                           ratetable = survexp.us[,,"2014"],
                           scale.ratetable = 365.25,
                           newdata = data
  ) %>% 
    bind_cols(dist = dist) %>%
    bind_cols("Relative Survival" = TRUE) %>% 
    bind_cols(Cure = TRUE)
  ss.surv <- ss.surv %>% 
    bind_rows(ss.surv.new)
  ss.surv
}

km1 = data.frame(time = hazard$time, 
                 at1 = hazard$estimate, 
                 dist = "km",
                 `Relative Survival` = FALSE, 
                 Cure = FALSE, 
                 Method = "a) Standard parametric models",
                 check.names = F
) %>% as.tibble()
km2 = data.frame(time = hazard$time, 
                 at1 = hazard$estimate,  
                 dist = "km",
                 `Relative Survival` = TRUE, 
                 Cure = FALSE, 
                 Method = "b) Excess hazard (no cure) models",
                 check.names = F
) %>% as.tibble()

km3 = data.frame(time = hazard$time, 
                 at1 = hazard$estimate, 
                 dist = "km",
                 `Relative Survival` = TRUE, 
                 Cure = TRUE, 
                 Method = "c) Excess hazard (cure) models",
                 check.names = F
) %>% as.tibble()

ss.surv  = ss.surv %>% 
  bind_rows(km1,km2,km3) %>% 
  mutate(dist=factor(dist)) %>% 
  mutate(dist=fct_relevel(dist,c("km","exp","weibull","gamma","gompertz","llogis","lnorm","gengamma")))


colnames(ss.surv)[3] = "Distribution"

ss.surv <- ss.surv %>% 
  mutate(Method = ifelse(`Relative Survival`==FALSE & Cure==FALSE,
                         "a) Standard parametric models",
                         ifelse(Cure==FALSE, 
                                "b) Excess hazard (no cure) models", "c) Excess hazard (cure) models")))
ss.surv$Method <- factor(ss.surv$Method, 
                         levels = c("a) Standard parametric models", 
                                    "b) Excess hazard (no cure) models", 
                                    "c) Excess hazard (cure) models"))

col  = RColorBrewer::brewer.pal(7,"Set3")
plot1 = ggplot(ss.surv) + geom_line(aes(x=time,y=at1,color=Distribution), lwd=0.5) + 
  facet_wrap(~ Method, ncol=1) + 
  theme_bw() +
  ylab("Hazard") +
  xlab("Time (years)")+
  scale_color_manual(labels = c("Kaplan-Meier","Exponential","Weibull","Gamma","Gompertz","Log-logistic","Log-normal","Generalized Gamma"),
                     values=c("black", col[1:7]))+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent'),
    strip.background = element_rect(
      color="black", fill="transparent", size=0, linetype="solid"
    ),
    strip.text = element_text(hjust = 0, size=10),
    legend.position="right",
    legend.title=element_blank()
  )
ggsave('myplot4.png', plot1, bg='transparent',width = 32,height = 16,units = "cm")

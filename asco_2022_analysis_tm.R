# This is a loose first analysis of the 2022 asco-bacteriacide
# osmia cornifrosn experiments
# run in the spring of 2022


# first a few packages
library(tidyverse) # for everything
library(rstatix) # for pipe friendly stats
library(lme4) # for models

library(survival)
library(survminer)
library(broom)

rm(list=ls()) # then clean the environment


# then read in the data
summary_df <- read_csv("input/asco_plate_maps_summary.csv")

summary_df$end_status <- as.factor(summary_df$end_status)
summary_df$plate <- as.factor(summary_df$plate)
summary_df$well <- as.factor(summary_df$well)
summary_df$larva_nest <- as.factor(summary_df$larva_nest)
summary_df$larva_number <- as.factor(summary_df$larva_number)
summary_df$treatment <- as.factor(summary_df$treatment)
summary_df$larva_stage <- as.factor(summary_df$larva_stage)


df<- summary_df%>%
  filter(end_status != "removed_dead")

# then a table perchance?
# a little graffy waffy?


df %>%
  ggplot()+
  geom_bar(aes(x=end_status, fill=treatment))

# and then just looking at the dead ones
df %>%
  filter(status=="dead")%>%
  ggplot()+
  geom_bar(aes(x=end_status, fill=treatment))


def_to_cocoon<- df %>%
  ggplot()+
  geom_boxplot(aes(x = treatment, y=def_to_cocoon, fill=treatment))+
  geom_point(aes(x = treatment, y=def_to_cocoon, fill=treatment))

ggsave(plot= def_to_cocoon, "output/def_to_cocoon.jpeg")

df %>%
  kruskal_test(def_to_cocoon~treatment)

df %>%
  dunn_test(def_to_cocoon~treatment)





def_to_cocoon_model <- lm(data=df, def_to_cocoon~treatment)
summary(def_to_cocoon_model)                   
                             
write.csv(tidy(def_to_cocoon_model), "output/def_to_cocoon_model.csv")



death_model <- glmer(data=df, status~treatment + (1|plate), family=binomial,
                     glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



summary(glm(data=df, end_status~treatment, family=binomial))


# test of a cox model


a.cox <- coxph(Surv(time_to_event, status)~treatment, data=df)


new_df <- with(df,
               data.frame(treatment = c("LbA","Lb","Hb","HbA","A","C") )
)
new_df


ggsurvplot(surv_fit(a.cox, data = df), pallette="#2E9FDF", ggtheme = theme_minimal())



fit <- survfit(a.cox, newdata = new_df)
ggsurvplot(fit, conf.int = TRUE, palette = "Dark2", 
           censor = FALSE, surv.median.line = "hv")









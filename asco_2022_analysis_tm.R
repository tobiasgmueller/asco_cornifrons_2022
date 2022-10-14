# This is a loose first analysis of the 2022 asco-bacteriacide
# osmia cornifrosn experiments
# run in the spring of 2022


# packages ####
library(tidyverse) # for everything
library(rstatix) # for pipe friendly stats
library(lme4) # for models

library(survival)
library(survminer)
library(broom)

rm(list=ls()) # then clean the environment


# read in the data ####
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
  filter(end_status=="dead")%>%
  ggplot()+
  geom_bar(aes(x=end_status, fill=treatment))+
  ggtitle("break down of just dead")


def_to_cocoon<- df %>%
  ggplot(aes(x = treatment, y=def_to_cocoon, fill=treatment))+
  geom_boxplot()+
  geom_point(position = position_jitter(w = 0.1, h = 0))

ggsave(plot= def_to_cocoon, "output/def_to_cocoon.jpeg")

df %>%
  kruskal_test(def_to_cocoon~treatment)

df %>%
  dunn_test(def_to_cocoon~treatment)





def_to_cocoon_model <- lm(data=df, def_to_cocoon~treatment)
summary(def_to_cocoon_model)                   
                             
write.csv(tidy(def_to_cocoon_model), "output/def_to_cocoon_model.csv")



death_model <- glmer(data=df, end_status~treatment + (1|plate), family=binomial,
                     glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



summary(glm(data=df, end_status~treatment, family=binomial))

# remaining pollen ####

df %>%
  drop_na(remaining_pollen)%>%
  ggplot(aes(x=remaining_pollen, fill=treatment))+
  geom_bar(position="fill")


remaining_pollen <- df %>%
  drop_na(remaining_pollen)%>%
  ggplot(aes(x=treatment, y=remaining_pollen, fill=treatment))+
  geom_boxplot(alpha=.6)+
  geom_point(position = position_jitter(w = 0.1, h = 0))
ggsave(plot= remaining_pollen, "output/remaining_pollen.jpeg", width = 10, height = 12)

df %>%
  drop_na(remaining_pollen)%>%
  kruskal_test(remaining_pollen~treatment)

df %>%
  drop_na(remaining_pollen)%>%
  dunn_test(remaining_pollen~treatment)


# prewinter weight ####

pre_winter_weight<-df %>%
  drop_na(prewinter_cocoon_weight)%>%
  ggplot(aes(x=treatment, y=prewinter_cocoon_weight, fill=treatment))+
  geom_boxplot(alpha=.6)+
  geom_point(position = position_jitter(w = 0.1, h = 0))
pre_winter_weight

ggsave(plot= pre_winter_weight, "output/pre_winter_weight.jpeg", width = 10, height = 12)


df %>%
  drop_na(prewinter_cocoon_weight)%>%
  kruskal_test(prewinter_cocoon_weight~treatment)



# not significant difference in pre winter weight but trending same way as everything else




# test of a cox model


a.cox <- coxph(Surv(time_to_event, status)~treatment, data=df)


ggsurvplot(surv_fit(a.cox), color="#2E9FDF", ggtheme = theme_minimal())













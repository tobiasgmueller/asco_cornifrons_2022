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
  filter(end_status_long != "removed_dead")

# then a table perchance?
# a little graffy waffy?


df %>%
  ggplot()+
  geom_bar(aes(x=end_status, fill=treatment))


# proportion death ####
prop_death <- df %>%
  filter(end_status_long=="dead")%>%
  ggplot()+
  geom_bar(aes(x=end_status_long, fill=treatment))+
  ggtitle("break down of death")+
  facet_wrap(~larva_stage)
prop_death
ggsave(plot= prop_death, "output/prop_death.jpeg")






# time to cocoon ####
def_to_cocoon<- df %>%
  ggplot(aes(x = treatment, y=def_to_cocoon, fill=treatment))+
  geom_boxplot()+
  geom_point(position = position_jitter(w = 0.1, h = 0), alpha = .5)+
  ylab("days from defecation until cocoon")+
  facet_wrap(~larva_stage)
def_to_cocoon
ggsave(plot= def_to_cocoon, "output/def_to_cocoon_facet.jpeg", width = 5, height = 6)

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
  geom_point(position = position_jitter(w = 0.2, h = .1))
remaining_pollen
ggsave(plot= remaining_pollen, "output/remaining_pollen.jpeg", width = 5, height = 6)

df %>%
  drop_na(remaining_pollen)%>%
  kruskal_test(remaining_pollen~treatment)

df %>%
  drop_na(remaining_pollen)%>%
  dunn_test(remaining_pollen~treatment,p.adjust.method = "holm")


# prewinter weight ####

pre_winter_weight<-df %>%
  drop_na(prewinter_cocoon_weight)%>%
  ggplot(aes(x=treatment, y=prewinter_cocoon_weight, fill=treatment))+
  geom_boxplot(alpha=.6)+
  geom_point(position = position_jitter(w = 0.1, h = 0))+
  facet_wrap(~larva_stage)
pre_winter_weight

ggsave(plot= pre_winter_weight, "output/pre_winter_weight_facet.jpeg", width = 5, height = 6)


df %>%
  drop_na(prewinter_cocoon_weight)%>%
  kruskal_test(prewinter_cocoon_weight~treatment)

df %>%
  drop_na(prewinter_cocoon_weight)%>%
  dunn_test(prewinter_cocoon_weight~treatment)


# not significant difference in pre winter weight but trending same way as everything else




# test of a cox model


a.cox <- coxph(Surv(time_to_event, status)~treatment, data=df)


ggsurvplot(surv_fit(a.cox), color="#2E9FDF", ggtheme = theme_minimal())













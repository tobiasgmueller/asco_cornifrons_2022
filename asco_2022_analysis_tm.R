# This is a loose first analysis of the 2022 asco-bacteriacide
# osmia cornifrosn experiments
# run in the spring of 2022


# first a few packages
library(tidyverse)
rm(list=ls()) # then clean the environment


# then read in the data
summary_df <- read_csv("input/asco_plate_maps_summary.csv")


# then a table perchance?
# a little graffy waffy?


summary_df %>%
  ggplot()+
  geom_bar(aes(x=status, fill=treatment))

# and then just looking at the dead ones
summary_df %>%
  filter(status=="dead")%>%
  ggplot()+
  geom_bar(aes(x=status, fill=treatment))


















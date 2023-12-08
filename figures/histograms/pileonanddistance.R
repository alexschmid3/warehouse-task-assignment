
#install.packages(c("tidyverse", "shiny")) #, "leaflet", "DT", "hash"))
library(shiny)
library(tidyverse)
library(dplyr)

#---------------------------------------------------------------#

#Read in data
data <- read_csv('pickdistrib.csv')

ddata <- read_csv('distdistrib.csv')

#---------------------------------------------------------------#

means <- data %>% group_by(picktype) %>% 
  summarise(mean_picks=mean(picks))

#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 10))

dev.new(width=8, height=4)

png(file="orderpickdistrib.png",
    width=500, height=600)

#Make the histogram
ggplot(data, aes(x=picks, fill=picktype)) +
  geom_density(adjust=3, alpha=0.5, color=NA) +
  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
  geom_vline(xintercept = deframe(means[1, 'mean_picks']), linetype="dashed", color = "#0072B2", size=1) +
  geom_vline(xintercept = deframe(means[2, 'mean_picks']), linetype="dashed", color = "#D55E00", size=1) +
  geom_vline(xintercept = 1, linetype="dashed", color = "#000000", size=1) +
  xlab(" Picks per pod ") +
  xlim(1,5) +
  ylab(" Density ") +
  labs(title="") +
  #theme(legend.position = c(0.4,.9), legend.direction = "horizontal")
  theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) 

dev.off()

#---------------------------------------------------------------#

#Group the pods by picks
bigpickmin = 3

aggdata_stg <- ddata %>%
  select(disttraveled, itempicks) %>%
  mutate(pickgroup = if_else(itempicks >= bigpickmin, "Multi-item picks", "Single-item picks"))

aggdata <- aggdata_stg %>%
  select(pickgroup, disttraveled)

means <- aggdata_stg %>% group_by(pickgroup) %>% 
  summarise(mean_distance=mean(disttraveled))

#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 10))

dev.new(width=8, height=4)

png(file="distdistrib.png",
    width=500, height=600)

#Make the histogram
ggplot(aggdata, aes(x=disttraveled, fill=pickgroup)) +
  geom_density(alpha=0.5, color=NA) +
  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
  geom_vline(xintercept = deframe(means[1, 'mean_distance']), linetype="dashed", color = "#0072B2", size=1) +
  geom_vline(xintercept = deframe(means[2, 'mean_distance']), linetype="dashed", color = "#D55E00", size=1) +
  xlab(" Pod distance travelled ") +
  ylab(" Density ") +
  #xlim(0,3) +
  labs(title="") +
  #theme(legend.position = c(0.4,.9), legend.direction = "horizontal")
  theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) 

dev.off()



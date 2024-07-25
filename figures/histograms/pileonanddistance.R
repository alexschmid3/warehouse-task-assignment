
#install.packages(c("tidyverse", "shiny")) #, "leaflet", "DT", "hash"))
library(shiny)
library(tidyverse)
library(dplyr)

#---------------------------------------------------------------#

#Read in data
data <- read_csv('pickdistrib_new.csv')
ddata <- read_csv('distdistrib_new.csv')

#---------------------------------------------------------------#

means <- data %>% group_by(picktype) %>% 
  summarise(mean_picks=mean(picks))
#means[2,2] = 1.96

#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 10))

dev.new(width=8, height=4)

png(file="orderpickdistrib_new.png",
    width=1000, height=1200)

#Make the histogram
ggplot(data, aes(x=picks, fill=picktype, color=picktype)) +
  geom_density(adjust=3, alpha=0.5, size=1.5) +
  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
  scale_colour_manual(values=c("#0072B2", "#D55E00")) +
  geom_vline(xintercept = deframe(means[1, 'mean_picks']), linetype="dashed", color = "#0072B2", size=2.5) +
  geom_vline(xintercept = deframe(means[2, 'mean_picks']), linetype="dashed", color = "#D55E00", size=2.5) +
  geom_vline(xintercept = 1, linetype="dashed", color = "#000000", size=2.5) +
  xlab(" Picks per pod ") +
  xlim(1,5) +
  ylab(" Density ") +
  labs(title="") +
  #theme(legend.position = c(0.4,.9), legend.direction = "horizontal")
  theme(legend.position="bottom", legend.text = element_text(size = 36)) +
  theme(axis.text = element_text(size = 30)) +
  theme(axis.title = element_text(size = 36)) 

dev.off()


#---------------------------------------------------------------#

#Group the pods by picks
aggdata_stg <- ddata %>%
  select(disttraveled, itempicks) %>%
  mutate(pickgroup = if_else(itempicks >= 4, "4-item picks", if_else(itempicks >= 3, "3-item picks", if_else(itempicks >= 2, "2-item picks","1-item picks"))))

aggdata <- aggdata_stg %>%
  select(pickgroup, disttraveled)

means <- aggdata_stg %>% group_by(pickgroup) %>% 
  summarise(mean_distance=mean(disttraveled))

#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 20))

dev.new(width=8, height=4)

png(file="distdistrib.png",
    width=1000, height=1200)

#Make the histogram
ggplot(aggdata, aes(x=disttraveled, fill=pickgroup, color=pickgroup)) +
  geom_density(alpha=0.3, size=1.5) +
  scale_fill_manual(values=c("#0072B2", "#533CBF", "#D55E00", "#FFB000")) +
  scale_colour_manual(values=c("#0072B2", "#533CBF", "#D55E00", "#FFB000")) +
  geom_vline(xintercept = deframe(means[1, 'mean_distance']), linetype="dashed", color = "#2D63F1", size=3) +
  geom_vline(xintercept = deframe(means[2, 'mean_distance']), linetype="dashed", color = "#533CBF", size=3) +
  geom_vline(xintercept = deframe(means[3, 'mean_distance']), linetype="dashed", color = "#FE6100", size=3) +
  geom_vline(xintercept = deframe(means[4, 'mean_distance']), linetype="dashed", color = "#FFB000", size=3) +
   xlab(" Pod distance travelled ") +
  ylab(" Density ") +
  xlim(0,300) +
  labs(title="") +
  #theme(legend.position = c(0.4,.9), legend.direction = "horizontal")
  theme(legend.position="bottom", legend.text = element_text(size = 36)) +
  theme(axis.text = element_text(size = 30)) +
  theme(axis.title = element_text(size = 36)) 

dev.off()



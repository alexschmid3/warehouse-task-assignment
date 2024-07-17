install.packages(c("tidyverse", "shiny", "leaflet", "DT", "hash"))
library(shiny)
library(DT)
library(tidyverse)
library(dplyr)

#---------------------------------------------------------------#
#---------------------------------------------------------------#

#Read in detour and delay distribution data
congdata <- read_csv('congestionsolution.csv')
nocongdata <- read_csv('nocongestionsolution.csv')

#---------------------------------------------------------------#
#---------------------------------------------------------------#

vec_cong <- rep("Congestion-aware", nrow(congdata))
vec_nocong <- rep("Congestion-blind", nrow(nocongdata))

congdata <- congdata %>%
  mutate(datalabel = vec_cong)
nocongdata <- nocongdata %>%
  mutate(datalabel = vec_nocong)

alldata = rbind(congdata, nocongdata)

#---------------------------------------------------------------#

#Detour
update_geom_defaults("text", list(size = 48))

dev.new(width=8, height=4)

png(file="congestionhistogram.png", width=900, height=1200)
# width=450, height=600)

alldata %>%
  ggplot(aes(x=congutil, color=datalabel, fill=datalabel)) +
  geom_density(alpha=0.3,size=1.5)+ 
  geom_vline(data = mean_detour, aes(xintercept = 1), size=1.5)+
  #xlim(0, 0.00000001)+
  #ylim(0, 100)+
  scale_colour_manual("", 
                      breaks = c("Congestion-aware", "Congestion-blind"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("Congestion-aware", "Congestion-blind"),
                    values = c("#DA853E", "#649BCB")) +
  labs(x="Congestion utilization (as fraction of capacity)", y="Density")+
  #scale_x_continuous(limits=c(0,0.15), labels = scales::percent_format(accuracy = 1)) +
  #theme(legend.position="none") +
  theme(legend.position="bottom", legend.text = element_text(size = 40)) +
  theme(axis.text = element_text(size = 40)) +
  theme(axis.title = element_text(size = 40))

dev.off()


#---------------------------------------------------------------#

#Detour
update_geom_defaults("text", list(size = 48))

dev.new(width=8, height=4)

png(file="congestionhistogram.png", width=900, height=1200)
# width=450, height=600)

alldata %>%
  ggplot(aes(x=congutil, color=datalabel, fill=datalabel)) +
  geom_histogram(alpha=0.3,size=1.5)+ 
  geom_vline(data = mean_detour, aes(xintercept = 1.05), size=1.5)+
  #xlim(0, 0.00000001)+
  #ylim(0, 100)+
  scale_colour_manual("", 
                      breaks = c("Congestion-aware", "Congestion-blind"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("Congestion-aware", "Congestion-blind"),
                    values = c("#DA853E", "#649BCB")) +
  labs(x="Congestion utilization", y="Intersection-time count")+
  #scale_x_continuous(limits=c(0,0.15), labels = scales::percent_format(accuracy = 1)) +
  #theme(legend.position="none") +
  theme(legend.position="bottom", legend.text = element_text(size = 40)) +
  theme(axis.text = element_text(size = 40)) +
  theme(axis.title = element_text(size = 40))

dev.off()

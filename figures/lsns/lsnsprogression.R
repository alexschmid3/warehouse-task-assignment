
#install.packages(c("tidyverse", "shiny")) #, "leaflet", "DT", "hash"))
library(shiny)
library(tidyverse)
library(dplyr)

#---------------------------------------------------------------#

#Read in data
lsnsdata <- read_csv('lsnsprogression.csv')
lsnsdata$scheme <- paste(lsnsdata$method, lsnsdata$initialization)

lsnsdata <- lsnsdata %>%
  filter(initialization == "(cold start)")


#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 10))

dev.new(width=8, height=4)

png(file="lsnsiter.png",
    width=1000, height=1200)

ggplot(lsnsdata, aes(x = lsnsiteration, y = objective, colour =scheme, group = scheme)) + 
  geom_line(size=2.5)+
  scale_color_manual(breaks = c("Learn-then-optimize (warm start)", "Learn-then-optimize (cold start)", "Learning-enhanced (warm start)", "Learning-enhanced (cold start)","Domain-based (warm start)", "Domain-based (cold start)","Randomized (warm start)", "Randomized (cold start)"),
                     values = c("#FFB385", "#FE6100", "#EF99C4", "#DC267F", "#B9CDFF", "#648FFF", "#A1A1A1", "#5A5A5A")) +
  xlim(0,50) +
  xlab(" \nLSNS iteration ") +
  ylab(" Throughput\n ") +
  labs(title="") +
  #theme(legend.position = c(0.4,.9), legend.direction = "horizontal")
  #theme(legend.position="bottom", legend.text = element_text(face="bold", size = 30)) +
  theme(legend.position="none")+
  theme(axis.text = element_text(size = 25)) +
  theme(axis.title = element_text(face="bold", size = 35)) +
  guides(color=guide_legend(nrow=4, byrow=TRUE))

dev.off()

#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 10))

dev.new(width=8, height=4)

png(file="lsnstime.png",
    width=1000, height=1200)

ggplot(lsnsdata, aes(x = cumul_time, y = objective, colour =scheme, group = scheme)) + 
  geom_line(size=2.5)+
  scale_color_manual(breaks = c("Learn-then-optimize (warm start)", "Learn-then-optimize (cold start)", "Learning-enhanced (warm start)", "Learning-enhanced (cold start)","Domain-based (warm start)", "Domain-based (cold start)","Randomized (warm start)", "Randomized (cold start)"),
                     values = c("#FFB385", "#FE6100", "#EF99C4", "#DC267F", "#B9CDFF", "#648FFF", "#A1A1A1", "#5A5A5A")) +
  xlim(0,400) +
  xlab(" \nTime (s) ") +
  ylab(" Throughput\n ") +
  labs(title="") +
  #theme(legend.position = c(0.4,.9), legend.direction = "horizontal")
  #theme(legend.position="bottom", legend.text = element_text(size = 24)) +
  theme(legend.position="none")+
  theme(axis.text = element_text(size = 25)) +
  theme(axis.title = element_text(face="bold", size = 35)) +
  guides(color=guide_legend(nrow=4, byrow=TRUE))

dev.off()

#---------------------------------------------------------------#

# breaks = c("Learn-then-optimize (warm start)", "Learn-then-optimize (cold start)", "Learning-enhanced (warm start)", "Learning-enhanced (cold start)","Domain-based (warm start)", "Domain-based (cold start)","Randomized (warm start)", "Randomized (cold start)"),
# values = c("#FE6100", "#FFB385", "#DC267F", "#EF99C4", "#648FFF", "#B9CDFF", "#5A5A5A", "#A1A1A1")

# breaks = c("Learn-then-optimize (warm start)", "Learn-then-optimize (cold start)", "Learning-enhanced (warm start)", "Learning-enhanced (cold start)","Domain-based (warm start)", "Domain-based (cold start)","Randomized (warm start)", "Randomized (cold start)"),
# values = c("#FFB385", "#FE6100", "#EF99C4", "#DC267F", "#B9CDFF", "#648FFF", "#A1A1A1", "#5A5A5A")



# Blue:    "#648FFF", "#B9CDFF"
# Purples: "#785EF0", "#B1A3F7"
# Pinks:   "#DC267F", "#EF99C4" 
# Oranges: "#FE6100", "#FFB385"
# Yellows: "#FFB000", "#FFE4A7"

# values = c("#E27532", "#2D5C87", "#DC267F"))

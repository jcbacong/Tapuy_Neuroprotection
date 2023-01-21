## Install packages for Kaplan-Meier plots
# install.packages("survivalAnalysis")
# install.packages("survminer")

## Import packages
library(survival)
library("survminer")
library(dplyr)
library("ggplot2")
library(ggsignif)
library(ggpubr)
library(rstatix)
library(cowplot)
library("lemon")
library(reprex)
library(AICcmodavg)
library(DescTools)
library("plotrix")
library(data.table)
library("rockchalk")
library(truncnorm)
library("BSDA")

## Set working directory to the file location using Rstudio API
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Prepare group name
rename_wine_or_lees = function(x) {
  if(substr(x,start=1,stop=2) == "FW"){
    return("Field Wine")
  } else if (substr(x,start=1,stop=2) == "FL") {
    return("Field Lees")
  } else if(substr(x,start=1,stop=2) == "LW") {
    return("Lab Wine")
  } else if(substr(x,start=1,stop=2) == "LL") {
    return("Lab Lees")
  } else {
    return("Control")
  }
}

rename_concentration = function(x) {
  if(substr(x,start=3,stop=3) == "L"){
    return("Low")
  } else if (substr(x,start=3,stop=3) == "M") {
    return("Mid")
  } else if(substr(x,start=3,stop=3) == "H") {
    return("High")
  } else {
    return(x)
  }
}

##########################################################################################################
########################################### SHORT TERM MEMORY ############################################
##########################################################################################################

## Load survival csv file
stmraw.csv = read.csv("Learning & Memory - STM Raw.csv", header = T)

## 1. Without toxin 
stmraw.csv$group = lapply(stmraw.csv$treatment, rename_wine_or_lees)
stmraw.csv[stmraw.csv == "Untreated"] = "OP50"
stmraw.csv$CI_index = (stmraw.csv$N_butanone - stmraw.csv$N_EtoH) / ((stmraw.csv$N_total - stmraw.csv$N_origin) + 1e-6)

## Factor categories
stmraw.csv$treatment <- factor(stmraw.csv$treatment, levels = c('OP50','EGCG','DMSO',  'Metab',
                                                                'FWL', 'FWM', 'FWH',
                                                                'LWL', 'LWM', 'LWH',
                                                                'FLL', 'FLM', 'FLH',
                                                                'LLL', 'LLM', 'LLH'))
stmraw.csv$group <- factor(stmraw.csv$group, levels = c('Control','Field Wine', 'Lab Wine',
                                                        'Field Lees', 'Lab Lees'))

# healthy.baseline.mean = mean(stmraw.csv$CI_index[stmraw.csv == 'U'])
# healthy.baseline.std = sd(stmraw.csv$CI_index[stmraw.csv == 'U'])
stmraw.csv$type = "(-) 12% NT"


## Visualize boxplots
color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")

## 2. With toxin
stmraw.toxin.csv = read.csv("Learning & Memory - STMTox Raw.csv", header = T)
stmraw.toxin.csv[stmraw.toxin.csv == "Untreated w/ toxin"] = "OP50"
stmraw.toxin.csv[stmraw.toxin.csv == "Untreated w/o toxin"] = "U"
stmraw.toxin.csv$group = lapply(stmraw.toxin.csv$treatment, rename_wine_or_lees)

stmraw.toxin.csv$CI_index = (stmraw.toxin.csv$N_butanone - stmraw.toxin.csv$N_EtoH) / ((stmraw.toxin.csv$N_total - stmraw.toxin.csv$N_origin) + 1e-6)

## Factor categories
stmraw.toxin.csv$treatment <- factor(stmraw.toxin.csv$treatment, levels = c('U', 'OP50','EGCG','DMSO',  'Metab',
                                                                            'FWL', 'FWM', 'FWH',
                                                                            'LWL', 'LWM', 'LWH',
                                                                            'FLL', 'FLM', 'FLH',
                                                                            'LLL', 'LLM', 'LLH'))
stmraw.toxin.csv$group <- factor(stmraw.toxin.csv$group, levels = c('Control','Field Wine', 'Lab Wine',
                                                                    'Field Lees', 'Lab Lees'))

# healthy.baseline.mean = mean(stmraw.toxin.csv$CI_index[stmraw.toxin.csv == 'U'])
# healthy.baseline.std = sd(stmraw.toxin.csv$CI_index[stmraw.toxin.csv == 'U'])

color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")
stmraw.toxin.csv = subset(stmraw.toxin.csv, treatment != 'U')
stmraw.toxin.csv$type = "(+) 12% NT"


## Short-term controls
short.term.csv = rbind(stmraw.csv, stmraw.toxin.csv)

## Box plot
short.term.summary <- short.term.csv %>% 
  group_by(treatment, group, type) %>%
  summarise(mean.CI_index = mean(CI_index),
            sem.CI_index = sd(CI_index),
            obs.raw = n()) %>%
  as.data.frame()

## Control only
short.term.summary.control = short.term.summary[short.term.summary$group == "Control",]

## Barplot
short.term.summary.control = short.term.summary.control[short.term.summary.control$treatment != "Metab",]
x = ggplot(short.term.summary.control, aes(fill=type, y=mean.CI_index, x=treatment)) +
  theme_classic() + labs(x = '', y = 'Chemotaxis Index (C.I)') +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0,1.0)) +
  scale_x_discrete(labels = c("OP50", "0.01% EGCG","0.1% DMSO")) +
  geom_errorbar(aes(ymin=mean.CI_index, ymax = mean.CI_index + (sem.CI_index/sqrt(obs.raw))), width=.2,position=position_dodge(.9)) +
  scale_fill_manual('',values=c("#f59541","#c750c1")) +
  theme(axis.text.x = element_text(color="black", size=12),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=13),
        axis.title.y = element_text(size=15)
        
  ) 
x

######## STATISTICS FOR SHORT TERM MEMORY ########
## Unpaired t-test between +/- 12% NT
## A. OP50
tox_norm1 <- c(rnorm(10^5, mean = round(subset(short.term.summary.control, treatment == "OP50" & type == "(+) 12% NT")[,"mean.CI_index"],2),
                        sd = round(subset(short.term.summary.control, treatment == "OP50" & type == "(+) 12% NT")[,"sem.CI_index"],2)
                    )) ## 0.16 +/- 0.23

notox_norm1 <- c(rnorm(10^5, mean = round(subset(short.term.summary.control, treatment == "OP50" & type == "(-) 12% NT")[,"mean.CI_index"],2),
                    sd = round(subset(short.term.summary.control, treatment == "OP50" & type == "(-) 12% NT")[,"sem.CI_index"],2)
                    ))  ## 0.81 +/- 0.01

t.test(tox_norm1, notox_norm1) ##  p < 0.0001
# t.test(subset(short.term.csv, treatment == "OP50" & type == "(+) 12% NT")[,"CI_index"],
#        subset(short.term.csv, treatment == "OP50" & type == "(-) 12% NT")[,"CI_index"])


## B. DMSO
# tox_norm2 <- c(rnorm(10^5, mean = round(subset(short.term.summary.control, treatment == "DMSO" & type == "(+) 12% NT")[,"mean.CI_index"],2),
#                     sd = round(subset(short.term.summary.control, treatment == "DMSO" & type == "(+) 12% NT")[,"sem.CI_index"],2)
# )) ## 0.46 +/- 0.-4
# 
# notox_norm2 <- c(rnorm(10^5, mean = round(subset(short.term.summary.control, treatment == "DMSO" & type == "(-) 12% NT")[,"mean.CI_index"],2),
#                       sd = round(subset(short.term.summary.control, treatment == "DMSO" & type == "(-) 12% NT")[,"sem.CI_index"],2)
# ))  ## 0.72 +/- 0.02
# t.test(tox_norm2, notox_norm2) ##  p < 0.0001
t.test(subset(short.term.csv, treatment == "DMSO" & type == "(+) 12% NT")[,"CI_index"],
       subset(short.term.csv, treatment == "DMSO" & type == "(-) 12% NT")[,"CI_index"])


## C. EGCG
t.test(subset(short.term.csv, treatment == "EGCG" & type == "(+) 12% NT")[,"CI_index"],
       subset(short.term.csv, treatment == "EGCG" & type == "(-) 12% NT")[,"CI_index"])
# tox_norm3 <- c(rnorm(10^5, mean = round(subset(short.term.summary.control, treatment == "EGCG" & type == "(+) 12% NT")[,"mean.CI_index"],2),
#                     sd = round(subset(short.term.summary.control, treatment == "EGCG" & type == "(+) 12% NT")[,"sem.CI_index"],2)
# )) ## 0.64 +/- 0.06
# 
# notox_norm3 <- c(rnorm(10^5, mean = round(subset(short.term.summary.control, treatment == "EGCG" & type == "(-) 12% NT")[,"mean.CI_index"],2),
#                       sd = round(subset(short.term.summary.control, treatment == "EGCG" & type == "(-) 12% NT")[,"sem.CI_index"],2)
# ))  ## 0.611 +/- 0.41
# t.test(tox_norm3, notox_norm3) ##  p < 0.0001


## Significance bars
## Only with ***
op50 <- tibble(
  x = c(0.77, 0.77, 1.23, 1.23),
  y = c(0.83, 0.87, 0.87, 0.38)
) ## p < 0.001

dmso <- tibble(
  x = c(2.77, 2.77, 3.23, 3.23),
  y = c(0.77, 0.83, 0.83, 0.52)
) ## p < 0.001


x1 = x + 
  annotate("text", x = 1.0, y = 0.88, label = "***", size = 6, color = "black") +
  geom_line(data = op50, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 3.0, y = 0.85, label = "***", size = 6, color = "black") +
  geom_line(data = dmso, aes(x= x, y = y, group=1), inherit.aes = F)

## Dunnet multiple comparison test among +12% NT
control.tox = subset(short.term.csv, type == "(+) 12% NT" & group == "Control" & treatment != "Metab")
summary(aov(formula = CI_index ~ treatment, data = control.tox))
DunnettTest(x=control.tox$CI_index, g=control.tox$treatment)

t.test(subset(short.term.csv, treatment == "EGCG" & type == "(+) 12% NT")[,"CI_index"],
       subset(short.term.csv, treatment == "FLL" & type == "(+) 12% NT")[,"CI_index"])

t.test(subset(short.term.csv, treatment == "EGCG" & type == "(+) 12% NT")[,"CI_index"],
       subset(short.term.csv, treatment == "LLL" & type == "(+) 12% NT")[,"CI_index"])

## What if -12% NT
control.notox = subset(short.term.csv, type == "(-) 12% NT" & group == "Control" & treatment != "Metab")
summary(aov(formula = CI_index ~ treatment, data = control.notox))
DunnettTest(x=control.notox$CI_index, g=control.notox$treatment)

egcg <- tibble(
  x = c(1.23, 1.23, 2.23, 2.23),
  y = c(0.91, 0.98, 0.98, 0.72)
) ## p < 0.009 **

x1 = x1 + 
  annotate("text", x = 1.75, y = 1., label = "**", size = 6, color = "red") +
  geom_line(data = egcg, aes(x= x, y = y, group=1), inherit.aes = F) 
x1


png(filename = file.path(getwd(),"/Paper-short-term.png"), width = 6, height = 4, units = "in", res = 600)
# short.term.boxplot = ggarrange(plotlist = list(short.term.x, short.term.y), ncol = 1, nrow = 2,
#                                labels = c("A", "B"))
plot(x1)
dev.off()
print(paste("Saved at ", getwd()))


## Tapuy extracts
treatment.tox = subset(short.term.csv, type == "(+) 12% NT" &  !(treatment %in% c('Metab', 'DMSO')) )
summary(aov(formula = CI_index ~ treatment, data = treatment.tox))
DunnettTest(x=treatment.tox$CI_index, g=treatment.tox$treatment)

short.term.summary.treatment = subset(short.term.summary, !(treatment %in% c('Metab', 'DMSO')) & type == "(+) 12% NT")

x2 = ggplot(short.term.summary.treatment, aes(y=mean.CI_index, x=treatment, fill=group)) +
  theme_classic() + labs(x = '', y = 'Chemotaxis Index (C.I)') +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0,1.0)) +
  scale_x_discrete(labels = c("OP50", "0.01%\nEGCG",
                              "Low", "Mid", "High",
                              "Low", "Mid", "High",
                              "Low", "Mid", "High",
                              "Low", "Mid", "High")) +
  geom_errorbar(aes(ymin=mean.CI_index, ymax = mean.CI_index + (sem.CI_index/sqrt(obs.raw))), width=.2,position=position_dodge(.9)) +
  # scale_fill_manual('',values=c("#c750c1", "#f7f00c","#3392d1", "#84d94f","#f25c78" )) +
  # scale_fill_manual('',values=c("#c750c1",  "#00888b","#2a4858","#f5fa6e","#61c888")) +
  # scale_fill_manual('',values=c("#c750c1", "#e65089","#ff736b", "#ffa451","#ffd84f" )) +
  scale_fill_manual('',values=c("#c750c1", "#e65089","#ff9e5f", "#ff7271","#ffcd5f" )) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(color="black", size=11),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        axis.title.y = element_text(size=15),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=13)
        
  ) 
x2
x2 = x2 + 
  annotate("text", x = 2, y = 0.70, label = "**", size = 6, color = "red") +
  annotate("text", x = 9, y = 0.87, label = "**", size = 6, color = "red") +
  annotate("text", x = 12, y = 0.83, label = "**", size = 6, color = "red") 
x2


## Control
short.term.summary.control= subset(short.term.summary, (treatment %in% c('OP50', 'EGCG')) & type == "(+) 12% NT")
x2 = ggplot(short.term.summary.control, aes(y=mean.CI_index, x=treatment, fill=group)) +
  theme_classic() + labs(x = '', y = 'Chemotaxis Index (C.I)',title="Control") +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0,1.0)) +
  scale_x_discrete(labels = c("OP50", "0.01%\nEGCG")) +
  geom_errorbar(aes(ymin=mean.CI_index, ymax = mean.CI_index + (sem.CI_index/sqrt(obs.raw))), width=.2,position=position_dodge(.9)) +
  # scale_fill_manual('',values=c("#c750c1", "#f7f00c","#3392d1", "#84d94f","#f25c78" )) +
  # scale_fill_manual('',values=c("#c750c1",  "#00888b","#2a4858","#f5fa6e","#61c888")) +
  # scale_fill_manual('',values=c("#c750c1", "#e65089","#ff736b", "#ffa451","#ffd84f" )) +
  scale_fill_manual('',values=c("#c750c1")) +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color="black", size=11),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        axis.title.y = element_text(size=15)
  ) +
  annotate("text", x = 2, y = 0.70, label = "**", size = 6, color = "red") 
x2

short.term.summary.fw= subset(short.term.summary, (treatment %in% c('FWL', 'FWM', 'FWH')) & type == "(+) 12% NT")
x3 = ggplot(short.term.summary.fw, aes(y=mean.CI_index, x=treatment, fill=group)) +
  theme_classic() + labs(x = '', y = '',title="Field Wine") +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "", "" ),
                     limits = c(0,1.0)) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  geom_errorbar(aes(ymin=mean.CI_index, ymax = mean.CI_index + (sem.CI_index/sqrt(obs.raw))), width=.2,position=position_dodge(.9)) +
  # scale_fill_manual('',values=c("#c750c1", "#f7f00c","#3392d1", "#84d94f","#f25c78" )) +
  # scale_fill_manual('',values=c("#c750c1",  "#00888b","#2a4858","#f5fa6e","#61c888")) +
  # scale_fill_manual('',values=c("#c750c1", "#e65089","#ff736b", "#ffa451","#ffd84f" )) +
  scale_fill_manual('',values=c("#e65089")) +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color="black", size=11),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        axis.title.y = element_text(size=15)
  ) 
x3

short.term.summary.lw= subset(short.term.summary, (treatment %in% c('LWL', 'LWM', 'LWH')) & type == "(+) 12% NT")
x4 = ggplot(short.term.summary.lw, aes(y=mean.CI_index, x=treatment, fill=group)) +
  theme_classic() + labs(x = '', y = '',title="Lab Wine") +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "", "" ),
                     limits = c(0,1.0)) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  geom_errorbar(aes(ymin=mean.CI_index, ymax = mean.CI_index + (sem.CI_index/sqrt(obs.raw))), width=.2,position=position_dodge(.9)) +
  # scale_fill_manual('',values=c("#c750c1", "#f7f00c","#3392d1", "#84d94f","#f25c78" )) +
  # scale_fill_manual('',values=c("#c750c1",  "#00888b","#2a4858","#f5fa6e","#61c888")) +
  # scale_fill_manual('',values=c("#c750c1", "#e65089","#ff736b", "#ffa451","#ffd84f" )) +
  scale_fill_manual('',values=c("#ff9e5f")) +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color="black", size=11),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        axis.title.y = element_text(size=15)
  ) 
x4

short.term.summary.fl= subset(short.term.summary, (treatment %in% c('FLL', 'FLM', 'FLH')) & type == "(+) 12% NT")
x5 = ggplot(short.term.summary.fl, aes(y=mean.CI_index, x=treatment, fill=group)) +
  theme_classic() + labs(x = '', y = '',title="Field Lees") +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "", "" ),
                     limits = c(0,1.0)) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  geom_errorbar(aes(ymin=mean.CI_index, ymax = mean.CI_index + (sem.CI_index/sqrt(obs.raw))), width=.2,position=position_dodge(.9)) +
  # scale_fill_manual('',values=c("#c750c1", "#f7f00c","#3392d1", "#84d94f","#f25c78" )) +
  # scale_fill_manual('',values=c("#c750c1",  "#00888b","#2a4858","#f5fa6e","#61c888")) +
  # scale_fill_manual('',values=c("#c750c1", "#e65089","#ff736b", "#ffa451","#ffd84f" )) +
  scale_fill_manual('',values=c("#ff7271")) +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color="black", size=11),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        axis.title.y = element_text(size=15)
  ) +
  annotate("text", x = 1, y = 0.87, label = "**", size = 6, color = "red") 
x5

short.term.summary.ll = subset(short.term.summary, (treatment %in% c('LLL', 'LLM', 'LLH')) & type == "(+) 12% NT")
x6 = ggplot(short.term.summary.ll, aes(y=mean.CI_index, x=treatment, fill=group)) +
  theme_classic() + labs(x = '', y = '',title="Lab Lees") +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "", "" ),
                     limits = c(0,1.0)) +
  scale_x_discrete(labels = c("L", "M", "H")) +
  geom_errorbar(aes(ymin=mean.CI_index, ymax = mean.CI_index + (sem.CI_index/sqrt(obs.raw))), width=.2,position=position_dodge(.9)) +
  # scale_fill_manual('',values=c("#c750c1", "#f7f00c","#3392d1", "#84d94f","#f25c78" )) +
  # scale_fill_manual('',values=c("#c750c1",  "#00888b","#2a4858","#f5fa6e","#61c888")) +
  # scale_fill_manual('',values=c("#c750c1", "#e65089","#ff736b", "#ffa451","#ffd84f" )) +
  scale_fill_manual('',values=c("#ffcd5f")) +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color="black", size=11),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        axis.title.y = element_text(size=15)
  ) +
  annotate("text", x = 1, y = 0.83, label = "**", size = 6, color = "red") 
x6
combplot = plot_grid(x2, x3, x4, x5, x6, align = "h", ncol = 5,rel_widths = c(0.20, 0.20, 0.20, 0.20, 0.20 ))
combplot






png(filename = file.path(getwd(),"/Paper-short-term-treatment.png"), width = 8, height = 4, units = "in", res = 600)
# short.term.boxplot = ggarrange(plotlist = list(short.term.x, short.term.y), ncol = 1, nrow = 2,
#                                labels = c("A", "B"))
plot(x2)
dev.off()
print(paste("Saved at ", getwd()))


##########################################################################################################
########################################### LONG TERM MEMORY #############################################
##########################################################################################################

## Load survival csv file
ltmraw.csv = read.csv("Learning & Memory - LTM Raw.csv", header = T)


## 1. Without toxin 
ltmraw.csv$group = lapply(ltmraw.csv$treatment, rename_wine_or_lees)
ltmraw.csv[ltmraw.csv == "Untreated"] = "OP50"
ltmraw.csv$CI_index = (ltmraw.csv$N_butanone - ltmraw.csv$N_EtoH) / ((ltmraw.csv$N_total - ltmraw.csv$N_origin) + 1e-6)

## Factor categories
ltmraw.csv$treatment <- factor(ltmraw.csv$treatment, levels = c('OP50', 'EGCG', 'DMSO','Metab',
                                                                'FWL', 'FWM', 'FWH',
                                                                'LWL', 'LWM', 'LWH',
                                                                'FLL', 'FLM', 'FLH',
                                                                'LLL', 'LLM', 'LLH'))
ltmraw.csv$group <- factor(ltmraw.csv$group, levels = c('Control','Field Wine', 'Lab Wine',
                                                        'Field Lees', 'Lab Lees'))

# healthy.baseline.mean = mean(ltmraw.csv$CI_index[ltmraw.csv == 'U'])
# healthy.baseline.std = sd(ltmraw.csv$CI_index[ltmraw.csv == 'U'])
ltmraw.csv$type = "(-) 12% NT"

## Visualize boxplots
color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")


## 2. With toxin
ltmraw.toxin.csv = read.csv("Learning & Memory - LTMTox Raw.csv", header = T)
ltmraw.toxin.csv[ltmraw.toxin.csv == "Untreated w/ toxin"] = "OP50"
ltmraw.toxin.csv[ltmraw.toxin.csv == "Untreated w/o toxin"] = "U"
ltmraw.toxin.csv$group = lapply(ltmraw.toxin.csv$treatment, rename_wine_or_lees)

ltmraw.toxin.csv$CI_index = (ltmraw.toxin.csv$N_butanone - ltmraw.toxin.csv$N_EtoH) / ((ltmraw.toxin.csv$N_total - ltmraw.toxin.csv$N_origin) + 1e-6)

## Factor categories
ltmraw.toxin.csv$treatment <- factor(ltmraw.toxin.csv$treatment, levels = c('U','OP50','EGCG', 'DMSO',  'Metab',
                                                                            'FWL', 'FWM', 'FWH',
                                                                            'LWL', 'LWM', 'LWH',
                                                                            'FLL', 'FLM', 'FLH',
                                                                            'LLL', 'LLM', 'LLH'))
ltmraw.toxin.csv$group <- factor(ltmraw.toxin.csv$group, levels = c('Control','Field Wine', 'Lab Wine',
                                                                    'Field Lees', 'Lab Lees'))

# healthy.baseline.mean = mean(ltmraw.toxin.csv$CI_index[ltmraw.toxin.csv == 'U'])
# healthy.baseline.std = sd(ltmraw.toxin.csv$CI_index[ltmraw.toxin.csv == 'U'])

color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")
ltmraw.toxin.csv = subset(ltmraw.toxin.csv, treatment != 'U')
ltmraw.toxin.csv$type = "(+) 12% NT"

## Long-term controls
long.term.csv = rbind(ltmraw.csv, ltmraw.toxin.csv)

## Box plot
long.term.summary <- long.term.csv %>% 
  group_by(treatment, group, type) %>%
  summarise(mean.CI_index = mean(CI_index),
            sem.CI_index = sd(CI_index),
            obs.raw = n()) %>%
  as.data.frame()

## Control only
long.term.summary.control = long.term.summary[long.term.summary$group == "Control",]

## Barplot
long.term.summary.control = long.term.summary.control[long.term.summary.control$treatment != "Metab",]
y = ggplot(long.term.summary.control, aes(fill=type, y=mean.CI_index, x=treatment)) +
  theme_classic() + labs(x = '', y = 'Chemotaxis Index (C.I)') +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0,1.0)) +
  scale_x_discrete(labels = c("OP50", "0.01% EGCG", "0.1% DMSO")) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean.CI_index, ymax = mean.CI_index + sem.CI_index/sqrt(obs.raw)), width=.2,position=position_dodge(.9)) +
  scale_fill_manual('',values=c("#f59541","#c750c1")) +
  theme(axis.text.x = element_text(color="black", size=12),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=13),
        axis.title.y = element_text(size=15)
        
  ) 
y

######## STATISTICS FOR SHORT TERM MEMORY ########
## Unpaired t-test between +/- 12% NT
## A. OP50
tox_norm1 <- c(rnorm(10^5, mean = round(subset(long.term.summary.control, treatment == "OP50" & type == "(+) 12% NT")[,"mean.CI_index"],2),
                     sd = round(subset(long.term.summary.control, treatment == "OP50" & type == "(+) 12% NT")[,"sem.CI_index"],2)
)) ## 0.16 +/- 0.23

notox_norm1 <- c(rnorm(10^5, mean = round(subset(long.term.summary.control, treatment == "OP50" & type == "(-) 12% NT")[,"mean.CI_index"],2),
                       sd = round(subset(long.term.summary.control, treatment == "OP50" & type == "(-) 12% NT")[,"sem.CI_index"],2)
))  ## 0.81 +/- 0.01
# t.test(tox_norm1, notox_norm1) ##  p < 0.0001
t.test(subset(long.term.csv, treatment == "OP50" & type == "(+) 12% NT")[,"CI_index"],
       subset(long.term.csv, treatment == "OP50" & type == "(-) 12% NT")[,"CI_index"])


## B. DMSO
# tox_norm2 <- c(rnorm(10^5, mean = round(subset(short.term.summary.control, treatment == "DMSO" & type == "(+) 12% NT")[,"mean.CI_index"],2),
#                     sd = round(subset(short.term.summary.control, treatment == "DMSO" & type == "(+) 12% NT")[,"sem.CI_index"],2)
# )) ## 0.46 +/- 0.-4
# 
# notox_norm2 <- c(rnorm(10^5, mean = round(subset(short.term.summary.control, treatment == "DMSO" & type == "(-) 12% NT")[,"mean.CI_index"],2),
#                       sd = round(subset(short.term.summary.control, treatment == "DMSO" & type == "(-) 12% NT")[,"sem.CI_index"],2)
# ))  ## 0.72 +/- 0.02
# t.test(tox_norm2, notox_norm2) ##  p < 0.0001
t.test(subset(short.term.csv, treatment == "DMSO" & type == "(+) 12% NT")[,"CI_index"],
       subset(short.term.csv, treatment == "DMSO" & type == "(-) 12% NT")[,"CI_index"])


## C. EGCG
t.test(subset(short.term.csv, treatment == "EGCG" & type == "(+) 12% NT")[,"CI_index"],
       subset(short.term.csv, treatment == "EGCG" & type == "(-) 12% NT")[,"CI_index"])

dmso <- tibble(
  x = c(2.77, 2.77, 3.23, 3.23),
  y = c(0.78, 0.83, 0.83, 0.55)
) ## p < 0.01


y1 = y + 
  annotate("text", x = 3.0, y = 0.85, label = "**", size = 6, color = "black") +
  geom_line(data = dmso, aes(x= x, y = y, group=1), inherit.aes = F)
y1


## Dunnet multiple comparison test among +12% NT
control.tox = subset(long.term.csv, type == "(+) 12% NT" & group == "Control" & treatment != "Metab")
summary(aov(formula = CI_index ~ treatment, data = control.tox))
DunnettTest(x=control.tox$CI_index, g=control.tox$treatment) ## not significant

## What if -12% NT
control.notox = subset(long.term.csv, type == "(-) 12% NT" & group == "Control" & treatment != "Metab")
summary(aov(formula = CI_index ~ treatment, data = control.notox))
DunnettTest(x=control.notox$CI_index, g=control.notox$treatment) ## not significant

png(filename = file.path(getwd(),"/Paper-long-term.png"), width = 6, height = 4, units = "in", res = 600)
# short.term.boxplot = ggarrange(plotlist = list(short.term.x, short.term.y), ncol = 1, nrow = 2,
#                                labels = c("A", "B"))
plot(y1)
dev.off()
print(paste("Saved at ", getwd()))

## Tapuy extracts
treatment.tox = subset(long.term.csv, type == "(+) 12% NT" &  !(treatment %in% c('Metab', 'DMSO')) )
summary(aov(formula = CI_index ~ treatment, data = treatment.tox))
DunnettTest(x=treatment.tox$CI_index, g=treatment.tox$treatment)

long.term.summary.treatment = subset(long.term.summary, !(treatment %in% c('Metab', 'DMSO')) & type == "(+) 12% NT")

y2 = ggplot(long.term.summary.treatment, aes(y=mean.CI_index, x=treatment, fill=group)) +
  theme_classic() + labs(x = '', y = 'Chemotaxis Index (C.I)') +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0,1.0)) +
  scale_x_discrete(labels = c("OP50", "0.01%\nEGCG",
                              "Low", "Mid", "High",
                              "Low", "Mid", "High",
                              "Low", "Mid", "High",
                              "Low", "Mid", "High")) +
  geom_errorbar(aes(ymin=mean.CI_index, ymax = mean.CI_index + (sem.CI_index/sqrt(obs.raw))), width=.2,position=position_dodge(.9)) +
  scale_fill_manual('',values=c("#c750c1", "#e65089","#ff9e5f", "#ff7271","#ffcd5f" )) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(color="black", size=11),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        axis.title.y = element_text(size=15),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=13)
        
  ) 
y2

# x2 = x2 + 
#   annotate("text", x = 2, y = 0.70, label = "**", size = 6, color = "red") +
#   annotate("text", x = 9, y = 0.87, label = "**", size = 6, color = "red") +
#   annotate("text", x = 12, y = 0.83, label = "**", size = 6, color = "red") 
# x2

png(filename = file.path(getwd(),"/Paper-long-term-treatment.png"), width = 8, height = 4, units = "in", res = 600)
# short.term.boxplot = ggarrange(plotlist = list(short.term.x, short.term.y), ncol = 1, nrow = 2,
#                                labels = c("A", "B"))
plot(y2)
dev.off()
print(paste("Saved at ", getwd()))

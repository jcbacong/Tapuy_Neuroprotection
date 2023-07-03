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
library(ggtext)
library(grid)


## Set working directory to the file location using Rstudio API
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Prepare group name
rename_wine_or_lees = function(x) {
  if(substr(x,start=1,stop=2) == "FW"){
    return("Wine")
  } else if (substr(x,start=1,stop=2) == "FL") {
    return("Lees")
  } else if(substr(x,start=1,stop=2) == "LW") {
    return("Wine")
  } else if(substr(x,start=1,stop=2) == "LL") {
    return("Lees")
  } else {
    return(x)
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
    return("Control")
  }
}

t_test_means = function(df.norm.6h, X, Y="OP50", alternative="less") {
  x.mean = df.norm.6h[df.norm.6h$Treatment==X, "mean.CI"]
  x.sem = df.norm.6h[df.norm.6h$Treatment==X,"std.CI"]
  x.obs = 3
  
  y.mean = df.norm.6h[df.norm.6h$Treatment==Y,"mean.CI"]
  y.sem = df.norm.6h[df.norm.6h$Treatment==Y,"std.CI"]
  y.obs = 3
  
  x = tsum.test(x.mean, x.sem, x.obs,
                y.mean, y.sem, y.obs,
                alternative = alternative)
  return(x)
}

stat_results_df = function(df.t, Y="OP50", alternative="less") {
  results_df <- data.frame()
  for (x in df.t$Treatment) {
    if(x != Y) {
      t_test_result = t_test_means(df.t, x, Y=Y, alternative=alternative)
      p_value <- t_test_result$p.value
      statistics <- t_test_result$statistic
      estimate <- t_test_result$estimate
      new_row <- data.frame(group = paste(x,"-", Y),  t = statistics, P_Value = p_value)
      results_df <- rbind(results_df, new_row)
    }
  }
  return(results_df)
}


##########################################################################################################
########################################### LEARNING and MEMORY  ############################################
##########################################################################################################

## Load survival csv file
learn.csv = read.csv("Learning & Memory - Learning.csv", header = T)
memory.csv = read.csv("Learning & Memory - Memory.csv", header = T)

## Learning
learn.csv$CI_index = (learn.csv$N_butanone - learn.csv$N_EtoH) / ((learn.csv$N_total - learn.csv$N_origin) + 1e-6)
learn.csv$Type = "Learning"

## Memory
memory.csv$CI_index = (memory.csv$N_butanone - memory.csv$N_EtoH) / ((memory.csv$N_total - memory.csv$N_origin) + 1e-6)
memory.csv$Type = "Memory"

## Summarize
lam.csv = rbind(learn.csv,memory.csv)

lam.csv$group = lapply(lam.csv$Treatment, rename_wine_or_lees)
lam.csv$group = factor(lam.csv$group, 
                       levels = c("OP50", "DMSO", "EGCG","Wine", "Lees"),
                       labels = c("OP50", "DMSO", "EGCG","Wine", "Lees"))

lam.csv$Treatment = factor(lam.csv$Treatment, 
                           levels = c("OP50", "DMSO", "EGCG","FWH", "LWH", "FLH", "LLH"),
                           labels = c("OP50", "DMSO", "EGCG","Field Wine", "Lab Wine", "Field Lees", "Lab Lees"))

lam.csv$Type = factor(lam.csv$Type,
                      levels=c("Learning", "Memory"),
                      labels=c("Learning", "Memory"))
lam.csv$Toxin = factor(lam.csv$Toxin,
                       levels=c("without toxin", "with toxin"),
                       labels=c("-12% NT", "+12% NT"))


lam.df <- lam.csv %>% 
  group_by(Treatment, group, Type, Toxin) %>%
  summarise(mean.CI = mean(CI_index) + 1e-6,
            std.CI = std.error(CI_index),
            count = n()) %>%
  as.data.frame()
lam.df

write.csv(lam.df, "L&M-values.csv")
lam.df = subset(lam.df, Treatment != "DMSO")

## Statistics
##########################################################################################################
########################################### LEARNING############################################
##########################################################################################################
## 1. Between -12 NT and +12 NT
dat.tox = subset(lam.df, Treatment == "OP50" & Type == "Learning")
## x = -12NT, y = +12NT, alternative = greater
x.mean = dat.tox[dat.tox$Toxin == "-12% NT", "mean.CI"]
x.std = dat.tox[dat.tox$Toxin == "-12% NT", "std.CI"]
x.obs= dat.tox[dat.tox$Toxin == "-12% NT", "count"]

y.mean = dat.tox[dat.tox$Toxin == "+12% NT", "mean.CI"]
y.std = dat.tox[dat.tox$Toxin == "+12% NT", "std.CI"]
y.obs= dat.tox[dat.tox$Toxin == "+12% NT", "count"]

res.learn.tox = tsum.test(x.mean, x.std, x.obs,
                          y.mean, y.std, y.obs,
                          alternative = "greater")
res.learn.tox


### Visualization
color_scheme = c("#f25c78",  "#84d94e","#ffa600", "#b64cdf")
lam.df.learn = subset(lam.df, Type == "Learning")
lam.df.learn$mean.CI[4] = subset(subset(learn.csv, Treatment == "FWH"), CI_index != 0)$CI_index
lam.df.learn$name = paste(lam.df.learn$Treatment, lam.df.learn$Toxin)
lam.df.learn$name = factor(lam.df.learn$name, 
                           levels = c("OP50 -12% NT", "OP50 +12% NT",
                                      "EGCG +12% NT",
                                      "Field Wine +12% NT","Lab Wine +12% NT",
                                      "Field Lees +12% NT", "Lab Lees +12% NT"))

## 2. Across treatment
dat.treat = subset(lam.csv, Toxin == "+12% NT" & Type == "Learning")
res.learn.treat = DunnettTest(x=dat.treat$CI_index, g=dat.treat$Treatment)
# res.learn.treat = stat_results_df(subset(lam.df.learn, Toxin !="-12% NT"), alternative = "greater")
res.learn.treat

p = ggplot(lam.df.learn, aes(name,mean.CI,fill=group)) +
  theme_classic() + labs(x = '', y = 'Learning Index') + 
  theme(axis.title.y = element_text(size = 14), axis.text.x=element_text(size=13), axis.text.y=element_text(size=12)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI + std.CI), width=.2,position=position_dodge(.9)) +
  # facet_grid('. ~ Toxin', scales="free", switch = "x", space='free') +
  scale_fill_manual(values=color_scheme) +
  theme(legend.position="none") + 
  theme(strip.placement = "outside") +
  theme(axis.title.x = element_blank(),
        axis.title.y =element_text(size=16),
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_markdown(size=10),
        axis.text.y = element_text(size=11),
        legend.position="none",
        strip.text.x = element_text(size=13, face = "bold"),
        strip.background = element_blank()) +
  scale_x_discrete(labels=c("OP50 <br><i>E. coli</i>", "OP50 <br><i>E. coli</i>",
                            "200 uM <br> EGCG",
                            "Field Wine", "Lab Wine", "Field Lees",
                            "Lab Lees")) +
  theme(plot.margin = unit(c(0,0,2,0), "lines")) +
  coord_cartesian(ylim=c(0,1), clip="off") 
p

tox.vs.notox <- tibble(
  x = c(1, 1, 2,2),
  y = c(0.84, 0.92, 0.92, 0.62))

p = p + 
  annotate("text", x = 1.5, y = 0.94, label = "*", size = 7, color = "black") +
  geom_line(data = tox.vs.notox, aes(x= x, y = y, group=1), inherit.aes = F) 
  # annotate("text", x = 3, y = 0.687, label = "***", size = 6, color = "#960101")  
  # annotate("text", x = 6, y = 0.475, label = "*", size = 6, color = "#960101") +
  # annotate("text", x = 7, y = 0.519, label = "**", size = 6, color = "#960101")
p



##########################################################################################################
########################################### MEMORY  ############################################
##########################################################################################################

## 1. Between -12 NT and +12 NT
dat.tox = subset(lam.df, Treatment == "OP50" & Type == "Memory")
## x = -12NT, y = +12NT, alternative = greater
x.mean = dat.tox[dat.tox$Toxin == "-12% NT", "mean.CI"]
x.std = dat.tox[dat.tox$Toxin == "-12% NT", "std.CI"]
x.obs= dat.tox[dat.tox$Toxin == "-12% NT", "count"]

y.mean = dat.tox[dat.tox$Toxin == "+12% NT", "mean.CI"]
y.std = dat.tox[dat.tox$Toxin == "+12% NT", "std.CI"]
y.obs= dat.tox[dat.tox$Toxin == "+12% NT", "count"]

res.mem.tox = tsum.test(x.mean, x.std, x.obs,
                          y.mean, y.std, y.obs,
                          alternative = "greater")
res.mem.tox ## p = 0.1535

## 2. Across treatment
dat.treat = subset(lam.csv, Toxin == "+12% NT" & Type == "Memory")
res.mem.treat = DunnettTest(x=dat.treat$CI_index, g=dat.treat$Treatment)
res.mem.treat


### Visualization
lam.df.mem = subset(lam.df, Type == "Memory")
lam.df.mem$name = paste(lam.df.mem$Treatment, lam.df.mem$Toxin)
lam.df.mem$name = factor(lam.df.mem$name, 
                           levels = c("OP50 -12% NT", "OP50 +12% NT",
                                      "EGCG +12% NT",
                                      "Field Wine +12% NT","Lab Wine +12% NT",
                                      "Field Lees +12% NT", "Lab Lees +12% NT"))

# res.mem.treat = stat_results_df(subset(lam.df.mem, Toxin !="-12% NT"), alternative = "greater")
# res.mem.treat

q = ggplot(lam.df.mem, aes(name,mean.CI,fill=group)) +
  theme_classic() + labs(x = '', y = 'Memory Index') + 
  theme(axis.title.y = element_text(size = 14), axis.text.x=element_text(size=13), axis.text.y=element_text(size=12)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI + std.CI), width=.2,position=position_dodge(.9)) +
  # facet_grid('. ~ Toxin', scales="free", switch = "x", space='free') +
  scale_fill_manual(values=color_scheme) +
  theme(legend.position="none") + 
  theme(strip.placement = "outside") +
  theme(axis.title.x = element_blank(),
        axis.title.y =element_text(size=16),
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_markdown(size=10),
        axis.text.y = element_text(size=11),
        legend.position="none",
        strip.text.x = element_text(size=13, face = "bold"),
        strip.background = element_blank()) +
  scale_x_discrete(labels=c("OP50 <br><i>E. coli</i>", "OP50 <br><i>E. coli</i>",
                            "200 uM <br> EGCG",
                            "Field Wine", "Lab Wine", "Field Lees",
                            "Lab Lees")) +
  theme(plot.margin = unit(c(0,0,2,0), "lines")) +
  coord_cartesian(xlim=c(1,7), ylim=c(0,1), clip="off") 
q

tox.vs.notox <- tibble(
  x = c(1, 1, 2,2),
  y = c(0.77, 0.88, 0.88, 0.75))

q = q + 
  annotate("text", x = 1.5, y = 0.93, label = "n.s", size = 5, color = "black") +
  geom_line(data = tox.vs.notox, aes(x= x, y = y, group=1), inherit.aes = F)  
  # annotate("segment", x = 5.5, xend = 8, y = -0.21, yend = -0.21) +
  # annotate("text", x = 5, y = -0.21, label = "+ 12% NT", fontface =2, size=5)  +
  # annotate("text", x = 3, y = 0.98, label = "***", size=6, color = "#960101")
q

ggsave(plot=p, filename="Learning.png", width=8, height=4, dpi=600)
# png(filename = file.path(getwd(),"/Learning.png"), width = 8, height = 3, units = "in", res = 600)
# plot(p)
# dev.off()
print(paste("Saved at ", getwd()))

ggsave(plot=q, filename="Memory.png", width=8, height=4, dpi=600)
# png(filename = file.path(getwd(),"/Memory.png"), width = 8, height = 3, units = "in", res = 600)
# plot(q)
# dev.off()
print(paste("Saved at ", getwd()))

all = ggarrange(plotlist = list(p, q), ncol = 1, nrow = 2, labels = c("A", "B"))
ggsave(plot=all, filename="Combined.png", width=8, height=4, dpi=600)








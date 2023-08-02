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


####### WITH THE TREATMENT #########
## Load survival csv file
df.csv = read.csv("With Glyphosate + Treatment - Survival Assay 3.csv", header = T)
df.csv = subset(df.csv, Treatment != "DMSO")

df.csv <- df.csv %>%
  mutate(X = case_when(
    grepl("DMSO|EGCG|E\\. coli", Treatment) ~ Treatment,
    grepl("Lab Wine|Field Wine", Treatment) ~ "Wine",
    grepl("Lab Lees|Field Lees", Treatment) ~ "Lees"
  ),
  Y = case_when(
    grepl("DMSO|EGCG|E\\. coli", Treatment) ~ "Control",
    grepl("Lab Wine|Field Wine", Treatment) ~ sub(".* ", "", Treatment),
    grepl("Lab Lees|Field Lees", Treatment) ~ sub(".* ", "", Treatment)
  ))

## Stat
df.res = subset(newdat, Toxin == "+12% NT")
df.res$lifetime = df.res$Time -1
res = DunnettTest(x=df.res$lifetime, g=df.res$Treatment)
print(res)


newdat$name = paste(newdat$Treatment, newdat$Toxin)
newdat$name = factor(newdat$name,
                     levels = c("E. coli -12% NT","E. coli +12% NT",
                                "DMSO +12% NT","EGCG +12% NT",
                                "Field Wine +12% NT","Field Lees +12% NT",
                                "Lab Wine +12% NT", "Lab Lees +12% NT"))
newdat

df = newdat %>% group_by(Treatment, X, Y, Toxin) %>% summarise(mean = mean(lifetime), stdv = std.error(lifetime), count = n())
df$upper = df$mean + df$stdv
df$lower = df$mean - df$stdv


df$Treatment= factor(df$Treatment, levels=c("E. coli", "DMSO", "EGCG",
                                            "Field Wine", "Lab Wine", "Field Lees",
                                            "Lab Lees"),
                     labels=c("E. coli", "DMSO", "EGCG",
                              "Field Wine", "Lab Wine", "Field Lees",
                              "Lab Lees"))


color_scheme = c("#f25c78",  "#84d94e","#ffa600", "#b64cdf")
df$name = paste(df$Treatment, df$Toxin)
df$name= factor(df$name, levels=c("E. coli -12% NT",
                                       "E. coli +12% NT",
                                       "DMSO +12% NT",
                                       "EGCG +12% NT",
                                      "Field Wine +12% NT", "Lab Wine +12% NT",
                                      "Field Lees +12% NT","Lab Lees +12% NT"))

tox.res = subset(newdat, Treatment == "E. coli")
res2 = t.test(lifetime ~ Toxin, tox.res)
res2

## Normalized
df$per.rel = df$stdv/df$mean
df$norm.mean = df$mean/6.37
df$norm.low = df$norm.mean*(1-df$per.rel)
df$norm.upp = df$norm.mean*(1+df$per.rel)

df = subset(df, name != "E. coli -12% NT")
n = ggplot(df, aes(name,norm.mean,fill=X)) + 
  theme_classic() + labs(x = '', y = 'Relative lifespan of <i>C. elegans</i> with GlyP') + 
  theme(axis.title.y = element_text(size = 14), axis.text.x=element_text(size=13), axis.text.y=element_text(size=12)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=norm.low, ymax=norm.upp), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=color_scheme) +
  theme(legend.position="none") + 
  theme(strip.placement = "outside") +
  theme(axis.title.x = element_blank(),
        axis.title.y =element_markdown(size=14),
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_markdown(size=12),
        axis.text.y = element_text(size=11),
        legend.position="none",
        strip.text.x = element_text(size=13, face = "bold"),
        strip.background = element_blank()) +
  scale_x_discrete(labels=c("Untreated",
                            "200 uM <br> EGCG",
                            "Traditional<br> Wine", "Lab Wine", "Traditional<br>Lees",
                            "Lab Lees")) +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  coord_cartesian(xlim=c(1,6), ylim=c(0,1), clip="off") 

n = n + annotate("text", x = 6, y = 0.75, label = "**", size = 7, color = "#960101") 

plot(n)
ggsave(plot=n, filename="Final-survival-barplot.png", width=6.5, height=4, dpi=600)



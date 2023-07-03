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



to_survival_from_table = function(dat) {
  newdat = data.frame(dat[rep(seq_len(dim(dat)[1]), dat$Dead), c('Time','Dead','Treatment', 'Trial', 'X', 'Y', 'Toxin'), drop = FALSE], row.names=NULL)
  newdat$Treatment= factor(newdat$Treatment, levels=c("E. coli", "DMSO", "EGCG",
                                               "Field Wine", "Lab Wine", "Field Lees",
                                               "Lab Lees"),
                        labels=c("E. coli", "DMSO", "EGCG",
                                 "Field Wine", "Lab Wine", "Field Lees",
                                 "Lab Lees"))
  newdat$X = factor(newdat$X, levels=c("E. coli", "DMSO", "EGCG","Wine", "Lees"),
                           labels=c("OP50 <br><i>E. coli</i>", "0.1% <br>DMSO", "200 uM <br> EGCG", "Wine", "Lees"))
  newdat$Y = factor(newdat$Y, levels=c("Control", "Wine", "Lees"),
                  labels=c("Control", "Wine", "Lees"))
  newdat$Toxin = factor(newdat$Toxin, levels = c("without toxin","with toxin"),
                        labels = c("-12% NT", "+12% NT"))
  newdat$status = 1
  return(newdat)
}

# to_survival_from_table_2 = function(dat) {
#   newdat = data.frame(dat[rep(seq_len(dim(dat)[1]), dat$Dead), c('Time','Dead','Concentration','Treatment', 'Trial'), drop = FALSE], row.names=NULL)
#   newdat$Treatment= factor(newdat$Treatment, levels=c("E. coli", "DMSO", "EGCG",
#                                                       "Field Wine", "Lab Wine", "Field Lees",
#                                                       "Lab Lees"), 
#                            labels=c("E. coli", "DMSO", "EGCG",
#                                     "Field Wine", "Lab Wine", "Field Lees",
#                                     "Lab Lees"))
#   # newdat$Concentration= factor(newdat$Concentration, levels=c("Control", "Low", "Mid", "High"), 
#   #                              labels=c("Control", "Low", "Mid", "High"))
#   newdat$status = 1
#   return(newdat)
# }

newdat = to_survival_from_table(df.csv)


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
newdat$lifetime = newdat$Time -1
df = newdat %>% group_by(Treatment, X, Y, Toxin) %>% summarise(mean = mean(lifetime), stdv = std.error(lifetime), count = n())
df$upper = df$mean + df$stdv
df$lower = df$mean - df$stdv


df$Treatment= factor(df$Treatment, levels=c("E. coli", "DMSO", "EGCG",
                                            "Field Wine", "Lab Wine", "Field Lees",
                                            "Lab Lees"),
                     labels=c("E. coli", "DMSO", "EGCG",
                              "Field Wine", "Lab Wine", "Field Lees",
                              "Lab Lees"))
# labels=c("OP50 <br><i>E. coli</i>", "0.1% <br>DMSO", "200 uM <br> EGCG",
#          "Field Wine", "Lab Wine", "Field Lees",
#          "Lab Lees")
## Fit survival function
fit <- survfit(Surv(newdat$Time, newdat$status) ~ name, data = newdat)
fit

## Graph the K-M plot.
# color_scheme = c("#f29c50",  "#53518c", "#fa5cf2","#c7522a")
g = ggsurvplot(fit, 
               # risk.table = TRUE, # Add risk table
               # risk.table.col = "group", # Change risk table color by groups
               # linetype = "strata", # Change line type by groups
               # ggtheme = theme_classic(), # Change ggplot2 theme
               # palette = color_scheme,
               xlab="Days",
               size = 0.8,
               legend = c(0.9, 0.7),
               # linetype = c(1,2),
               linetype = "strata", # Set linetype to "strata" to use different line colors for different groups
               # palette = c("#fa5cf2", "#fa5cf2", "#ffa600", "#84d94e","#f25c78","#f25c78", "#b64cdf","#b64cdf"),
               # palette = "jco",
               # legend.title = element_blank(),
               legend.labs = c("OP50 <i>E. coli </i>", "OP50 <i>E. coli</i>",
                               "200 uM EGCG",
                               "Field Wine", "Lab Wine", "Field Lees",
                               "Lab Lees"),
               xlim = c(0,11),
               font.legend = c(13, "plain", "black"),
               font.x = c(15, 'plain', 'black'),
               font.y = c(15, 'plain', 'black'),
               font.tickslab = c(11, 'plain', 'black'),
               surv.scale = 'percent'
) 
g$plot

g$plot = g$plot + 
  theme(legend.text = element_markdown(),
        legend.title = element_blank(),
                        legend.background = element_blank(),
                        legend.key.width = unit(1.8,"lines"),
                        legend.key.height = unit(1.4,"lines"),
                        legend.key = element_blank()) +
  scale_linetype_manual(values = c("dotdash", "solid","solid","solid",
                                   "dashed","solid","dashed")) +
  # "#ffa600",  -> yellow
  # "#f25c78","#f25c78" --? redish pinks
  scale_colour_manual(values = c("#f25c78", "#f25c78", "#84d94e","#ffa600","#ffa600",
                                  "#b64cdf","#b64cdf")) +
  annotate("segment", x = 8.9, xend = 8.9, y = 0.42, yend = 0.91) + 
  annotate("text", x = 8.05,  y = 0.66, label = "+ 12% NT", size=5, fontface=2) + 
  annotate("segment", x = 8.9, xend = 8.9, y = 0.93, yend = 1) + 
  annotate("text", x = 8.05,  y = 0.975, label = "\u2013 12% NT", size=5, fontface=2) 
g$plot

ggsave(plot=g$plot, filename="KM-plot.png", width=8, height=4, dpi=600)


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


p = ggplot(df, aes(name,mean,fill=X)) + 
  theme_classic() + labs(x = '', y = 'Mean Lifespan (Days)') + 
  theme(axis.title.y = element_text(size = 14), axis.text.x=element_text(size=13), axis.text.y=element_text(size=12)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(.9)) +
  # facet_grid('. ~ Toxin', scales="free", switch = "x", space='free') +
  scale_fill_manual(values=color_scheme) +
  theme(legend.position="none") + 
  theme(strip.placement = "outside") +
  theme(axis.title.x = element_blank(),
        axis.title.y =element_text(size=15),
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_markdown(size=12),
        axis.text.y = element_text(size=11),
        legend.position="none",
        strip.text.x = element_text(size=13, face = "bold"),
        strip.background = element_blank()) +
  scale_x_discrete(labels=c("OP50 <br><i>E. coli</i>", "OP50 <br><i>E. coli</i>",
                            "200 uM <br> EGCG",
                            "Field Wine", "Lab Wine", "Field Lees",
                            "Lab Lees")) +
  theme(plot.margin = unit(c(0,0,2,0), "lines")) +
  coord_cartesian(xlim=c(1,7), ylim=c(0,8), clip="off") 

tox.vs.notox <- tibble(
  x = c(1, 1, 2,2),
  y = c(7, 7.6, 7.6, 3))

p = p + annotate("text", x = 7, y = 5.0, label = "**", size = 7, color = "#960101") +
  geom_line(data = tox.vs.notox, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 1.5, y = 7.8, label = "***", size = 7, color = "black") 
  # annotate("segment", x = 2, xend = 4.4, y = -1.8, yend = -1.8) +
  # annotate("segment", x = 5.5, xend = 8, y = -1.8, yend = -1.8) +
  # annotate("text", x = 5, y = -1.8, label = "+ 12% NT", fontface =2, size=5)  +
  # annotate("text", x = 1, y = -1.8, label = "- 12% NT", fontface =2, size=5)
plot(p)

ggsave(plot=p, filename="Bar-plot.png", width=8, height=4, dpi=600)


png(filename = file.path(getwd(),"/Combined-A-B.png"), width = 8, height = 7, units = "in", res = 600)
kmbar = ggarrange(plotlist = list(g$plot, p), ncol = 1, nrow = 2, labels = c("A", "B"))
plot(kmbar)

dev.off()

print(paste("Saved at ", getwd()))



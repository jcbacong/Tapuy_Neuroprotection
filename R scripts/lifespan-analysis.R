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

## Set working directory to the file location using Rstudio API
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Load survival csv file
data = read.csv("survival-data.csv", header = T)

## Define a function that will create a survival/event table
## Per row will be a biological sample that either survived (0) or dead (1).
## No censored will be applied here. 

to_survival_from_table = function(dat) {
  newdat = data.frame(dat[rep(seq_len(dim(dat)[1]), dat$dead), c('day','dead','group', 'replicate'), drop = FALSE], row.names=NULL)
  newdat$group = factor(newdat$group, levels=c("F", "DF", "FT", "DFT"), labels=c("Untreated", "0.1% DMSO", "12% NT", "0.1% DMSO + 12% NT"))
  newdat$status = 1
  return(newdat)
}

newdat = to_survival_from_table(data)
newdat

## Fit survival function
fit <- survfit(Surv(newdat$day, newdat$status) ~ group, data = newdat)
fit

## Graph the K-M plot.
color_scheme = c("#f29c50",  "#53518c", "#fa5cf2","#c7522a")
g = ggsurvplot(fit, 
           # risk.table = TRUE, # Add risk table
           # risk.table.col = "group", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           ggtheme = theme_classic(), # Change ggplot2 theme
           palette = color_scheme,
           xlab="Time in days",
           size = 0.8,
           legend = c(0.24, 0.25),
           legend.title = element_blank(),
           legend.labs = c("Untreated",  "0.1% DMSO", "12% NT", "0.1% DMSO + 12% NT"),
           xlim = c(0,11),
           font.legend = c(11, "plain", "black"),
           font.x = c(15, 'plain', 'black'),
           font.y = c(15, 'plain', 'black'),
           font.tickslab = c(11, 'plain', 'black'),
           surv.scale = 'percent'
          ) + guides(colour = guide_legend(nrow = 4)) 
g


## Plot mean lifetime
## Prepare dataframe with mean and standard dev
newdat
newdat$lifetime = newdat$day -1
df = newdat %>% group_by(group) %>% summarise(mean = mean(lifetime), stdv = sd(lifetime), count = n())
df$upper = df$mean + df$stdv
df$lower = df$mean - df$stdv
# df$group = as.factor(df$group, levels=c("Untreated", "DMSO", "Toxin", "DMSO & Toxin" ), labels=c("Untreated", "DMSO", "Toxin", "DMSO & Toxin" ) )


## Visualize plot
## Using bar plots of means. Test for significance will be t-test
## 1. Checked anova. p-val = 1.34e-6 
anova.res = oneway.test(lifetime ~ group, data = newdat,var.equal = FALSE) # assuming equal variances
anova.res

## 2. Used pairwise independent two-tailed t-test with hochberg correction
stat.test = newdat %>% t_test(lifetime ~ group, p.adjust.method = "hochberg", var.equal = F)
stat.test

## 1. One way to visualize, automatic
# bp <- ggbarplot(
#   newdat, x = "group", y = "lifetime", add = "mean_sd", 
#   color= "group", palette = color_scheme,
#   position = position_dodge(0.8)
# )
# bp
# stat.test <- stat.test %>%
#   add_xy_position(fun = "mean_sd", x = "group", dodge = 0.8) 
# bp + stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01)
# bp


## 1. Manual addition of significance bars (preferred)
p = ggplot(df, aes(group,mean,fill=group)) + 
    theme_classic() + labs(x = '', y = 'Mean Lifespan') + 
    theme(axis.title.y = element_text(size = 14), axis.text.x=element_text(size=13), axis.text.y=element_text(size=12)) +
    geom_bar(stat="identity", color="black",position=position_dodge()) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(.9)) +
    scale_fill_manual(values=color_scheme) +
    theme(legend.position="none") + ylim(c(0,12)) +
    scale_x_discrete(labels = c("Untreated", "0.1% DMSO", "12% NT", "0.1% DMSO\n+ 12% NT")) 
p

## Significance bars
## Only with ***
untreated.vs.dmso.p <- tibble(
  x = c(1, 1, 1.98, 1.98),
  y = c(9.5, 10, 10, 9.5)
) ## n.s

untreated.vs.toxin.p <- tibble(
  x = c(1, 1, 3, 3),
  y = c(10.2, 10.7, 10.7, 10.2)
) ## p = 0.000081

dmso.vs.toxin.p <- tibble(
  x = c(2.02, 2.02, 3, 3),
  y = c(9.5, 10, 10, 9.5)
) ## p = 0.0000259

untreated.vs.dmso.toxin.p <- tibble(
  x = c(1, 1, 4, 4),
  y = c(10.8, 11.3, 11.3, 10.8)
) ## p = 0.0000259

toxin.vs.dmso.toxin.p <- tibble(
  x = c(3, 3, 4, 4),
  y = c(8, 8.3, 8.3, 8)
) ## p = 0.035

p = p + 
    geom_line(data = untreated.vs.dmso.p, aes(x= x, y = y, group=1), inherit.aes = F) +
    annotate("text", x = 1.5, y = 10.4, label = "n.s", size = 4, color = "#22292F") +
    geom_line(data = untreated.vs.toxin.p, aes(x= x, y = y, group=1), inherit.aes = F) +
    annotate("text", x = 2, y = 10.8, label = "***", size = 5, color = "#22292F") +
    geom_line(data = dmso.vs.toxin.p, aes(x= x, y = y, group=1), inherit.aes = F) +
    annotate("text", x = 2.5, y = 10.1, label = "***", size = 5, color = "#22292F") +
    geom_line(data = untreated.vs.dmso.toxin.p, aes(x= x, y = y, group=1), inherit.aes = F) +
    annotate("text", x = 2.5, y = 11.4, label = "*", size = 5, color = "#22292F") +
    geom_line(data = toxin.vs.dmso.toxin.p, aes(x= x, y = y, group=1), inherit.aes = F) +
    annotate("text", x = 3.5, y = 8.4, label = "*", size = 5, color = "#22292F") 
p

## 3. Using boxplots for non-parametric consideration. Mann whitney u test will be used.
# p = ggplot(newdat, aes(x = group, y = lifetime, fill =group)) +                 # Draw ggplot2 boxplot
#   geom_boxplot() +
#   theme_classic() + labs(x = '', y = 'Mean Lifespan') + theme(axis.title.y = element_text(size = 14), axis.text.x=element_text(size=13)) +
#   stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
#   scale_fill_manual(values=color_scheme) +
#   theme(legend.position="none") +
#   stat_compare_means(comparisons = list(c("Untreated", "DMSO"),
#                                  c("Untreated", "Toxin")),
#               label.y = c(11, 12),
#               method = 'wilcox.test',
#               map_signif_level=FALSE,
#               symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
# p

tiff(filename = file.path(getwd(),"/KM-Bar-Plot.tiff"), width = 6, height = 7, units = "in", res = 600)
png(filename = file.path(getwd(),"/KM-Bar-Plot.png"), width = 6, height = 7, units = "in", res = 600)
kmbar = ggarrange(plotlist = list(g$plot, p), ncol = 1, nrow = 2)
                  # labels = c("A", "B"))
plot(kmbar)


# png(filename = file.path(getwd(),"/A-plot.png"), res = 600)
# plot(g$plot)
# 

# png(filename = file.path(getwd(),"/B-plot.png"), width = 6, height = 3.5, units = "in", res = 600)
# plot(p)
dev.off()

print(paste("Saved at ", getwd()))





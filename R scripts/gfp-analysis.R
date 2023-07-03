  
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
  if( toupper(substr(x,start=1,stop=2)) == "FW"){
    return("Field Wine")
  } else if ( toupper(substr(x,start=1,stop=2)) == "FL") {
    return("Field Lees")
  } else if( toupper(substr(x,start=1,stop=2)) == "LW") {
    return("Lab Wine")
  } else if( toupper(substr(x,start=1,stop=2)) == "LL") {
    return("Lab Lees")
  } else {
    return(x)
  }
}

## Prepare group name
rename_fill = function(x) {
  if( toupper(substr(x,start=1,stop=2)) == "FW"){
    return("Wine")
  } else if ( toupper(substr(x,start=1,stop=2)) == "FL") {
    return("Lees")
  } else if( toupper(substr(x,start=1,stop=2)) == "LW") {
    return("Wine")
  } else if( toupper(substr(x,start=1,stop=2)) == "LL") {
    return("Lees")
  } else {
    return(x)
  }
}


rename_concentration = function(x) {
  if( toupper(substr(x,start=3,stop=3)) == "L"){
    return("Low")
  } else if ( toupper(substr(x,start=3,stop=3)) == "M") {
    return("Mid")
  } else if( toupper(substr(x,start=3,stop=3)) == "H") {
    return("High")
  } else {
    return("Control")
  }
}

strain_summary_normalize = function(strain.summary, rel.time = "t0") {
  strain.summary$rel.error = strain.summary$sem.rawint/strain.summary$mean.rawint
  setDT(strain.summary)
  set(strain.summary, j = "mean.rawint", value = as.numeric(strain.summary[["mean.rawint"]]))
  set(strain.summary, j = "rel.error", value = as.numeric(strain.summary[["rel.error"]]))
  set(strain.summary, j = "obs.rawint", value = as.numeric(strain.summary[["obs.rawint"]]))
  
  strain.summary.norm = strain.summary[,`:=`(mean.rawint = mean.rawint / mean.rawint[time == rel.time],
                                             rel.error = rel.error + rel.error[time == rel.time],
                                             obs.rawint = obs.rawint + obs.rawint[time == rel.time]),
                                       by = control]
  return(strain.summary.norm)
}

## Within griup mean difference
## Use the following to simulate a gaussian curve
t_test_means = function(df.norm.6h, X, Y="T", alternative="less") {
  x.mean = df.norm.6h[df.norm.6h$control==X,mean.rawint]
  x.sem = df.norm.6h[df.norm.6h$control==X,mean.rawint]*df.norm.6h[df.norm.6h$control==X,rel.error]
  x.obs = df.norm.6h[df.norm.6h$control==X,obs.rawint]
  
  y.mean = df.norm.6h[df.norm.6h$control=="T",mean.rawint]
  y.sem = df.norm.6h[df.norm.6h$control=="T",mean.rawint]*df.norm.6h[df.norm.6h$control=="T",rel.error]
  y.obs = df.norm.6h[df.norm.6h$control==X,obs.rawint]
  
  x = tsum.test(x.mean, x.sem, x.obs,
                y.mean, y.sem, y.obs,
                alternative = alternative)
  return(x)
}

stat_results_df = function(df.t, Y="T", alternative="less") {
  results_df <- data.frame()
  for (x in df.t$control) {
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

######################### READ FILES #############################
### OH441
strain = "OH441"
controls.csv = read.csv(file.path(getwd(), strain,"controls-2.csv"), header = T)
treatment.csv = read.csv(file.path(getwd(), strain,"treatment-2.csv"), header = T)
treatment.csv$control = lapply(treatment.csv$control, toupper)
df = rbind(controls.csv, treatment.csv)

df$concentration = lapply(df$control, rename_concentration)
df$group = lapply(df$control, rename_wine_or_lees)
df$name = lapply(df$control, rename_fill)
df[df == "OP50"] = "T"
df[df == "Neg"] = "U"

## Factor categories
df$control <- factor(df$control, levels = c('U', 'T', 'DMSO', 'EGCG', 'Metab',
                                                                'FWL', 'FWM', 'FWH',
                                                                'LWL', 'LWM', 'LWH',
                                                                'FLL', 'FLM', 'FLH',
                                                                'LLL', 'LLM', 'LLH'))

df$concentration <- factor(df$concentration, levels = c('Control','Low', 'Mid','High'))

df$group <- factor(df$group, levels = c('T', "DMSO", "EGCG",'Field Wine', 'Lab Wine',
                                                        'Field Lees', 'Lab Lees'))
df$name = factor(df$name, levels = c('T', "DMSO", "EGCG",'Wine', 'Lees'))
df$neuron = strain 
df = subset(df, concentration %in% c("Control", "High"))
df = subset(df, !(control %in% c("Metab", "DMSO")))

df.summary <- df %>% 
                group_by(control, time, group, neuron, concentration, name) %>%
                summarise(mean.rawint = mean(RawIntDen),
                          sem.rawint = std.error(RawIntDen),
                          obs.rawint = n()) %>%
                as.data.frame()



## Statistics
df.norm.6h = subset(strain_summary_normalize(df.summary, rel.time = "t6"),time == "t0")
res1 = stat_results_df(df.norm.6h)
res1


### Plot indivudals
color_scheme = c("#f25c78",  "#84d94e","#ffa600", "#b64cdf")

y1 = ggplot(df.norm.6h, aes(group, mean.rawint, fill=name)) +
  theme_classic() + labs(x = '', y = 'Neurodegeneration after 6h') +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin= mean.rawint, 
                    ymax= mean.rawint * (1+rel.error)),
                width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>", "200 uM <br> EGCG", 
                              "Field\n Wine","Lab\n Wine", 
                              "Field\n Lees", "Lab\nLees")) +
  scale_fill_manual(values=color_scheme) +
  geom_hline(aes(yintercept = df.norm.6h$mean.rawint[1], color = "T"), linetype='dashed', color='#960101') +
  theme(axis.text.x = element_markdown(color="black", size=12),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=12),
        axis.title.y = element_text(size=15),
        plot.margin=unit(c(0,0,0,0), 'cm'),
        legend.position="none"
        
  ) + 
  scale_y_continuous(limits = c(0,4)) 
# geom_hline(aes(yintercept = toxic.recovery.mean, color = "U"), linetype='dashed', color='#960101') +
# scale_color_manual(values="#960101")
y1

h1 = y1 + 
  annotate("text", x = 0.73, y = df.norm.6h$mean.rawint[1] + 0.2, label = paste(c("(",round(df.norm.6h$mean.rawint[1],2),")"), collapse=""), size = 3, color = "#960101") +
  annotate("text", x=4, y = 1.6, label = "***", size = 6, color = "#960101") + 
  annotate("text", x=5, y = 2.5, label = "*", size = 6, color = "#960101") +
  annotate("text", x=6, y = 1.2, label = "***", size = 6, color = "#960101")
h1



df.norm.9h = subset(strain_summary_normalize(df.summary, rel.time = "t9"),time == "t6")
res2 = stat_results_df(df.norm.9h)
res2

y2 = ggplot(df.norm.9h, aes(group, mean.rawint, fill=name)) +
  theme_classic() + labs(x = '', y = 'Neurodegeneration after 9h') +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin= mean.rawint, 
                    ymax= mean.rawint * (1+rel.error)),
                width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>", "200 uM <br> EGCG", 
                              "Field\n Wine","Lab\n Wine", 
                              "Field\n Lees", "Lab\nLees")) +
  scale_fill_manual(values=color_scheme) +
  geom_hline(aes(yintercept = df.norm.9h$mean.rawint[1], color = "T"), linetype='dashed', color='#960101') +
  theme(axis.text.x = element_markdown(color="black", size=12),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=12),
        axis.title.y = element_text(size=15),
        plot.margin=unit(c(0,0,0,0), 'cm'),
        legend.position="none"
        
  ) +
  scale_y_continuous(limits=c(0,4))
# geom_hline(aes(yintercept = toxic.recovery.mean, color = "U"), linetype='dashed', color='#960101') +
# scale_color_manual(values="#960101")
y2

h2 = y2 + 
  annotate("text", x = 0.73, y = df.norm.9h$mean.rawint[1] + 0.2, label = paste(c("(",round(df.norm.9h$mean.rawint[1],2),")"), collapse=""), size = 3, color = "#960101")
h2

ggsave(plot=h2, filename="BCH205-OH441-9h.png", width=5.8, height=3.5, dpi=600)
ggsave(plot=h1, filename="BCH205-OH441-6h.png", width=5.8, height=3.5, dpi=600)
print(paste("Saved at ", getwd()))



# ### LX929
strain = "LX929"
controls.csv = read.csv(file.path(getwd(), strain,"controls-2.csv"), header = T)
treatment.csv = read.csv(file.path(getwd(), strain,"treatment-2.csv"), header = T)
treatment.csv$control = lapply(treatment.csv$control, toupper)
df = rbind(controls.csv, treatment.csv)

df$concentration = lapply(df$control, rename_concentration)
df$group = lapply(df$control, rename_wine_or_lees)
df$name = lapply(df$control, rename_fill)
df[df == "OP50"] = "T"
df[df == "Neg"] = "U"

## Factor categories
df$control <- factor(df$control, levels = c('U', 'T', 'DMSO', 'EGCG', 'Metab',
                                            'FWL', 'FWM', 'FWH',
                                            'LWL', 'LWM', 'LWH',
                                            'FLL', 'FLM', 'FLH',
                                            'LLL', 'LLM', 'LLH'))

df$concentration <- factor(df$concentration, levels = c('Control','Low', 'Mid','High'))

df$group <- factor(df$group, levels = c('T', "DMSO", "EGCG",'Field Wine', 'Lab Wine',
                                        'Field Lees', 'Lab Lees'))
df$name = factor(df$name, levels = c('T', "DMSO", "EGCG",'Wine', 'Lees'))
df$neuron = strain 
df = subset(df, concentration %in% c("Control", "High"))
df = subset(df, !(control %in% c("Metab", "DMSO")))

df.summary <- df %>% 
  group_by(control, time, group, neuron, concentration, name) %>%
  summarise(mean.rawint = mean(RawIntDen),
            sem.rawint = std.error(RawIntDen),
            obs.rawint = n()) %>%
  as.data.frame()



## Statistics
df.norm.6h = subset(strain_summary_normalize(df.summary, rel.time = "t6"),time == "t0")
res1 = stat_results_df(df.norm.6h)
res1


### Plot indivudals
color_scheme = c("#f25c78",  "#84d94e","#ffa600", "#b64cdf")

y1 = ggplot(df.norm.6h, aes(group, mean.rawint, fill=name)) +
  theme_classic() + labs(x = '', y = 'Neurodegeneration after 6h') +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin= mean.rawint, 
                    ymax= mean.rawint * (1+rel.error)),
                width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>", "200 uM <br> EGCG", 
                              "Field\n Wine","Lab\n Wine", 
                              "Field\n Lees", "Lab\nLees")) +
  scale_fill_manual(values=color_scheme) +
  geom_hline(aes(yintercept = df.norm.6h$mean.rawint[1], color = "T"), linetype='dashed', color='#960101') +
  theme(axis.text.x = element_markdown(color="black", size=12),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=12),
        axis.title.y = element_text(size=15),
        plot.margin=unit(c(0,0,0,0), 'cm'),
        legend.position="none"
        
  ) + 
  scale_y_continuous(limits = c(0,2)) 
# geom_hline(aes(yintercept = toxic.recovery.mean, color = "U"), linetype='dashed', color='#960101') +
# scale_color_manual(values="#960101")
y1

h1 = y1 + 
  annotate("text", x = 0.73, y = df.norm.6h$mean.rawint[1] + 0.07, label = paste(c("(",round(df.norm.6h$mean.rawint[1],2),")"), collapse=""), size = 3, color = "#960101") 
h1

df.norm.9h = subset(strain_summary_normalize(df.summary, rel.time = "t9"),time == "t6")
res2 = stat_results_df(df.norm.9h)
res2

y2 = ggplot(df.norm.9h, aes(group, mean.rawint, fill=name)) +
  theme_classic() + labs(x = '', y = 'Neurodegeneration after 9h') +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin= mean.rawint, 
                    ymax= mean.rawint * (1+rel.error)),
                width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>", "200 uM <br> EGCG", 
                              "Field\n Wine","Lab\n Wine", 
                              "Field\n Lees", "Lab\nLees")) +
  scale_fill_manual(values=color_scheme) +
  geom_hline(aes(yintercept = df.norm.9h$mean.rawint[1], color = "T"), linetype='dashed', color='#960101') +
  theme(axis.text.x = element_markdown(color="black", size=12),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=12),
        axis.title.y = element_text(size=15),
        plot.margin=unit(c(0,0,0,0), 'cm'),
        legend.position="none"
        
  ) +
  scale_y_continuous(limits=c(0,2))
# geom_hline(aes(yintercept = toxic.recovery.mean, color = "U"), linetype='dashed', color='#960101') +
# scale_color_manual(values="#960101")
y2

h2 = y2 + 
  annotate("text", x = 0.73, y = df.norm.9h$mean.rawint[1] + 0.07, label = paste(c("(",round(df.norm.9h$mean.rawint[1],2),")"), collapse=""), size = 3, color = "#960101") + 
  annotate("text", x=6, y = 0.83, label="**", size=6, color="#960101")
h2

ggsave(plot=h2, filename="BCH205-LX929-9h.png", width=5.8, height=3.5, dpi=600)
ggsave(plot=h1, filename="BCH205-LX929-6h.png", width=5.8, height=3.5, dpi=600)
print(paste("Saved at ", getwd()))


### BZ555
strain = "BZ555"
controls.csv = read.csv(file.path(getwd(), strain,"controls-2.csv"), header = T)
treatment.csv = read.csv(file.path(getwd(), strain,"treatment-2.csv"), header = T)
treatment.csv$control = lapply(treatment.csv$control, toupper)
df = rbind(controls.csv, treatment.csv)

df$concentration = lapply(df$control, rename_concentration)
df$group = lapply(df$control, rename_wine_or_lees)
df$name = lapply(df$control, rename_fill)
df[df == "OP50"] = "T"
df[df == "Neg"] = "U"

## Factor categories
df$control <- factor(df$control, levels = c('U', 'T', 'DMSO', 'EGCG', 'Metab',
                                            'FWL', 'FWM', 'FWH',
                                            'LWL', 'LWM', 'LWH',
                                            'FLL', 'FLM', 'FLH',
                                            'LLL', 'LLM', 'LLH'))

df$concentration <- factor(df$concentration, levels = c('Control','Low', 'Mid','High'))

df$group <- factor(df$group, levels = c('T', "DMSO", "EGCG",'Field Wine', 'Lab Wine',
                                        'Field Lees', 'Lab Lees'))
df$name = factor(df$name, levels = c('T', "DMSO", "EGCG",'Wine', 'Lees'))
df$neuron = strain 
df = subset(df, concentration %in% c("Control", "High"))
df = subset(df, !(control %in% c("Metab", "DMSO")))

df.summary <- df %>% 
  group_by(control, time, group, neuron, concentration, name) %>%
  summarise(mean.rawint = mean(RawIntDen),
            sem.rawint = std.error(RawIntDen),
            obs.rawint = n()) %>%
  as.data.frame()


## Statistics
df.norm.6h = subset(strain_summary_normalize(df.summary, rel.time = "t6"),time == "t0")
res1 = stat_results_df(df.norm.6h)
res1


### Plot indivudals
color_scheme = c("#f25c78",  "#84d94e","#ffa600", "#b64cdf")

y1 = ggplot(df.norm.6h, aes(group, mean.rawint, fill=name)) +
  theme_classic() + labs(x = '', y = 'Neurodegeneration after 6h') +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin= mean.rawint, 
                    ymax= mean.rawint * (1+rel.error)),
                width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>", "200 uM <br> EGCG", 
                              "Field\n Wine","Lab\n Wine", 
                              "Field\n Lees", "Lab\nLees")) +
  scale_fill_manual(values=color_scheme) +
  geom_hline(aes(yintercept = df.norm.6h$mean.rawint[1], color = "T"), linetype='dashed', color='#960101') +
  theme(axis.text.x = element_markdown(color="black", size=12),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=12),
        axis.title.y = element_text(size=15),
        plot.margin=unit(c(0,0,0,0), 'cm'),
        legend.position="none"
        
  ) + 
  scale_y_continuous(limits = c(0,2)) 
# geom_hline(aes(yintercept = toxic.recovery.mean, color = "U"), linetype='dashed', color='#960101') +
# scale_color_manual(values="#960101")
y1

h1 = y1 + 
  annotate("text", x = 0.73, y = df.norm.6h$mean.rawint[1] + 0.07, label = paste(c("(",round(df.norm.6h$mean.rawint[1],2),")"), collapse=""), size = 3, color = "#960101") +
  annotate("text", x=4, y=1.38, label="***", size=6, color = "#960101" ) +
  annotate("text", x=5, y=0.22, label="***", size=6, color = "#960101" ) +
  annotate("text", x=6, y=0.18, label="***", size=6, color = "#960101" ) 
h1

df.norm.9h = subset(strain_summary_normalize(df.summary, rel.time = "t9"),time == "t6")
res2 = stat_results_df(df.norm.9h)
res2

y2 = ggplot(df.norm.9h, aes(group, mean.rawint, fill=name)) +
  theme_classic() + labs(x = '', y = 'Neurodegeneration after 9h') +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin= mean.rawint, 
                    ymax= mean.rawint * (1+rel.error)),
                width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>", "200 uM <br> EGCG", 
                              "Field\n Wine","Lab\n Wine", 
                              "Field\n Lees", "Lab\nLees")) +
  scale_fill_manual(values=color_scheme) +
  geom_hline(aes(yintercept = df.norm.9h$mean.rawint[1], color = "T"), linetype='dashed', color='#960101') +
  theme(axis.text.x = element_markdown(color="black", size=12),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=12),
        axis.title.y = element_text(size=15),
        plot.margin=unit(c(0,0,0,0), 'cm'),
        legend.position="none"
        
  ) +
  scale_y_continuous(limits=c(0,2))
# geom_hline(aes(yintercept = toxic.recovery.mean, color = "U"), linetype='dashed', color='#960101') +
# scale_color_manual(values="#960101")
y2

h2 = y2 + 
  annotate("text", x = 0.73, y = df.norm.9h$mean.rawint[1] + 0.07, label = paste(c("(",round(df.norm.9h$mean.rawint[1],2),")"), collapse=""), size = 3, color = "#960101") 
h2

ggsave(plot=h2, filename="BCH205-BZ555-9h.png", width=5.8, height=3.5, dpi=600)
ggsave(plot=h1, filename="BCH205-BZ555-6h.png", width=5.8, height=3.5, dpi=600)
print(paste("Saved at ", getwd()))


### EG1285
strain = "EG1285"
controls.csv = read.csv(file.path(getwd(), strain,"controls-2.csv"), header = T)
treatment.csv = read.csv(file.path(getwd(), strain,"treatment-2.csv"), header = T)
treatment.csv$control = lapply(treatment.csv$control, toupper)
df = rbind(controls.csv, treatment.csv)

df$concentration = lapply(df$control, rename_concentration)
df$group = lapply(df$control, rename_wine_or_lees)
df$name = lapply(df$control, rename_fill)
df[df == "OP50"] = "T"
df[df == "Neg"] = "U"

## Factor categories
df$control <- factor(df$control, levels = c('U', 'T', 'DMSO', 'EGCG', 'Metab',
                                            'FWL', 'FWM', 'FWH',
                                            'LWL', 'LWM', 'LWH',
                                            'FLL', 'FLM', 'FLH',
                                            'LLL', 'LLM', 'LLH'))

df$concentration <- factor(df$concentration, levels = c('Control','Low', 'Mid','High'))

df$group <- factor(df$group, levels = c('T', "DMSO", "EGCG",'Field Wine', 'Lab Wine',
                                        'Field Lees', 'Lab Lees'))
df$name = factor(df$name, levels = c('T', "DMSO", "EGCG",'Wine', 'Lees'))
df$neuron = strain 
df = subset(df, concentration %in% c("Control", "High"))
df = subset(df, !(control %in% c("Metab", "DMSO")))

df.summary <- df %>% 
  group_by(control, time, group, neuron, concentration, name) %>%
  summarise(mean.rawint = mean(RawIntDen),
            sem.rawint = std.error(RawIntDen),
            obs.rawint = n()) %>%
  as.data.frame()


## Statistics
df.norm.6h = subset(strain_summary_normalize(df.summary, rel.time = "t6"),time == "t0")
res1 = stat_results_df(df.norm.6h)
res1

### Plot indivudals
color_scheme = c("#f25c78",  "#84d94e","#ffa600", "#b64cdf")

y1 = ggplot(df.norm.6h, aes(group, mean.rawint, fill=name)) +
  theme_classic() + labs(x = '', y = 'Neurodegeneration after 6h') +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin= mean.rawint, 
                    ymax= mean.rawint * (1+rel.error)),
                width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>", "200 uM <br> EGCG", 
                              "Field\n Wine","Lab\n Wine", 
                              "Field\n Lees", "Lab\nLees")) +
  scale_fill_manual(values=color_scheme) +
  geom_hline(aes(yintercept = df.norm.6h$mean.rawint[1], color = "T"), linetype='dashed', color='#960101') +
  theme(axis.text.x = element_markdown(color="black", size=12),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=12),
        axis.title.y = element_text(size=15),
        plot.margin=unit(c(0,0,0,0), 'cm'),
        legend.position="none"
        
  ) + 
  scale_y_continuous(limits = c(0,5)) 
# geom_hline(aes(yintercept = toxic.recovery.mean, color = "U"), linetype='dashed', color='#960101') +
# scale_color_manual(values="#960101")
y1

h1 = y1 + 
  annotate("text", x = 0.73, y = df.norm.6h$mean.rawint[1] + 0.2, label = paste(c("(",round(df.norm.6h$mean.rawint[1],2),")"), collapse=""), size = 3, color = "#960101") +
  annotate("text", x=2, y=1.65, label="***", size=6, color = "#960101" ) +
  annotate("text", x=3, y=2.8, label="***", size=6, color = "#960101" ) +
  annotate("text", x=4, y=2.83, label="***", size=6, color = "#960101" ) +
  annotate("text", x=5, y=2.75, label="***", size=6, color = "#960101" ) +
  annotate("text", x=6, y=0.85, label="***", size=6, color = "#960101" ) 
h1

df.norm.9h = subset(strain_summary_normalize(df.summary, rel.time = "t9"),time == "t6")
res2 = stat_results_df(df.norm.9h)
res2

y2 = ggplot(df.norm.9h, aes(group, mean.rawint, fill=name)) +
  theme_classic() + labs(x = '', y = 'Neurodegeneration after 9h') +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin= mean.rawint, 
                    ymax= mean.rawint * (1+rel.error)),
                width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>", "200 uM <br> EGCG", 
                              "Field\n Wine","Lab\n Wine", 
                              "Field\n Lees", "Lab\nLees")) +
  scale_fill_manual(values=color_scheme) +
  geom_hline(aes(yintercept = df.norm.9h$mean.rawint[1], color = "T"), linetype='dashed', color='#960101') +
  theme(axis.text.x = element_markdown(color="black", size=12),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = 'black', size=10),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=12),
        axis.title.y = element_text(size=15),
        plot.margin=unit(c(0,0,0,0), 'cm'),
        legend.position="none"
        
  ) +
  scale_y_continuous(limits=c(0,5))
# geom_hline(aes(yintercept = toxic.recovery.mean, color = "U"), linetype='dashed', color='#960101') +
# scale_color_manual(values="#960101")
y2

h2 = y2 + 
  annotate("text", x = 0.73, y = df.norm.9h$mean.rawint[1] + 0.18, label = paste(c("(",round(df.norm.9h$mean.rawint[1],2),")"), collapse=""), size = 3, color = "#960101") +
  annotate("text", x=3, y=0.3, label="***", size=6, color = "#960101" ) +
  annotate("text", x=5, y=0.31, label="***", size=6, color = "#960101" ) 
h2

ggsave(plot=h2, filename="BCH205-EG1285-9h.png", width=5.8, height=3.5, dpi=600)
ggsave(plot=h1, filename="BCH205-EG1285-6h.png", width=5.8, height=3.5, dpi=600)
print(paste("Saved at ", getwd()))

##########################################################################################################
################################## NEURONAL TOXICITY  ####################################################
##########################################################################################################

neuronal.toxicity.s = rbind(oh441.summary[oh441.summary$control == "T",],
                          bz555.summary[bz555.summary$control == "T",],
                          eg1285.summary[eg1285.summary$control == "T",],
                          lx929.summary[lx929.summary$control == "T",]
                          )
neuronal.toxicity.s$rel.error = neuronal.toxicity.s$sem.rawint/neuronal.toxicity.s$mean.rawint

setDT(neuronal.toxicity.s)
set(neuronal.toxicity.s, j = "mean.rawint", value = as.numeric(neuronal.toxicity.s[["mean.rawint"]]))
set(neuronal.toxicity.s, j = "rel.error", value = as.numeric(neuronal.toxicity.s[["rel.error"]]))

## Metric for toxicity is t0/t6. Higher value mens higher toxicity.
neuronal.toxicity.s.summary = neuronal.toxicity.s[,`:=`(mean.rawint = mean.rawint/mean.rawint[time == "t6"],
                                                         rel.error = rel.error + rel.error[time == "t6"],
                                                        obs.rawint = obs.rawint + obs.rawint[time == "t6"]),
                                                  by = neuron]

neuronal.toxicity.s.summary = subset(neuronal.toxicity.s.summary, time == "t0", select=-c(time, control, concentration, group))
neuronal.toxicity.s.summary$upper.bound = neuronal.toxicity.s.summary$mean.rawint * (1+neuronal.toxicity.s.summary$rel.error)
neuronal.toxicity.s.summary$lower.bound = neuronal.toxicity.s.summary$mean.rawint * (1-neuronal.toxicity.s.summary$rel.error)
neuronal.toxicity.s.summary$sem = neuronal.toxicity.s.summary$mean.rawint * neuronal.toxicity.s.summary$rel.error

## Statistics for neurodegeneration
## Using bar plots of means. Test for significance using one sample t-test using the
## the mean and standard error of the 
## ratio of GFP expression of pre- to post-6h toxin.
## Null hypothesis: Sample mean is less than 1
## Alternative hypothesis: Sample mean is greater than 1

oh441.norm = c(rtruncnorm(10e6,
                          a=0, b=100,  
                          mean = neuronal.toxicity.s.summary[neuronal.toxicity.s.summary$neuron=="OH441",mean.rawint],
                          sd = neuronal.toxicity.s.summary[neuronal.toxicity.s.summary$neuron=="OH441",sem]))
t.test(oh441.norm, mu = 1, alternative = "greater") ## p = 2.2e-16. Reject null, OH441 has greater mean compared to 1

lx929.norm = c(rtruncnorm(10e6,
                          a=0, b=100, 
                          mean = neuronal.toxicity.s.summary[neuronal.toxicity.s.summary$neuron=="LX929",mean.rawint],
                          sd = neuronal.toxicity.s.summary[neuronal.toxicity.s.summary$neuron=="LX929",sem]))
t.test(lx929.norm, mu = 1, alternative = "greater") ## p = 1. Reject null. LX929 has lower mean.

bz555.norm = c(rtruncnorm(10e6,
                          a=0, b=100, 
                          mean = neuronal.toxicity.s.summary[neuronal.toxicity.s.summary$neuron=="BZ555",mean.rawint],
                          sd = neuronal.toxicity.s.summary[neuronal.toxicity.s.summary$neuron=="BZ555",sem]))
t.test(bz555.norm, mu = 1, alternative = "greater") ## p = 2.2e-16. Reject null, BZ555 has greater mean compared to 1

eg1285.norm = c(rtruncnorm(10e6,
                           a=0, b=100, 
                           mean = neuronal.toxicity.s.summary[neuronal.toxicity.s.summary$neuron=="EG1285",mean.rawint],
                           sd = neuronal.toxicity.s.summary[neuronal.toxicity.s.summary$neuron=="EG1285",sem]))
t.test(eg1285.norm, mu = 1, alternative = "greater") ## p = 2.2e-16. Reject null. EG1285 has higher mean compared to 1


## Visualize
## 

neuronal.toxicity.s.summary$neuron = factor(neuronal.toxicity.s.summary$neuron,
                                     level = c( "OH441", "LX929" ,"BZ555" , "EG1285"),
                                     labels = c("Pan-neuronal", "Cholinergic", "Dopaminergic", "GABAergic"))

color_scheme = c("#9b59b6",  "#f2d427","#3498db", "#1abc9c")
x1 = ggplot(neuronal.toxicity.s.summary, aes(neuron, mean.rawint, fill=neuron)) +
  theme_classic() + labs(x = '', y = 'Neurodegeneration after 6h') +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=lower.bound, ymax=upper.bound),
                width=.2,
                position=position_dodge(.9)) +
  # scale_x_discrete(labels = c("Pan-neuronal", "Cholinergic", "Dopaminergic", "GABA-ergic")) +
  ylim(0,5) +
  scale_fill_manual(values=color_scheme) +
  # geom_hline(aes(yintercept = healthy.baseline.mean, color = "U"), linetype='dotted', color='red') +
  # scale_color_manual(values="red")+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(color = c("black", "#960101", "black", "black", "black"), size=12),
        axis.text.x = element_text(size=13),
        axis.title.y = element_text(size=15)
  ) +
  
  geom_hline(aes(yintercept = 1.0, color = "U"), linetype='dashed', color='#960101') +
  annotate("text", x = 1.0, y = 2.85, label = "***", size = 6, color = "#960101") +
  annotate("text", x = 3.0, y = 1.95, label = "***", size = 6, color = "#960101") +
  annotate("text", x = 4.0, y = 4.5, label = "***", size = 6, color = "#960101")
x1

tiff(filename = file.path(getwd(),"/Neuronal-Toxicity.tiff"), width = 5.5, height = 3.5, units = "in", res = 600)
plot(x1)
dev.off()
print(paste("Saved at ", getwd()))



##### SKIP THIS LINE ####
# 
# ##########################################################################################################
# ##################################### NEURODEGENERATION  ###############################################
# ##########################################################################################################
# 
# ### Determine the the trend during t9 and compare to the baseline toxic level. 
# 
# ### NORMALIZATION 
# strain_summary_normalize = function(strain.summary, rel.time = "t0") {
#   strain.summary$rel.error = strain.summary$sem.rawint/strain.summary$mean.rawint
#   setDT(strain.summary)
#   set(strain.summary, j = "mean.rawint", value = as.numeric(strain.summary[["mean.rawint"]]))
#   set(strain.summary, j = "rel.error", value = as.numeric(strain.summary[["rel.error"]]))
#   set(strain.summary, j = "obs.rawint", value = as.numeric(strain.summary[["obs.rawint"]]))
#   
#   strain.summary.norm = strain.summary[,`:=`(mean.rawint = mean.rawint / mean.rawint[time == rel.time],
#                                            rel.error = rel.error + rel.error[time == rel.time],
#                                            obs.rawint = obs.rawint + obs.rawint[time == rel.time]),
#                                      by = control]
#   return(strain.summary.norm)
# }
# 
# ## Select Low concentrations only
# merge_norm_6h_9h = function(df_6h, df_9h, conc = "High") {
#   a = subset(df_6h, concentration == conc | control %in% c("T", "DMSO", "EGCG"))
#   b = subset(df_9h, concentration == conc | control %in% c("T", "DMSO", "EGCG"))
#   a$time = "6h"
#   b$time = "9h"
#   x = rbind(a,b)
#   x$sem = x$mean.rawint * x$rel.error
#   return(x)
# }
# 
# ###############################
# ### FINAL PLOTS
# ### Summary dataframe
# ## 1.1 Normalizing each strain (default is rel_time = t0)
# oh441.summary.norm = strain_summary_normalize(oh441.summary, rel.time = "t6")
# lx929.summary.norm = strain_summary_normalize(lx929.summary, rel.time = "t6")
# bz555.summary.norm = strain_summary_normalize(bz555.summary, rel.time = "t6")
# eg1285.summary.norm = strain_summary_normalize(eg1285.summary, rel.time = "t6")
# 
# ## 1.2. Get Post-6h (Normalized t0/t6)
# oh441.summary.norm.6h = oh441.summary.norm[oh441.summary.norm$time == "t0" & oh441.summary.norm$control != "U"]
# lx929.summary.norm.6h = lx929.summary.norm[lx929.summary.norm$time == "t0"]
# bz555.summary.norm.6h = bz555.summary.norm[bz555.summary.norm$time == "t0"]
# eg1285.summary.norm.6h = eg1285.summary.norm[eg1285.summary.norm$time == "t0"]
# 
# 
# ## 1.1 Normalizing each strain (default is rel_time = t0)
# oh441.summary.norm = strain_summary_normalize(oh441.summary, rel.time = "t9")
# lx929.summary.norm = strain_summary_normalize(lx929.summary, rel.time = "t9")
# bz555.summary.norm = strain_summary_normalize(bz555.summary, rel.time = "t9")
# eg1285.summary.norm = strain_summary_normalize(eg1285.summary, rel.time = "t9")
# 
# ## 1.2. Get Post-6h (Normalized t0/t6)
# oh441.summary.norm.9h = oh441.summary.norm[oh441.summary.norm$time == "t0" & oh441.summary.norm$control != "U"]
# lx929.summary.norm.9h = lx929.summary.norm[lx929.summary.norm$time == "t0"]
# bz555.summary.norm.9h = bz555.summary.norm[bz555.summary.norm$time == "t0"]
# eg1285.summary.norm.9h = eg1285.summary.norm[eg1285.summary.norm$time == "t0"]
# 
# ####################### FOR High CONCENTRATION ##################
# conc = "High"
# oh441.low.norm = merge_norm_6h_9h(oh441.summary.norm.6h, oh441.summary.norm.9h)
# lx929.low.norm = merge_norm_6h_9h(lx929.summary.norm.6h, lx929.summary.norm.9h)
# bz555.low.norm = merge_norm_6h_9h(bz555.summary.norm.6h, bz555.summary.norm.9h)
# eg1285.low.norm = merge_norm_6h_9h(eg1285.summary.norm.6h, eg1285.summary.norm.9h)
# 
# 
# ####################################################
# #################### EG1285 ########################
# ####################################################
# ### For t=6h
# t0.eg1285 = subset(eg1285.summary, time=="t6")
# t0.eg1285$rep = 1
# eg1285.summary.raw <- subset(eg1285.df, time!= "t6") %>% 
#   group_by(control, time, group, neuron, concentration, rep) %>%
#   summarise(mean.rawint = mean(RawIntDen),
#             sem.rawint = std.error(RawIntDen),
#             obs.rawint = n()) %>%
#   as.data.frame()
# eg1285.summary.raw = rbind(eg1285.summary.raw, t0.eg1285)
# eg1285.raw.norm = strain_summary_normalize(eg1285.summary.raw, rel.time = "t6")
# eeg1285.raw.norm.6h = eg1285.raw.norm[eg1285.raw.norm$time == "t0" & (eg1285.raw.norm$concentration=="High" | 
#                                                                         eg1285.raw.norm$concentration=="Control") & eg1285.raw.norm$control != "Metab"]
# ### For t=9h
# t0.eg1285 = subset(eg1285.summary, time=="t9")
# t0.eg1285$rep = 1
# eg1285.summary.raw <- subset(eg1285.df, time!= "t9") %>% 
#   group_by(control, time, group, neuron, concentration, rep) %>%
#   summarise(mean.rawint = mean(RawIntDen),
#             sem.rawint = std.error(RawIntDen),
#             obs.rawint = n()) %>%
#   as.data.frame()
# eg1285.summary.raw = rbind(eg1285.summary.raw, t0.eg1285)
# eg1285.raw.norm = strain_summary_normalize(eg1285.summary.raw, rel.time = "t9")
# eeg1285.raw.norm.9h = eg1285.raw.norm[eg1285.raw.norm$time == "t0"]
# 
# eg1285.summary.norm.6h = subset(eg1285.summary.norm.6h, control!="Metab" & concentration %in% c("Control", "High") )
# 
# ## Check if there is difference among 6h
# summary(aov(formula = mean.rawint ~ control, data = eeg1285.raw.norm.6h)) ## ANOVA (p = 0.326)
# DunnettTest(x=eeg1285.raw.norm.6h$mean.rawint, g=eeg1285.raw.norm.6h$control) 
# 
# ## Within griup mean difference
# ## Use the following to simulate a gaussian curve
# t_test_means_6h_9h = function(control.norm.9h, ctrl.6h, ctrl.9h, num=10^5) {
#   mean.6h <- c(rnorm(num, mean = round(subset(control.norm.9h, control == ctrl.6h & time == "6h")[,"mean.rawint"],2)$mean.rawint,
#                      sd = round(subset(control.norm.9h, control == ctrl.6h & time == "6h")[,"sem"],2)$sem
#   ))
#   
#   mean.9h <- c(rnorm(num, mean = round(subset(control.norm.9h, control == ctrl.9h & time == "9h")[,"mean.rawint"],2)$mean.rawint,
#                      sd = round(subset(control.norm.9h, control == ctrl.9h & time == "9h")[,"sem"],2)$sem
#   ))  
#   x = t.test(mean.6h, mean.9h, alternative = "greater")
#   return(x)
# }
# 
# t_test_means_6h_9h(eg1285.low.norm, "T", "T") ## p < 0.001
# t.test(subset(eeg1285.raw.norm.6h, control == "DMSO")[,"mean.rawint"],
#        subset(eeg1285.raw.norm.9h, control == "DMSO")[,"mean.rawint"], alternative="two.sided") ## p = 0.568
# t.test(subset(eeg1285.raw.norm.6h, control == "EGCG")[,"mean.rawint"],
#        subset(eeg1285.raw.norm.9h, control == "EGCG")[,"mean.rawint"], alternative="two.sided") ## p = 0.8793
# t.test(subset(eeg1285.raw.norm.6h, control == "FWH")[,"mean.rawint"],
#        subset(eeg1285.raw.norm.9h, control == "FWH")[,"mean.rawint"], alternative="two.sided") ## p = 0.06356
# t.test(subset(eeg1285.raw.norm.6h, control == "LWH")[,"mean.rawint"],
#        subset(eeg1285.raw.norm.9h, control == "LWH")[,"mean.rawint"], alternative="two.sided") ## p = 0.5545
# # t.test(subset(eeg1285.raw.norm.6h, control == "FLH")[,"mean.rawint"],
# #        subset(eeg1285.raw.norm.9h, control == "FLH")[,"mean.rawint"], alternative="two.sided") ## p = 0.5545
# # t.test(subset(eeg1285.raw.norm.6h, control == "LLH")[,"mean.rawint"],
# #        subset(eeg1285.raw.norm.9h, control == "LLH")[,"mean.rawint"], alternative="two.sided") ## p = 0.5545
# 
# # t_test_means_6h_9h(eg1285.low.norm, "LWH", "LWH") ## p < 0.001
# t_test_means_6h_9h(eg1285.low.norm, "FLH", "FLH") ## p < 0.001
# t_test_means_6h_9h(eg1285.low.norm, "LLH", "LLH") ## p = 1
# 
# control.norm.9h = eg1285.low.norm
# ff1 = control.norm.9h[control.norm.9h$control=="T" & control.norm.9h$time=="6h","mean.rawint"]$mean.rawint
# # control.norm.9h = subset(control.norm.9h, control != "DMSO")
# y1 = ggplot(control.norm.9h, aes(control, mean.rawint, fill=time)) +
#   theme_classic() + labs(x = '', y = 'Neurodegeneration') +
#   geom_bar(position=position_dodge(), stat="identity", colour='black') +
#   geom_errorbar(aes(ymin= mean.rawint * (1-rel.error), 
#                     ymax= mean.rawint * (1+rel.error)),
#                 width=.2,
#                 position=position_dodge(.9)) +
#   
#   scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>", "0.1% <br>DMSO", "200 uM <br> EGCG", 
#                               "Field\n Wine","Lab\n Wine", 
#                               "Field\n Lees", "Lab\nLees")) +
#   geom_hline(aes(yintercept = ff1,color = "U"), linetype='dashed', color='#960101') +
#   theme(axis.text.x = element_markdown(color="black", size=12),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.text.y = element_text(color = 'black', size=10),
#         legend.key.size = unit(0.8, 'cm'),
#         legend.text = element_text(size=12),
#         axis.title.y = element_text(size=15),
#         plot.margin=unit(c(1.5,0,0,0.3), 'cm'),
#         legend.position="none"
#         
#   ) +
#   scale_fill_manual("12% NT Exposure",values=c("#d570a8", '#f1ae77'), 
#                     labels = c('6 hrs', '9 hrs')) 
# # geom_hline(aes(yintercept = toxic.recovery.mean, color = "U"), linetype='dashed', color='#960101') +
# # scale_color_manual(values="#960101")
# y1
# 
# 
# op50 <- tibble(
#   x = c(0.77, 0.77, 1.23, 1.23),
#   y = c(4.5, 4.8, 4.8, 3.6)
# ) ## p < 0.01
# 
# fll <- tibble(
#   x = c(4.77, 4.77, 5.23, 5.23),
#   y = c(1.12, 1.38, 1.38, 0.4)
# ) ## p < 0.01
# 
# lll <- tibble(
#   x = c(5.77, 5.77, 6.23, 6.23),
#   y = c(4.29, 4.5,4.5, 1.4)
# ) ## p < 0.01
# 
# 
# y1 = y1 + 
#   annotate("text", x = 1, y = 4.9, label = "***", size = 6, color = "black") +
#   geom_line(data = op50, aes(x= x, y = y, group=1), inherit.aes = F) +
#   annotate("text", x = 5, y = 1.45, label = "***", size = 6, color = "black") +
#   geom_line(data = fll, aes(x= x, y = y, group=1), inherit.aes = F) +
#   annotate("text", x = 6, y = 4.58, label = "***", size = 6, color = "black") +
#   geom_line(data = lll, aes(x= x, y = y, group=1), inherit.aes = F) +
#   annotate("text", x = 1.5, y = ff1+0.4, label = format(round(ff1,2)),
#            size = 3.1, color = '#960101') 
# y1
# png(filename = file.path(getwd(),"/Paper-EG1285-6h-9h-2.png"), width = 7, height = 3, units = "in", res = 600)
# plot(y1)
# dev.off()
# print(paste("Saved at ", getwd()))
# 
# 
# color_scheme = c("#c750c1","#c750c1","#c750c1", "#e65089","#ff9e5f", "#ff7271","#ffcd5f" )
# control.norm.9h$control = factor(control.norm.9h$control, levels = c('T', 'DMSO', 'EGCG', 
#                                                                      'FWL', 'LWL','FLL', 'LLL'),
#                          labels = c("OP50", "0.1%\nDMSO", "0.01%\nEGCG", 
#                                     "Field\n Wine","Lab\n Wine", 
#                                     "Field\n Lees", "Lab\nLees")
# )
# q1 = ( ggplot(control.norm.9h, aes(x=time, y=mean.rawint, fill=control)) +
#          theme_classic()  + labs(x = '', y = 'Neurodegeneration') +
#          geom_bar(stat="identity", color="black", aes(alpha = time)) +
#          geom_errorbar(aes(ymin=mean.rawint-sem, ymax=mean.rawint+sem), width=.2,position=position_dodge(.9)) +
#          facet_grid('. ~ control', scales="free", switch = "x") + 
#          scale_fill_manual(values=color_scheme) +
#          theme(strip.placement = "outside") +
#          theme(axis.title.x = element_blank(),
#                axis.title.y =element_text(size=15),
#                axis.text.x = element_text(size=12),
#                axis.text.y = element_text(size=11),
#                legend.position="none",
#                
#                strip.text.x = element_text(size=12),
#                strip.background = element_blank()) +
#          
#          scale_alpha_discrete(range=c(c(0.9, 0.5))) +
#          guides(alpha="none")
#        
#        # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
#        #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
#        #       strip_text_x = element_text(size=12), strip_background = element_blank(),
#        #       axis_line_y = element_line(colour = 'black', linetype='solid'))
#        
#        # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
#        
# )
# q1
# 
# # op50 <- tibble(
# #   x = c(, 0.77, 1.23, 1.23),
# #   y = c(4.5, 4.8, 4.8, 3.6)
# # ) ## p < 0.01
# # 
# # fll <- tibble(
# #   x = c(5.77, 5.77, 6.23, 6.23),
# #   y = c(1.12, 1.38, 1.38, 0.4)
# # ) ## p < 0.01
# # 
# # lll <- tibble(
# #   x = c(6.77, 6.77, 7.23, 7.23),
# #   y = c(4.29, 4.5,4.5, 1.4)
# # ) ## p < 0.01
# # 
# # q1 = q1 + 
# #   annotate("text", x = 1, y = 4.9, label = "***", size = 6, color = "black") +
# #   geom_line(data = op50, aes(x= x, y = y, group=1), inherit.aes = F) +
# #   annotate("text", x = 6, y = 1.45, label = "***", size = 6, color = "black") +
# #   geom_line(data = fll, aes(x= x, y = y, group=1), inherit.aes = F) +
# #   annotate("text", x = 7, y = 4.58, label = "***", size = 6, color = "black") +
# #   geom_line(data = lll, aes(x= x, y = y, group=1), inherit.aes = F)
# # q1
# 
# ########################################################################################################
# 
# ####################################################
# #################### BZ555 ########################
# ####################################################
# ### For t=6h
# t0.bz555 = subset(bz555.summary, time=="t6")
# t0.bz555$rep = 1
# bz555.summary.raw <- subset(bz555.df, time!= "t6") %>% 
#   group_by(control, time, group, neuron, concentration, rep) %>%
#   summarise(mean.rawint = mean(RawIntDen),
#             sem.rawint = std.error(RawIntDen),
#             obs.rawint = n()) %>%
#   as.data.frame()
# bz555.summary.raw = rbind(bz555.summary.raw, t0.bz555)
# bz555.raw.norm = strain_summary_normalize(bz555.summary.raw, rel.time = "t6")
# bz555.raw.norm.6h = bz555.raw.norm[bz555.raw.norm$time == "t0" & (bz555.raw.norm$concentration=="Low" | 
#                                                                         bz555.raw.norm$concentration=="Control") & bz555.raw.norm$control != "Metab"]
# 
# ### For t=9h
# t0.bz555 = subset(bz555.summary, time=="t9")
# t0.bz555$rep = 1
# bz555.summary.raw <- subset(bz555.df, time!= "t9") %>% 
#   group_by(control, time, group, neuron, concentration, rep) %>%
#   summarise(mean.rawint = mean(RawIntDen),
#             sem.rawint = std.error(RawIntDen),
#             obs.rawint = n()) %>%
#   as.data.frame()
# bz555.summary.raw = rbind(bz555.summary.raw, t0.bz555)
# bz555.raw.norm = strain_summary_normalize(bz555.summary.raw, rel.time = "t9")
# bz555.raw.norm.9h = bz555.raw.norm[bz555.raw.norm$time == "t0"]
# 
# bz555.summary.norm.6h = subset(bz555.summary.norm.6h, control!="Metab" & concentration %in% c("Control", "Low") )
# 
# ## Check if there is difference among 6h
# summary(aov(formula = mean.rawint ~ control, data = bz555.raw.norm.6h)) ## ANOVA (p = 0.326)
# DunnettTest(x=bz555.raw.norm.6h$mean.rawint, g=bz555.raw.norm.6h$control) 
# 
# ## Within group mean difference
# t_test_means_6h_9h(bz555.low.norm, "T", "T") ## p < 0.001
# t.test(subset(bz555.raw.norm.6h, control == "DMSO")[,"mean.rawint"],
#        subset(bz555.raw.norm.9h, control == "DMSO")[,"mean.rawint"], alternative="two.sided") ## p = 0.08
# t_test_means_6h_9h(bz555.low.norm, "EGCG", "EGCG") ## p < 0.001
# t_test_means_6h_9h(bz555.low.norm, "EGCG", "EGCG") ## p < 0.001
# t_test_means_6h_9h(bz555.low.norm, "LWL", "LWL") ## p = 1.00
# t_test_means_6h_9h(bz555.low.norm, "FLL", "FLL") ## p < 0.001
# t_test_means_6h_9h(bz555.low.norm, "LLL", "LLL", alternative = "two.sided") ## p = 1.00
# 
# control.norm.9h = bz555.low.norm
# ff2 = control.norm.9h[control.norm.9h$control=="T" & control.norm.9h$time=="6h","mean.rawint"]$mean.rawint
# control.norm.9h$sem = control.norm.9h$mean.rawint*control.norm.9h$rel.error
# control.norm.9h = subset(control.norm.9h, control != "DMSO")
# y2 = ggplot(control.norm.9h, aes(control, mean.rawint, fill=time)) +
#   theme_classic() + labs(x = '', y = 'Neurodegeneration') +
#   geom_bar(position=position_dodge(), stat="identity", colour='black') +
#   geom_errorbar(aes(ymin= mean.rawint * (1-rel.error), 
#                     ymax= mean.rawint * (1+rel.error)),
#                 width=.2,
#                 position=position_dodge(.9)) +
#   geom_hline(aes(yintercept = ff2,color = "U"), linetype='dashed', color='#960101') +
#   
#   scale_x_discrete(labels = c("OP50", "0.01%\nEGCG", 
#                               "Field\n Wine","Lab\n Wine", 
#                               "Field\n Lees", "Lab\nLees")) +
#   theme(axis.text.x = element_text(color="black", size=12),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.text.y = element_text(color = 'black', size=10),
#         legend.key.size = unit(0.8, 'cm'),
#         legend.text = element_text(size=12),
#         axis.title.y = element_text(size=15),
#         plot.margin=unit(c(1.5,0,0,0.3), 'cm'),
#         legend.position="none"
#   ) +
#   scale_fill_manual("12% NT Exposure",values=c("#d570a8", '#f1ae77'), 
#                     labels = c('6 hrs', '9 hrs')) 
# # geom_hline(aes(yintercept = toxic.recovery.mean, color = "U"), linetype='dashed', color='#960101') +
# # scale_color_manual(values="#960101")
# y2
# 
# 
# op50 <- tibble(
#   x = c(0.77, 0.77, 1.23, 1.23),
#   y = c(2.0, 2.5, 2.5, 1.2)
# ) ## p < 0.001
# 
# egcg <- tibble(
#   x = c(1.77, 1.77, 2.23, 2.23),
#   y = c(2.1, 2.5, 2.5, 1.3)
# ) ## p < 0.001
# 
# fwl <- tibble(
#   x = c(2.77, 2.77, 3.23, 3.23),
#   y = c(2.5, 2.7, 2.7, 1.5)
# ) ## p < 0.001
# 
# fll <- tibble(
#   x = c(4.77, 4.77, 5.23, 5.23),
#   y = c(2.2, 2.8, 2.8, 1.4)
# ) ## p < 0.001
# 
# lwl <- tibble(
#   x = c(3.77, 3.77, 4.23, 4.23),
#   y = c(12.5, 12.7, 12.7, 1.2)
# ) ## p < 0.01
# 
# 
# y2 = y2 + 
#   annotate("text", x = 1, y = 2.7, label = "***", size = 6, color = "black") +
#   geom_line(data = op50, aes(x= x, y = y, group=1), inherit.aes = F) +
#   annotate("text", x = 2, y = 2.7, label = "***", size = 6, color = "black") +
#   geom_line(data = egcg, aes(x= x, y = y, group=1), inherit.aes = F) +
#   annotate("text", x = 3, y = 2.9, label = "***", size = 6, color = "black") +
#   geom_line(data = fwl, aes(x= x, y = y, group=1), inherit.aes = F) +
#   annotate("text", x = 5, y = 3.0, label = "***", size = 6, color = "black") +
#   geom_line(data = fll, aes(x= x, y = y, group=1), inherit.aes = F) +
#   annotate("text", x = 4, y = 12.9, label = "***", size = 6, color = "black") +
#   geom_line(data = lwl, aes(x= x, y = y, group=1), inherit.aes = F) +
#   annotate("text", x = 1.39, y = ff2+0.8, label = format(round(ff2,2)),
#            size = 3.1, color = '#960101') 
# 
# y2
# png(filename = file.path(getwd(),"/Paper-bz555-6h-9h.png"), width = 7, height = 3, units = "in", res = 600)
# plot(y2)
# dev.off()
# print(paste("Saved at ", getwd()))
# ########################################################################################################
# 
# ####################################################
# #################### LX929 ########################
# ####################################################
# ### For t=6h
# t0.lx929 = subset(lx929.summary, time=="t6")
# t0.lx929$rep = 1
# lx929.summary.raw <- subset(lx929.df, time!= "t6") %>% 
#   group_by(control, time, group, neuron, concentration, rep) %>%
#   summarise(mean.rawint = mean(RawIntDen),
#             sem.rawint = std.error(RawIntDen),
#             obs.rawint = n()) %>%
#   as.data.frame()
# lx929.summary.raw = rbind(lx929.summary.raw, t0.lx929)
# lx929.raw.norm = strain_summary_normalize(lx929.summary.raw, rel.time = "t6")
# lx929.raw.norm.6h = lx929.raw.norm[lx929.raw.norm$time == "t0" & (lx929.raw.norm$concentration=="Low" | 
#                                                                     lx929.raw.norm$concentration=="Control") & lx929.raw.norm$control != "Metab"]
# 
# ### For t=9h
# t0.lx929 = subset(lx929.summary, time=="t9")
# t0.lx929$rep = 1
# lx929.summary.raw <- subset(lx929.df, time!= "t9") %>% 
#   group_by(control, time, group, neuron, concentration, rep) %>%
#   summarise(mean.rawint = mean(RawIntDen),
#             sem.rawint = std.error(RawIntDen),
#             obs.rawint = n()) %>%
#   as.data.frame()
# lx929.summary.raw = rbind(lx929.summary.raw, t0.lx929)
# lx929.raw.norm = strain_summary_normalize(lx929.summary.raw, rel.time = "t9")
# lx929.raw.norm.9h = lx929.raw.norm[lx929.raw.norm$time == "t0"]
# 
# lx929.summary.norm.6h = subset(lx929.summary.norm.6h, control!="Metab" & concentration %in% c("Control", "Low") )
# 
# ## Check if there is difference among 6h
# summary(aov(formula = mean.rawint ~ control, data = lx929.raw.norm.6h)) ## ANOVA (p = 0.365)
# DunnettTest(x=lx929.raw.norm.6h$mean.rawint, g=lx929.raw.norm.6h$control) 
# 
# ## Within group mean difference
# t.test(subset(lx929.raw.norm.6h, control == "T")[,"mean.rawint"],
#        subset(lx929.raw.norm.9h, control == "T")[,"mean.rawint"], alternative="greater") ## p = 0.3498
# t.test(subset(lx929.raw.norm.6h, control == "DMSO")[,"mean.rawint"],
#        subset(lx929.raw.norm.9h, control == "DMSO")[,"mean.rawint"], alternative="greater") ## p = 0.7853
# t.test(subset(lx929.raw.norm.6h, control == "EGCG")[,"mean.rawint"],
#        subset(lx929.raw.norm.9h, control == "EGCG")[,"mean.rawint"], alternative="greater") ## p = 0.7679
# t.test(subset(lx929.raw.norm.6h, control == "FWL")[,"mean.rawint"],
#        subset(lx929.raw.norm.9h, control == "FWL")[,"mean.rawint"], alternative="greater") ## p = 0.9988
# t.test(subset(lx929.raw.norm.6h, control == "LWL")[,"mean.rawint"],
#        subset(lx929.raw.norm.9h, control == "LWL")[,"mean.rawint"], alternative="greater") ## p = 0.1266
# t.test(subset(lx929.raw.norm.6h, control == "FLL")[,"mean.rawint"],
#        subset(lx929.raw.norm.9h, control == "FLL")[,"mean.rawint"], alternative="two.sided") ## p = 0.4745
# t.test(subset(lx929.raw.norm.6h, control == "LLL")[,"mean.rawint"],
#        subset(lx929.raw.norm.9h, control == "LLL")[,"mean.rawint"], alternative="two.sided") ## p = 0.1386
# 
# control.norm.9h = lx929.low.norm
# ff3 = control.norm.9h[control.norm.9h$control=="T" & control.norm.9h$time=="6h","mean.rawint"]$mean.rawint
# control.norm.9h$sem = control.norm.9h$mean.rawint*control.norm.9h$rel.error
# control.norm.9h = subset(control.norm.9h, control != "DMSO")
# y3 = ggplot(control.norm.9h, aes(control, mean.rawint, fill=time)) +
#   theme_classic() + labs(x = '', y = 'Neurodegeneration') +
#   geom_bar(position=position_dodge(), stat="identity", colour='black') +
#   geom_errorbar(aes(ymin= mean.rawint * (1-rel.error), 
#                     ymax= mean.rawint * (1+rel.error)),
#                 width=.2,
#                 position=position_dodge(.9)) +
#   geom_hline(aes(yintercept = ff3,color = "U"), linetype='dashed', color='#960101') +
#   
#   scale_x_discrete(labels = c("OP50", "0.01%\nEGCG", 
#                               "Field\n Wine","Lab\n Wine", 
#                               "Field\n Lees", "Lab\nLees")) +
#   theme(axis.text.x = element_text(color="black", size=12),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.text.y = element_text(color = 'black', size=10),
#         legend.key.size = unit(0.8, 'cm'),
#         legend.text = element_text(size=12),
#         axis.title.y = element_text(size=15),
#         plot.margin=unit(c(1.5,0,0,0.3), 'cm'),
#         legend.position="none"
#   ) +
#   scale_fill_manual("12% NT Exposure",values=c("#d570a8", '#f1ae77'), 
#                     labels = c('6 hrs', '9 hrs')) +
#   annotate("text", x = 1.05, y = ff3+0.09, label = format(round(ff3,2)),
#            size = 3.1, color = '#960101') 
# # geom_hline(aes(yintercept = toxic.recovery.mean, color = "U"), linetype='dashed', color='#960101') +
# # scale_color_manual(values="#960101")
# y3
# 
# png(filename = file.path(getwd(),"/Paper-lx929-6h-9h-2.png"), width = 7, height = 3, units = "in", res = 600)
# plot(y3)
# dev.off()
# print(paste("Saved at ", getwd()))
# ########################################################################################################
# 
# ####################################################
# #################### OH441 ########################
# ####################################################
# ### For t=6h
# t0.oh441 = subset(oh441.summary, time=="t6")
# t0.oh441$rep = 1
# oh441.summary.raw <- subset(oh441.df, time!= "t6") %>% 
#   group_by(control, time, group, neuron, concentration, rep) %>%
#   summarise(mean.rawint = mean(RawIntDen),
#             sem.rawint = std.error(RawIntDen),
#             obs.rawint = n()) %>%
#   as.data.frame()
# oh441.summary.raw = rbind(oh441.summary.raw, t0.oh441)
# oh441.raw.norm = strain_summary_normalize(oh441.summary.raw, rel.time = "t6")
# oh441.raw.norm.6h = oh441.raw.norm[oh441.raw.norm$time == "t0" & (oh441.raw.norm$concentration=="Low" | 
#                                                                     oh441.raw.norm$concentration=="Control") & oh441.raw.norm$control != "Metab"]
# 
# ### For t=9h
# t0.oh441 = subset(oh441.summary, time=="t9")
# t0.oh441$rep = 1
# oh441.summary.raw <- subset(oh441.df, time!= "t9") %>% 
#   group_by(control, time, group, neuron, concentration, rep) %>%
#   summarise(mean.rawint = mean(RawIntDen),
#             sem.rawint = std.error(RawIntDen),
#             obs.rawint = n()) %>%
#   as.data.frame()
# oh441.summary.raw = rbind(oh441.summary.raw, t0.oh441)
# oh441.raw.norm = strain_summary_normalize(oh441.summary.raw, rel.time = "t9")
# oh441.raw.norm.9h = oh441.raw.norm[oh441.raw.norm$time == "t0"]
# 
# oh441.summary.norm.6h = subset(oh441.summary.norm.6h, control!="Metab" & concentration %in% c("Control", "Low") )
# 
# ## Check if there is difference among 6h
# summary(aov(formula = mean.rawint ~ control, data = oh441.raw.norm.6h)) ## ANOVA (p = 0.054)
# DunnettTest(x=oh441.raw.norm.6h$mean.rawint, g=oh441.raw.norm.6h$control) 
# 
# ## Within group mean difference
# t_test_means_6h_9h(oh441.low.norm, "T", "T", num=2) ## p = 0.02
# t.test(subset(oh441.raw.norm.6h, control == "DMSO")[,"mean.rawint"],
#        subset(oh441.raw.norm.9h, control == "DMSO")[,"mean.rawint"], alternative="greater") ## p = 0.1648
# t_test_means_6h_9h(oh441.low.norm, "EGCG", "EGCG", num=2) ## p = 0.244
# t.test(subset(oh441.raw.norm.6h, control == "FWL")[,"mean.rawint"],
#        subset(oh441.raw.norm.9h, control == "FWL")[,"mean.rawint"], alternative="greater") ## p = 0.077
# t.test(subset(oh441.raw.norm.6h, control == "LWL")[,"mean.rawint"],
#        subset(oh441.raw.norm.9h, control == "LWL")[,"mean.rawint"], alternative="greater") ## p = 0.09
# t.test(subset(oh441.raw.norm.6h, control == "FLL")[,"mean.rawint"],
#        subset(oh441.raw.norm.9h, control == "FLL")[,"mean.rawint"], alternative="two.sided") ## p = 0.191
# t.test(subset(oh441.raw.norm.6h, control == "LLL")[,"mean.rawint"],
#        subset(oh441.raw.norm.9h, control == "LLL")[,"mean.rawint"], alternative="two.sided") ## p = 0.1492
# 
# control.norm.9h = oh441.low.norm
# ff4 = control.norm.9h[control.norm.9h$control=="T" & control.norm.9h$time=="6h","mean.rawint"]$mean.rawint
# control.norm.9h$sem = control.norm.9h$mean.rawint*control.norm.9h$rel.error
# control.norm.9h = subset(control.norm.9h, control != "DMSO")
# y4 = ggplot(control.norm.9h, aes(control, mean.rawint, fill=time)) +
#   theme_classic() + labs(x = '', y = 'Neurodegeneration') +
#   geom_bar(position=position_dodge(), stat="identity", colour='black') +
#   geom_errorbar(aes(ymin= mean.rawint * (1-rel.error), 
#                     ymax= mean.rawint * (1+rel.error)),
#                 width=.2,
#                 position=position_dodge(.9)) +
#   geom_hline(aes(yintercept = ff4,color = "U"), linetype='dashed', color='#960101') +
#   scale_x_discrete(labels = c("OP50",  "0.01%\nEGCG", 
#                               "Field\n Wine","Lab\n Wine", 
#                               "Field\n Lees", "Lab\nLees")) +
#   theme(axis.text.x = element_text(color="black", size=12),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.text.y = element_text(color = 'black', size=10),
#         legend.key.size = unit(0.8, 'cm'),
#         legend.text = element_text(size=12),
#         axis.title.y = element_text(size=15),
#         plot.margin=unit(c(1.5,0,0,0.3), 'cm'),
#         legend.position="none"
#   ) +
#   scale_fill_manual("12% NT Exposure",values=c("#d570a8", '#f1ae77'), 
#                     labels = c('6 hrs', '9 hrs'))  +
#   annotate("text", x = 1.05, y = ff4+0.1, label = format(round(ff4,2)),
#            size = 3.1, color = '#960101') 
# # geom_hline(aes(yintercept = toxic.recovery.mean, color = "U"), linetype='dashed', color='#960101') +
# # scale_color_manual(values="#960101")
# y4
# 
# png(filename = file.path(getwd(),"/Paper-oh441-6h-9h-2.png"), width = 7, height = 3, units = "in", res = 600)
# plot(y4)
# dev.off()
# print(paste("Saved at ", getwd()))
# ########################################################################################################
# ########################################################################################################
# 
# ## ALL
# png(filename = file.path(getwd(),"/Paper-all-6h-9h-3.png"), width = 10, height = 8, units = "in", res = 600)
# strain.plots = ggarrange(plotlist = list(y4,y3,y2,y1), ncol = 2, nrow = 2)
# plot(strain.plots)
# dev.off()
# print(paste("Saved at ", getwd()))



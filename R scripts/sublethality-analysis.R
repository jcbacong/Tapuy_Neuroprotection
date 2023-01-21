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

######################### READ FILES #############################
sublethal.csv = read.csv("sublethal-data.csv", header = T)
sublethal.csv

sublethal.csv$Concentration <- factor(sublethal.csv$Concentration,
                                      levels = c(8, 10, 12),
                                      labels = c("8%", "10%", "12%"))

sublethal.csv$Day <- factor(sublethal.csv$Day, levels = c(1,2, 3),
                            labels = c("24h", "48h", "72h"))

sublethal.csv$Survival =  sublethal.csv$Alive/sublethal.csv$All

sublethal.csv.summary <- sublethal.csv %>% 
  group_by(Day, Concentration) %>%
  summarise(Survival.mean = mean(Survival),
            Survival.std = std.error(Survival),
            obs.rawint = n()) %>%
  as.data.frame()

sublethal.csv.summary

## Visualize
# color_scheme = c("#92ffc0",  "#5285d1", "#002661")
color_scheme = c("#62cff4",  "#548ecc", "#ff78d6") # Bubble gum
# color_scheme = c("#13c2fc",  "#7f40ff", "#ff78d6") # Bubble gum
q1 = ( ggplot(sublethal.csv.summary, aes(x=Day, y=Survival.mean, fill=Concentration)) +
        theme_classic()  + labs(x = '', y = 'Survival Rate') +
        geom_bar(stat="identity", color="black") +
        geom_errorbar(aes(ymin=Survival.mean-Survival.std, ymax=Survival.mean+Survival.std), width=.2,position=position_dodge(.9)) +
        facet_grid('. ~ Concentration', scales="free", switch = "x") + 
        geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
        scale_fill_manual(values=color_scheme) +
        theme(strip.placement = "outside") +
        theme(axis.title.x = element_blank(),
              axis.title.y =element_text(size=15),
              axis.text.x = element_text(size=11),
              axis.text.y = element_text(size=11, 
                                         color = c("black", "black", "black", "black", "#960101", "black")),
              legend.position="none",
              strip.text.x = element_text(size=13, face = "bold"),
              strip.background = element_blank()) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1)) 

       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
        
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q1

tiff(filename = file.path(getwd(),"/Sublethality-12per2.tiff"), width = 5.5, height = 3.5, units = "in", res = 600)
plot(q1)
dev.off()
print(paste("Saved at ", getwd()))

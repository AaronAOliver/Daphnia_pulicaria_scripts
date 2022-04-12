# phenotype

library(tidyverse)
library(grid)
library(gridExtra)
library(ggmap)
library(broom)
library(lsmeans)
library(cowplot)
library(doBy)

# Temperature profile ----------------------------------------------------------

lake.temp <- data.frame('depth' = c(0:13),
                        'blue' = c(14.1, 13.9, 13.8, 13.8, 13.8, 13.7, 13.6, 13.4,
                                   13, 12.5, 11.6, NA, NA, NA),
                        'gardisky' = c(14, 13.7, 13.6, 13.5, 13.4, 13.1, 12.4, 
                                       10.6, 8.6, 7.7, 7.2, 7, 6.8, 6.9))

plot.temp <- lake.temp %>% 
    gather(lake, temperature, - depth) %>% 
    ggplot() +
    geom_path(aes(temperature, depth, group = lake, color = lake), size = 1.5) +
    scale_y_reverse(lim=c(13, 0)) +
    ylab("Depth (m)") +
    xlab("Temperature (°C)") +
    scale_color_manual(name = "", values = c("steelblue3", "maroon2"), 
                       labels = c("Blue (3013 m)", "Gardisky (3204 m)")) +
    theme(panel.background = element_rect('white'),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey",size = rel(0.5)),
          axis.line = element_line('black'),
          legend.position = c(0.15, 0.9),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15))

plot.temp

png("lake_temp.png", width=5, height=5, res=300, units = 'in', type = 'cairo')
grid.arrange(plot.temp)
dev.off()

# Loading data -----------------------------------------------------------------

#pheno<-read.csv("phenotype.csv")

pheno <- read.csv("phenotype.csv") %>% 
    mutate_at(vars(lake, temperature, fish), funs(factor(.))) %>% 
    rowwise() %>% 
    mutate(avr.off = mean(c(off1, off2, off3), na.rm = TRUE)) %>% 
    mutate(sum.off = sum(c(off1, off2, off3)))
    
pheno <- pheno %>% as.data.frame()

for(i in 1:nrow(pheno)) {
    if(!is.na(pheno$off1[i])) {
        #i <- 44
        x <- c(pheno[i,'age1'], pheno[i,'age2'], pheno[i,'age3'])
        L <- rep(1, 3) 
        m <- c(pheno[i,'off1'], pheno[i,'off2'], pheno[i,'off3']) 
        
        eulerlotka <- function(r) sum(L * m * exp(-r * x)) - 1 
        
        res <- uniroot(f = eulerlotka, interval = c(0, 3), tol = 1e-8) 
        pheno$r[i] <- res$root
    } else {
        pheno$r[i] <- NA
    }
}

pheno <- pheno %>% as.tibble()

# Size  -----------------------------------------------------------------

## distribution
ggplot(pheno) +
    geom_histogram(aes(size))

ggplot(pheno) +
    geom_histogram(aes(log(size))) # log is better

md1 <- glm(log(size) ~ temperature * fish * lake, data = pheno)
anova(md1, test = "Chisq")
cld(lsmeans::lsmeans(md1, ~ fish * lake * temperature))

install.packages("lmerTest")
library(lme4)
library(lmerTest)
md1a<-lmer(log(size)~temperature + fish + temperature*fish + (1|lake), data=pheno) # random effect of lake
md1b<-lm(log(size)~temperature + fish + temperature*fish, data=pheno) # no random effect
anova(md1a, md1b) # testing which model is better, one with random effect is better
summary(md1a)
anova(md1a)


size.graph <- pheno %>% 
    filter(!is.na(size)) %>% 
    group_by(temperature, fish, lake) %>% 
    summarise(avr = mean(size), se = sd(size)/sqrt(length(size))) %>% 
    ggplot() +
    geom_point(aes(temperature, avr, color = fish, shape = fish), alpha = 0.5,
               size = 10, show.legend = FALSE) +
    geom_errorbar(aes(temperature, ymin = avr - se, ymax = avr + se, 
                      color = fish), width = 0, size = 2, show.legend = FALSE) +
    geom_line(aes(temperature, avr, color = fish, group = fish), 
              show.legend = FALSE) +
    facet_wrap(~lake) +
    ylab('Size at maturity') +
    xlab('') +
    scale_color_manual(name = "Fish cues", values = c("limegreen", "slateblue2"), 
                       labels = c("No", "Yes")) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey", size = rel(0.5),
                                          linetype = 2),
          panel.background = element_blank(),
          axis.line = element_line('black'),
          axis.title.y = element_text(size = 30),
          axis.text.y = element_text(size = 26),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(face = "bold", size = 17),
          strip.background = element_blank(),
          strip.text = element_blank()) 

size.graph


png("figure_size.png", width=7, height=4, res=300, units = 'in', 
    type = 'cairo')
grid.arrange(size.graph)
dev.off()

# Ctmax  -----------------------------------------------------------------

## distribution
ggplot(pheno) +
    geom_histogram(aes(tolerance))

ggplot(pheno) +
    geom_histogram(aes(log(tolerance))) # log is a bit better

md3 <- glm(log(tolerance) ~ log(size) * temperature * fish * lake, 
           data = pheno)
md3try <- glm(log(tolerance) ~ log(size) + temperature * fish * lake, 
           data = pheno)
anova(md3, md3try, test = "Chisq")
anova(md3try, test = "Chisq")

### Mixed effects models
md3a<-lmer(log(tolerance)~temperature + fish + temperature*fish + (1|lake), data=pheno) # random effect of lake
md3b<-lm(log(tolerance)~temperature + fish + temperature*fish, data=pheno) # no random effect
anova(md3a, md3b) # testing which model is better, one with random effect is better
summary(md3a)
anova(md3a)

cld(lsmeans::lsmeans(md3try, ~ fish * lake * temperature))


ctmax.graph <- pheno %>% 
    filter(test == "ctmax") %>% 
    group_by(temperature, fish, lake) %>% 
    summarise(avr = mean(tolerance), se = sd(tolerance)/sqrt(length(tolerance))) %>% 
    ggplot() +
    geom_point(aes(temperature, avr, color = fish, shape = fish), alpha = 0.5,
               size = 10, show.legend = TRUE) +
    geom_errorbar(aes(temperature, ymin = avr - se, ymax = avr + se, 
                      color = fish), width = 0, size = 2, show.legend = FALSE) + 
    geom_line(aes(temperature, avr, color = fish, group = fish), 
              show.legend = FALSE) +
    facet_wrap(~lake) +
    ylab(expression('CT'[max])) +
    xlab('Temperature (°C)') +
    scale_color_manual(name = "", values = c("limegreen", "slateblue2"), 
                       labels = c("N", "Y")) +
    scale_shape_manual(name = "", values = c(19, 17), 
                       labels = c("N", "Y")) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey", size = rel(0.5),
                                          linetype = 2),
          panel.background = element_blank(),
          axis.line = element_line('black'),
          legend.position = c(0, 0.9),
          legend.background = element_blank(),
          legend.key = element_rect(size = 3),
          legend.key.size = unit(2.5, 'lines'),
          legend.text = element_text(size = 26),
          axis.title = element_text(size = 30),
          axis.text = element_text(size = 26),
          plot.title = element_text(face = "bold", size = 17),
          strip.background = element_blank(),
          strip.text.x = element_blank()) 

ctmax.graph    

png("figure_ctmax.png", width=7, height=4, res=300, units = 'in', 
    type = 'cairo')
grid.arrange(ctmax.graph)
dev.off()


# Age  -----------------------------------------------------------------

## distribution
ggplot(pheno) +
    geom_histogram(aes(maturity))

ggplot(pheno) +
    geom_histogram(aes(log(maturity))) # log is just a bit better


md4 <- glm(log(maturity) ~ temperature * fish * lake, data = pheno)
anova(md4,  test = "Chisq")
cld(lsmeans::lsmeans(md4, ~ fish * lake * temperature))

### Mixed effects models
md4a<-lmer(log(maturity)~temperature + fish + temperature*fish + (1|lake), data=pheno) # random effect of lake
md4b<-lm(log(maturity)~temperature + fish + temperature*fish, data=pheno) # no random effect
anova(md4a, md4b) # testing which model is better, one with random effect is better
summary(md4a)
anova(md4a, test = "Chisq")

md4ag<-glmer(log(maturity)~temperature + fish + temperature*fish + (1|lake), data=pheno) # random effect of lake
summary(md4ag)
anova(md4ag)

age.graph <- pheno %>% 
    filter(!is.na(maturity)) %>% 
    group_by(temperature, fish, lake) %>% 
    summarise(avr = mean(maturity), se = sd(maturity)/sqrt(length(maturity))) %>% 
    ggplot() +
    geom_point(aes(temperature, avr, color = fish, shape = fish), alpha = 0.5,
               size = 10, show.legend = FALSE) +
    geom_errorbar(aes(temperature, ymin = avr - se, ymax = avr + se, 
                      color = fish), width = 0, show.legend = FALSE, size = 2) + 
    geom_line(aes(temperature, avr, color = fish, group = fish), 
              show.legend = FALSE) +
    facet_wrap(~lake) +
    ylab('Age at maturity') +
    ylim(c(4, 9.7)) +
    scale_color_manual(name = "Fish cues", values = c("limegreen", "slateblue2"), 
                       labels = c("No", "Yes")) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey",size = rel(0.5),
                                          linetype = 2),
          panel.background = element_blank(),
          axis.line = element_line('black'),
          axis.title.y = element_text(size = 30),
          axis.text.y = element_text(size = 26),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(face = "bold", size = 17),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    annotate('text', x = 1.5, y = 9.6, label = c('Blue', 'Gardisky'), size = 13,
             fontface = 'bold')

age.graph    

png("figure_age.png", width=7, height=4, res=300, units = 'in', 
    type = 'cairo')
grid.arrange(age.graph)
dev.off()


# r  -----------------------------------------------------------------

## distribution
ggplot(pheno) +
    geom_histogram(aes(r)) # r is already logged

md5 <- glm(r ~ temperature * fish * lake, data = pheno)
anova(md5, test = "Chisq")
cld(lsmeans::lsmeans(md5, ~ fish * lake * temperature))

## Mixed effects models
md5a<-lmer(r~temperature + fish + temperature*fish + (1|lake), data=pheno) # random effect of lake
md5b<-lm(r~temperature + fish + temperature*fish, data=pheno) # no random effect
anova(md5a, md5b) # testing which model is better, one with random effect is better
summary(md5a)
summary(md5b)

anova(md5a)

r.graph <- pheno %>% 
    filter(!is.na(r)) %>% 
    group_by(temperature, fish, lake) %>% 
    summarise(avr = mean(r), se = sd(r)/sqrt(length(r))) %>% 
    ggplot() +
    geom_point(aes(temperature, avr, color = fish,shape = fish), alpha = 0.5,
               size = 10, show.legend = FALSE) +
    geom_errorbar(aes(temperature, ymin = avr - se, ymax = avr + se, 
                      color = fish), width = 0, size = 2, show.legend = FALSE) + 
    geom_line(aes(temperature, avr, color = fish, group = fish), 
              show.legend = FALSE) +
    facet_wrap(~lake) +
    ylab('r') +
    xlab('Temperature (°C)') +
    scale_color_manual(name = "Fish cues", values = c("limegreen", "slateblue2"), 
                       labels = c("No", "Yes")) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey",size = rel(0.5),
                                          linetype = 2),
          panel.background = element_blank(),
          axis.line = element_line('black'),
          axis.title.y = element_text(size = 30),
          axis.text.y = element_text(size = 26),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank())

r.graph    

# Offspring  -----------------------------------------------------------------

## distribution
ggplot(pheno) +
    geom_histogram(aes(avr.off))

ggplot(pheno) +
    geom_histogram(aes(log(avr.off))) # log is better


# sum of the 3 clutches give the same result
md6 <- glm(log(avr.off) ~ temperature * fish * lake, data = pheno)
anova(md6, test = "Chisq")
cld(lsmeans::lsmeans(md6, ~ fish * lake * temperature))

## Mixed effects models
md6a<-lmer(log(avr.off)~temperature + fish + temperature*fish + (1|lake), data=pheno) # random effect of lake
md6b<-lm(log(avr.off)~temperature + fish + temperature*fish, data=pheno) # no random effect
anova(md6a, md6b) # testing which model is better, one with random effect is better
summary(md6a)
summary(md6b)

off.graph <- pheno %>% 
    filter(!is.na(avr.off)) %>% 
    group_by(temperature, fish, lake) %>% 
    summarise(avr = mean(avr.off), se = sd(avr.off)/sqrt(length(avr.off))) %>% 
    ggplot() +
    geom_point(aes(temperature, avr, color = fish, shape = fish), alpha = 0.5,
               size = 10, show.legend = FALSE) +
    geom_errorbar(aes(temperature, ymin = avr - se, ymax = avr + se, 
                      color = fish), width = 0, size = 2, show.legend = FALSE) + 
    geom_line(aes(temperature, avr, color = fish, group = fish), 
              show.legend = FALSE) +
    facet_wrap(~lake) +
    ylab('Number of offspring') +
    xlab('Temperature (°C)') +
    scale_color_manual(name = "Fish cues", values = c("limegreen", "slateblue2"), 
                       labels = c("No", "Yes")) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey",size = rel(0.5),
                                          linetype = 2),
          panel.background = element_blank(),
          axis.line = element_line('black'),
          axis.title.y = element_text(size = 30),
          axis.text.y = element_text(size = 26),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())

off.graph    

################################################################################

# Final figure phenotype
png("fig_phenotype.png", width = 10, height = 18, res = 300, units = 'in', 
    type = 'cairo')
plot_grid(age.graph, size.graph, off.graph, r.graph, ctmax.graph,
          ncol = 1, nrow = 5, align = 'v',
          rel_heights = c(rep(1, 4), 1.2))
dev.off()

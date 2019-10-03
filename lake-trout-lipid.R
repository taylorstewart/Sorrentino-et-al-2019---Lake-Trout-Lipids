##############################################################
##############################################################
##  Lake trout lipids (Sorrentino et al.) thesis work
##
##############################################################
##############################################################
## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names=TRUE))

## ===========================================================
## Load Packages
## ===========================================================
library(dplyr)      # manipulating data
library(magrittr)   # for %<>%
library(readxl)     # reading data
library(ggplot2)    # visualizations
library(FSA)        # ANCOVA functions
library(multcomp)   # for glht()
library(car)        # logit transformation


## ===========================================================
## lipid by region ANOVA 
## ===========================================================
ltlipids <- read_excel("lake-trout-lipids-2018.xlsx", sheet = "Data") %>% 
  ## remove problematic rows and hatchery fish
  filter(include == 'y') %>% 
  ## rename species/source variables, modify class to factor variables, and logit transform percent lipid
  mutate(location = factor(location, levels = c("North", "Central", "South", "Grand Isle Hatchery")),
         source = factor(ifelse(.$source == "wild", "Wild", "Stocked")),
         season = factor(season, levels = c("Pre-winter","Spring", "Summer", "Autumn")),
         age_class = ifelse(source == "Wild", age_class-1, age_class),
         lipid.logit = car::logit(avg_perc_lipid, percents = TRUE)) %>% 
  ## remove unused columns - clean data frame
  dplyr::select(sample_id, source, season, location, age_class, tl = length, wt = wet_wgt_stmch, avg_perc_lipid, lipid.logit, season)
str(ltlipids)


## ===========================================================
## Sample Size and Length Range Table
## ===========================================================
ltlipids.summary <- ltlipids %>% group_by(source, season, location) %>% 
  summarize(n = n(),
            min.tl = min(tl),
            max.tl = max(tl),
            mean.tl = mean(tl)
            )
print(ltlipids.summary)


## ===========================================================
## Fit model (lipids~tl)
## ===========================================================
## Remove hatchery fish
ltlipids.filt <- ltlipids %>% filter(location != "Grand Isle Hatchery")
lm.lipids <- lm(lipid.logit ~ tl * source * location * season, data = ltlipids.filt)

##------------------------------------------------------------
## Check Assumptions
##------------------------------------------------------------
##  Normality
shapiro.test(ltlipids.filt$lipid.logit) ## p = 0.1219; accept normality

## Homogeneity of Variance by source
bartlett.test(lipid.logit ~ interaction(source, location, season), data = ltlipids.filt)  ## p = 0.1804; equal variance


##------------------------------------------------------------
## ANOVA
##------------------------------------------------------------
## full model
summary(aov <- aov(lipid.logit ~ tl * source * location * season, data = ltlipids.filt))
      ## ATL = Significant
      ## Source = Significant
      ## Location = Not Significant
      ## Season = Significant
      ## TL x Source = Sig. Interaction
      ## TL x Location = No Interaction
      ## Source x Location = No Interaction
      ## TL x Season = Sig. Interaction
      ## Source x Season = No Interaction
      ## Location x Season = No Interaction
      ## TL x Source x Location = No Interaction
      ## TL x Source x Season = No Interaction
      ## TL x Location x Season = No Interaction
      ## Source x Location x Season = No Interaction
      ## TL x Source x Location x Season = No Interaction


## ===========================================================
## Make plot
## ===========================================================
ggplot(ltlipids, aes(x = tl, y = avg_perc_lipid, color = season, linetype = season)) +
  geom_text(aes(label = age_class), show.legend = FALSE, size = 5, alpha = 0.8) +
  geom_smooth(data = filter(ltlipids, location != "Grand Isle Hatchery"),
              aes(x = tl, y = avg_perc_lipid), 
              method = "lm", se = FALSE) +
  scale_x_continuous(limits = c(76, 340), breaks = seq(100, 300, 50), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 41), breaks = seq(0, 40, 10), expand = c(0, 0)) +
  scale_colour_manual(values = c("black", "#2b83ba", "#d7191c", "#fdae61")) + 
  scale_linetype_manual(values = c("solid", "solid", "dotdash", "dashed")) +
  labs(y = 'Mean % Lipid Content', x = 'Total Length (mm)') +
  theme(axis.text = element_text(size = 20), axis.line.x = element_line(), 
        axis.line.y = element_line(), 
        axis.title.x = element_text(size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(size = 22, margin = margin(0, 15, 0, 0)),
        axis.ticks.length = unit(2.5, 'mm'), 
        legend.text = element_text(size = 14), legend.title = element_blank(),
        legend.key.width = unit(1.2, 'cm'), legend.key.height = unit(0.8, 'cm'), 
        legend.key = element_blank(), legend.position = c(0.6, 0.88),
        legend.spacing.x = unit(0.2, 'cm'),
        strip.background = element_blank(), strip.text = element_text(size = 22),
        panel.background = element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(5, 10, 5, 5), "mm")) +
  facet_wrap(~source, scales = "free_y")

ggsave("figures/Sorrentino_et_al_Fig2.tiff", dpi = 300, width = 13, height = 7)
ggsave("figures/Sorrentino_et_al_Fig2_LowRes.tiff", dpi = 150, width = 13, height = 7)


ggplot(filter(ltlipids, location != "Grand Isle Hatchery"), aes(x = tl, y = avg_perc_lipid, color = source, linetype = source)) +
  geom_text(aes(label = age_class), show.legend = FALSE, size = 5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_manual(values = c("black", "#2b83ba")) + 
  scale_x_continuous(limits = c(76, 340), breaks = seq(100, 300, 50), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 41), breaks = seq(0, 40, 10), expand = c(0, 0)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(y = 'Mean % Lipid Content', x = 'Total Length (mm)') +
  theme(axis.text = element_text(size = 20), axis.line.x = element_line(), 
        axis.line.y = element_line(), 
        axis.title.x = element_text(size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(size = 22, margin = margin(0, 15, 0, 0)),
        axis.ticks.length = unit(2.5, 'mm'), 
        legend.text = element_text(size = 16), legend.title = element_blank(),
        legend.key.width = unit(1.2, 'cm'), legend.key.height = unit(0.8, 'cm'), 
        legend.key = element_blank(), legend.position = c(0.11, 0.94),
        legend.spacing.x = unit(0.3, 'cm'),
        panel.background = element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(5, 7.5, 2, 2), "mm"))

ggsave("figures/Sorrentino_et_al_Fig3.tiff", dpi = 300, width = 10, height = 7)
ggsave("figures/Sorrentino_et_al_Fig3_LowRes.tiff", dpi = 150, width = 10, height = 7)


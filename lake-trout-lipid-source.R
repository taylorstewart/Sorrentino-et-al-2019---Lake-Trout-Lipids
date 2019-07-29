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
library(car)        # leveneTest()
library(boot)


## ===========================================================
## lipid by region ANOVA 
## ===========================================================
ltlipids.source <- read_excel("lake-trout-lipids-2018.xlsx", sheet = "Data") %>% 
  ## remove problematic rows and hatchery fish
  filter(include == 'y',
         location != "Grand Isle Hatchery") %>% 
  ## rename species/source variables, modify class to factor variables, and logit transform percent lipid
  mutate(location = factor(location, levels = c("North", "Central", "South")),
         source = factor(ifelse(.$source == "wild", "Wild", "Stocked")),
         age_class = factor(age_class),
         season = factor(season, levels = c("Spring", "Summer", "Autumn")),
         group = interaction(location:source:age_class),
         lipid.logit = car::logit(avg_perc_lipid, percents = TRUE)) %>% 
  ## remove unused columns - clean data frame
  dplyr::select(sample_id, source, date_collected, date_processed, location, tl = length, wt = wet_wgt_stmch, avg_perc_lipid, lipid.logit, age_class, group, season)
str(ltlipids.source)


## ===========================================================
## Prelim Tables and Plots
## ===========================================================
##------------------------------------------------------------
## Sample Size, Length Range, and Weight Range Table
##------------------------------------------------------------
##
ltlipids.source.summary <- ltlipids.source %>% group_by(source) %>% 
  summarize(n = n(),
            mean.lipids = mean(avg_perc_lipid),
            sd.lipids = sd(avg_perc_lipid),
            se.lipids = sd.lipids/sqrt(n),
            min.tl = min(tl),
            mean.tl = mean(tl),
            max.tl = max(tl),
            #sd.tl = sd(tl),
            min.wt = min(wt),
            mean.wt = mean(wt),
            max.wt = max(wt),
            #sd.wt = sd(wt)
  )

ltlipids.source.loc.summary <- ltlipids.source %>% group_by(source, location) %>% 
  summarize(n = n(),
            mean.lipids = mean(avg_perc_lipid),
            sd.lipids = sd(avg_perc_lipid),
            se.lipids = sd.lipids/sqrt(n),
            min.tl = min(tl),
            mean.tl = mean(tl),
            max.tl = max(tl),
            #sd.tl = sd(tl),
            min.wt = min(wt),
            mean.wt = mean(wt),
            max.wt = max(wt),
            #sd.wt = sd(wt)
  )

## ===========================================================
## Fit model (lipids~tl)
## ===========================================================
lm.lipids <- lm(lipid.logit ~ tl + source * location, data = ltlipids.source)
summary(lm.lipids)

##------------------------------------------------------------
## Check Assumptions
##------------------------------------------------------------
##  Normality
shapiro.test(ltlipids.source$lipid.logit) ## p = 0.1219; accept normality

## Homogeneity of Variance by source
bartlett.test(lipid.logit ~ interaction(source, location), data = ltlipids.source)  ## p = 0.1804; equal variance


##------------------------------------------------------------
## ANOVA
##------------------------------------------------------------
## full model
summary(aov <- aov(lipid.logit ~ tl + source * location, data = ltlipids.source))
      ## Added interaction just to check. No interaction.

## Source Main Effect
mc.source <- glht(lm.lipids, mcp(source = "Tukey"))  ## only two groups...
summary(mc.source)

## Location Main Effect
mc.location <- glht(lm.lipids, mcp(location = "Tukey"))
summary(mc.location)

## Pairwise
ltlipids.source$group <- ltlipids.source$source:ltlipids.source$location
lm.pair <- lm(lipid.logit ~ tl + group, data = ltlipids.source)
mc.pair <- glht(lm.pair, mcp(group = "Tukey"))
summary(mc.pair)
mc.pair$df


## ===========================================================
## Make plot
## ===========================================================
ggplot(ltlipids.source.loc.summary, aes(x = location, y = mean.lipids, group = source)) +
  geom_bar(aes(fill = source), stat = "identity", position = position_dodge2(preserve = "single"), color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean.lipids + se.lipids, ymax = mean.lipids - se.lipids), width = 0.6, size = 0.75, position = position_dodge2(preserve = "single", padding = 0.6)) +
  annotate("text", x = 0.85, y = 0.9, label = "n=15", size = 5, vjust = 1.0) +
  annotate("text", x = 1.15, y = 0.9, label = "n=3 ", size = 5, vjust = 1.0) +
  annotate("text", x = 1.85, y = 0.9, label = "n=44", size = 5, vjust = 1.0) +
  annotate("text", x = 2.15, y = 0.9, label = "n=57", size = 5, vjust = 1.0) +
  annotate("text", x = 2.85, y = 0.9, label = "n=37", size = 5, vjust = 1.0) +
  annotate("text", x = 3.15, y = 0.9, label = "n=26", size = 5, vjust = 1.0) +
  scale_fill_manual(values = c("gray75", "white"), labels = c("Stocked", "Wild")) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.35)) +
  labs(y = 'Mean % Lipid Content', x = '') +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.line.x = element_line(), 
        axis.line.y = element_line(), axis.title.x = element_blank(),
        axis.title.y = element_text(size = 22, margin = margin(0, 15, 0, 0)),
        legend.position = c(0.9, 0.92), legend.text = element_text(size = 18), legend.title = element_blank(),
        legend.key.width = unit(1.0, 'cm'), legend.key.height = unit(1, 'cm'), 
        legend.spacing.x = unit(0.3, 'cm'),
        axis.ticks.length = unit(2.5, 'mm'), plot.margin = unit(c(5, 7.5, 2, 2), "mm"))

ggsave("figures/Sorrentino_et_al_Fig2.tiff", dpi = 300, width = 10, height = 8)
ggsave("figures/Sorrentino_et_al_Fig2_LowRes.tiff", dpi = 150, width = 10, height = 8)

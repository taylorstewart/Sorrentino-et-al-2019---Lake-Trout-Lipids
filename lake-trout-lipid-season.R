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
library(emmeans)    # calculating least-squares means
library(multcomp)   # for glht()
library(AICcmodavg) # for aictab()
library(car)        # leveneTest()


## ===========================================================
## lipid by region ANOVA 
## ===========================================================
ltlipids.season <- read_excel("lake-trout-lipids-2018.xlsx", sheet = "Data") %>% 
  ## remove problematic rows and hatchery fish
  filter(include == 'y', location %in% c("Central", "Grand Isle Hatchery")) %>% 
  ## rename species/source variables and modify class to factor variables
  mutate(location = factor(location),
         source = factor(ifelse(.$source == "wild", "Wild", "Stocked")),
         season = factor(season, levels = c("Pre-winter","Spring", "Summer", "Autumn")),
         lipid.logit = car::logit(avg_perc_lipid, percents = TRUE)) %>% 
  ## remove unused columns - clean data frame
  dplyr::select(sample_id, source, season, date_collected, date_processed, location, tl = length, wt = wet_wgt_stmch, avg_perc_lipid, lipid.logit)
str(ltlipids.season)



##------------------------------------------------------------
## Sample Size, Length Range, and Weight Range Table
##------------------------------------------------------------
ltlipids.season.summary <- ltlipids.season %>% group_by(source, season) %>% 
  summarize(n = n(),
            mean.lipids = mean(avg_perc_lipid),
            sd.lipids = sd(avg_perc_lipid),
            se.lipids = sd.lipids/sqrt(n)) %>% ungroup()


ltlipids.zero <- data.frame(source = "Wild",
                            season = "Pre-winter",
                            n = 0,
                            mean.lipids = 0,
                            sd.lipids = 0,
                            se.lipids = 0)
ltlipids.season.all <- bind_rows(ltlipids.season.summary, ltlipids.zero) %>% 
  mutate(source = factor(source),
         season = factor(season))


## ===========================================================
## Fit model (lipids ~ tl + )
## ===========================================================
## Remove hatchery fish
ltlipids.season.filt <- ltlipids.season %>% filter(location != "Grand Isle Hatchery")

lm.lipids <- lm(lipid.logit ~ tl + season * source, data = ltlipids.season.filt)
summary(lm.lipids)
fitPlot(lm.lipids, legend = "topleft")

##------------------------------------------------------------
## Check Assumptions
##------------------------------------------------------------
##  Normality
shapiro.test(ltlipids.season.filt$lipid.logit) ## p = 0.3209; accept normality

## Homogeneity of Variance by source
## Flinger-Killeen Test is a non-parametric test which is very robust against departures from normality.
fligner.test(lipid.logit ~ interaction(source, season), data = ltlipids.season.filt)  ## p = 0.1178; equal variance

##------------------------------------------------------------
## ANOVA
##------------------------------------------------------------
## full model
summary(aov <- aov(lipid.logit ~ tl + season * source, data = ltlipids.season.filt))
## Added interaction just to check. No interaction.

## Season Main Effect
mc.season <- glht(lm.lipids, mcp(season = "Tukey"))
summary(mc.season)

## Pairwise
ltlipids.season.filt$group <- ltlipids.season.filt$source:ltlipids.season.filt$season
lm.pair <- lm(lipid.logit ~ tl + group, data = ltlipids.season.filt)
mc.pair <- glht(lm.pair, mcp(group = "Tukey"))
summary(mc.pair)
mc.pair$df

## ===========================================================
## Make plot
## ===========================================================
ggplot(ltlipids.season.summary, aes(x = season, y = mean.lipids, group = source)) +
  geom_bar(aes(fill = source), stat = "identity", position = position_dodge2(preserve = "single"), color = "black", width = 0.82) +
  geom_errorbar(aes(ymin = mean.lipids + se.lipids, ymax = mean.lipids - se.lipids), width = 0.82, size = 0.75, position = position_dodge2(preserve = "single", padding = 0.6)) +
  annotate("text", x = 1.0, y = 1.2, label = "n=15", size = 5, vjust = 1.0) + 
  annotate("text", x = 1.78, y = 1.2, label = "n=13", size = 5, vjust = 1.0) +
  annotate("text", x = 2.22, y = 1.2, label = "n=15", size = 5, vjust = 1.0) +
  annotate("text", x = 2.78, y = 1.2, label = "n=16", size = 5, vjust = 1.0) +
  annotate("text", x = 3.22, y = 1.2, label = "n=27", size = 5, vjust = 1.0) +
  annotate("text", x = 3.78, y = 1.2, label = "n=15", size = 5, vjust = 1.0) +
  annotate("text", x = 4.22, y = 1.2, label = "n=15", size = 5, vjust = 1.0) +
  scale_fill_manual(values = c("gray75", "white"), labels = c("Stocked", "Wild")) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.33)) +
  labs(y = 'Mean % Lipid Content', x = '') +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.line.x = element_line(), 
        axis.line.y = element_line(), axis.title.x = element_blank(),
        axis.title.y = element_text(size = 22, margin = margin(0, 15, 0, 0)),
        legend.position = c(0.9, 0.92), legend.text = element_text(size = 18), legend.title = element_blank(),
        legend.key.width = unit(1.0, 'cm'), legend.key.height = unit(1, 'cm'), 
        legend.spacing.x = unit(0.3, 'cm'),
        axis.ticks.length = unit(2.5, 'mm'), plot.margin = unit(c(5, 7.5, 2, 2), "mm"))

ggsave("figures/Sorrentino_et_al_Fig3.tiff", dpi = 300, width = 10, height = 8)
ggsave("figures/Sorrentino_et_al_Fig3_LowRes.tiff", dpi = 150, width = 10, height = 8)

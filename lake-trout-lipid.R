##############################################################
## 
##  Lake trout lipids (Sorrentino et al.) thesis work
##
##############################################################

## CLEAR THE ENVIRONMENT FIRST ==================================

rm(list = ls(all.names=TRUE))


## SET SEED FOR REPRODUCIBILTY ==================================

set.seed(98024790)


## LOAD PACKAGES ================================================

library(dplyr)      # manipulating data
library(magrittr)   # for %<>%
library(readxl)     # reading data
library(ggplot2)    # visualizations
library(FSA)        # ANCOVA functions
library(multcomp)   # for glht()
library(car)        # logit transformation


## READ DATA AND MANIPULATE =====================================

ltlipids <- read_excel("lake-trout-lipids-2018.xlsx", sheet = "Data") %>% 
  ## remove problematic rows and hatchery fish
  filter(include == 'y') %>% 
  ## rename species/source variables, modify class to factor variables, and logit transform percent lipid
  mutate(location = factor(location, levels = c("South", "Central", "North", "Grand Isle Hatchery")),
         source = factor(ifelse(.$source == "wild", "Wild", "Stocked")),
         season = factor(season, levels = c("Pre-winter","Spring", "Summer", "Autumn")),
         age_class = ifelse(source == "Wild", age_class, age_class),
         lipids.prop = avg_perc_lipid/100,
         lipid.logit = car::logit(avg_perc_lipid, percents = TRUE)) %>% 
  ## remove unused columns - clean data frame
  dplyr::select(sample_id, source, season, location, age_class, tl = length, wt = wet_wgt_stmch, avg_perc_lipid, lipids.prop, lipid.logit, season)
str(ltlipids)


## SAMPLE SIZE AND LENGTH RANGE TABLE ===========================

ltlipids.tl.summary <- ltlipids %>% group_by(source, season, location) %>% 
  summarize(n = n(),
            min.tl = min(tl),
            max.tl = max(tl),
            mean.tl = mean(tl)
            )
print(ltlipids.tl.summary)

ltlipids.lipids.summary.age <- ltlipids %>% group_by(source, season, age_class) %>% 
  summarize(n = n(),
            min.lipids = min(avg_perc_lipid),
            max.lipids = max(avg_perc_lipid),
            mean.lipids = mean(avg_perc_lipid),
            sd.lipids = sd(avg_perc_lipid),
            se = sd.lipids/sqrt(n)
  )

ltlipids.lipids.summary <- filter(ltlipids, season != "Pre-winter") %>% group_by(source) %>% 
  summarize(n = n(),
            min.lipids = min(avg_perc_lipid),
            max.lipids = max(avg_perc_lipid),
            mean.lipids = mean(avg_perc_lipid),
            sd.lipids = sd(avg_perc_lipid),
            se = sd.lipids/sqrt(n)
  )


## FIT MODEL (lipids ~ tl * source * location * season) =========

## Remove hatchery fish first
ltlipids.filt <- ltlipids %>% filter(location != "Grand Isle Hatchery")
lm.lipids.logit <- lm(lipid.logit ~ tl * source * location * season, data = ltlipids.filt)


## CHECK ASSUMPTIONS ============================================

##  Normality
shapiro.test(ltlipids.filt$lipid.logit) 
  ## p = 0.122; accept normality

## Homogeneity of Variance by source
bartlett.test(lipid.logit ~ interaction(source, location, season), data = ltlipids.filt)
  ## p = 0.033; unequal variance!!


## BOOTSTRAPPING ================================================

# Define the number of bootstrap samples
boot_n <- 10000

## Run loop (be very patient!)
f.values.boot <- do.call(rbind, lapply(1:boot_n, function(i) {
      
  bootstrap.data <- do.call(rbind, lapply(1:nrow(ltlipids.filt), function(l) {
      ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(ltlipids.filt))
      lipid.logit <- sample(ltlipids.filt$lipid.logit, replace = T, size = 1)
      tl <- sample(ltlipids.filt$tl, replace = T, size = 1)
      source <- sample(ltlipids.filt$source, replace = T, size = 1)
      location <- sample(ltlipids.filt$location, replace = T, size = 1)
      season <- sample(ltlipids.filt$season, replace = T, size = 1)
      
      id.data <- data.frame(lipid.logit, tl, source, location, season)
  }))

  ## Calculate and extract F statistic from the bootstrapped data
  bootstrap.anova <- anova(lm(lipid.logit ~ tl * source * location * season, data = bootstrap.data))
  f.boot <- bootstrap.anova$`F value`
  f.boot.df <- data.frame("tl" = f.boot[1],
                          "source" = f.boot[2],
                          "location" = f.boot[3],
                          "season" = f.boot[4],
                          "tl:source" = f.boot[5],
                          "tl:location" = f.boot[6],
                          "source:location" = f.boot[7],
                          "tl:season" = f.boot[8],
                          "source:season" = f.boot[9],
                          "location:season" = f.boot[10],
                          "tl:source:location" = f.boot[11],
                          "tl:source:season" = f.boot[12],
                          "tl:location:season" = f.boot[13],
                          "source:location:season" = f.boot[14],
                          "tl:source:location:season" = f.boot[15])
}))

## Calculate the 95th percentile from bootstrapped distribution for each of the variables/interactions
fstat.boot.95perc <- do.call(rbind, lapply(1:ncol(f.values.boot), function(i) {
  perc <- quantile(f.values.boot[,i], probs = 0.95)
  data <- data.frame(variable = colnames(f.values.boot[i]),
                     perc95.f = perc)
}))

## Observed F statistics
obs.fstat <- data.frame(obs.f = anova(lm(lipid.logit ~ tl * source * location * season, data = ltlipids.filt))$`F value`[1:15],
                        variable = c("tl", "source", "location", "season", 
                                     "tl.source", "tl.location", "source.location", "tl.season",
                                     "source.season", "location.season", "tl.source.location",
                                     "tl.source.season", "tl.location.season", "source.location.season",
                                     "tl.source.location.season"))

## Calculate p-values from empirical distribution
p.value <- data.frame(p.value = c(
  1-ecdf(f.values.boot$tl)(filter(obs.fstat, variable == "tl")$obs.f),
  1-ecdf(f.values.boot$source)(filter(obs.fstat, variable == "source")$obs.f),
  1-ecdf(f.values.boot$location)(filter(obs.fstat, variable == "location")$obs.f),
  1-ecdf(f.values.boot$season)(filter(obs.fstat, variable == "season")$obs.f),
  1-ecdf(f.values.boot$tl.source)(filter(obs.fstat, variable == "tl.source")$obs.f),
  1-ecdf(f.values.boot$tl.location)(filter(obs.fstat, variable == "tl.location")$obs.f),
  1-ecdf(f.values.boot$source.location)(filter(obs.fstat, variable == "source.location")$obs.f),
  1-ecdf(f.values.boot$tl.season)(filter(obs.fstat, variable == "tl.season")$obs.f),
  1-ecdf(f.values.boot$source.season)(filter(obs.fstat, variable == "source.season")$obs.f),
  1-ecdf(f.values.boot$location.season)(filter(obs.fstat, variable == "location.season")$obs.f),
  1-ecdf(f.values.boot$tl.source.location)(filter(obs.fstat, variable == "tl.source.location")$obs.f),
  1-ecdf(f.values.boot$tl.source.season)(filter(obs.fstat, variable == "tl.source.season")$obs.f),
  1-ecdf(f.values.boot$tl.location.season)(filter(obs.fstat, variable == "tl.location.season")$obs.f),
  1-ecdf(f.values.boot$source.location.season)(filter(obs.fstat, variable == "source.location.season")$obs.f),
  1-ecdf(f.values.boot$tl.source.location.season)(filter(obs.fstat, variable == "tl.source.location.season")$obs.f)),
  variable = c("tl", "source", "location", "season", 
               "tl.source", "tl.location", "source.location", "tl.season",
               "source.season", "location.season", "tl.source.location",
               "tl.source.season", "tl.location.season", "source.location.season",
               "tl.source.location.season"))

## Compare oberserved F to boostrapped F 
f.test <- left_join(fstat.boot.95perc, obs.fstat) %>% 
  left_join(p.value) %>% 
  mutate(test = ifelse(obs.f > perc95.f, "Sig.", "Not Sig."))
## Significant variables and interactions:
## length, source, season
## tl:source, tl:season, location:season


## VISUALIZATION ================================================
## TL:Season Interaction; TL:Source Interaction
ggplot(ltlipids, aes(x = tl, y = avg_perc_lipid)) +
  geom_text(aes(label = age_class, color = season), show.legend = FALSE, size = 5.5, alpha = 0.8) +
  geom_smooth(data = filter(ltlipids, location != "Grand Isle Hatchery"),
              aes(linetype = season, color = season), 
              method = "lm", se = FALSE, size = 1) +
  geom_smooth(data = filter(ltlipids, location != "Grand Isle Hatchery"),
              method = "lm", se = FALSE, color = "black", size = 1, show.legend = FALSE) +
  scale_x_continuous(limits = c(75, 350), breaks = seq(100, 350, 50), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 45, 7.5), expand = c(0, 0)) +
  scale_colour_manual(values = c("black", "gray80", "gray60", "gray40", "black"),
  #scale_colour_manual(values = c("black", "#2b83ba", "#d7191c", "#fdae61", "black"), 
                      labels = c("", "", "Spring", "Summer", "Autumn")) + 
  scale_linetype_manual(values = c("dotted", "dotdash", "dashed", "solid")) + 
                        #labels = c( "Spring", "Summer", "Autumn")) +
  labs(y = 'Mean % Lipid Content', x = 'Total Length (mm)') +
  theme(axis.text = element_text(size = 20), axis.line.x = element_line(), 
        axis.line.y = element_line(), 
        axis.title.x = element_text(size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(size = 22, margin = margin(0, 15, 0, 0)),
        axis.ticks.length = unit(2.5, 'mm'), 
        legend.key = element_blank(), legend.position = "none",
        strip.background = element_blank(), strip.text = element_text(size = 22),
        panel.background = element_blank(), panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.margin = unit(c(5, 7.5, 2, 2), "mm")) +
  facet_wrap(~source, scales = "free_y")

ggsave("figures/Sorrentino_et_al_Fig2_grayscale.tiff", dpi = 300, width = 12, height = 7)
ggsave("figures/Sorrentino_et_al_Fig2_LowRes_grayscale.tiff", dpi = 100, width = 12, height = 7)



ggplot(filter(ltlipids, season != "Pre-winter"), aes(x = tl, y = avg_perc_lipid, linetype = season, color = season)) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_x_continuous(limits = c(75, 350), breaks = seq(100, 350, 50), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 45, 7.5), expand = c(0, 0)) +
  scale_colour_manual(values = c("gray80", "gray60", "gray40"), 
  #scale_colour_manual(values = c("#2b83ba", "#d7191c", "#fdae61"), 
                      labels = c("Spring", "Summer", "Autumn")) + 
  scale_linetype_manual(values = c("dotted", "dotdash", "dashed", "solid")) + 
  labs(y = 'Mean % Lipid Content', x = 'Total Length (mm)') +
  theme(axis.text = element_text(size = 20), axis.line.x = element_line(), 
        axis.line.y = element_line(), 
        axis.title.x = element_text(size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(size = 22, margin = margin(0, 15, 0, 0)),
        axis.ticks.length = unit(2.5, 'mm'), 
        legend.text = element_text(size = 16), legend.title = element_blank(),
        legend.key.width = unit(1.2, 'cm'), legend.key.height = unit(0.8, 'cm'), 
        legend.key = element_blank(), legend.position = c(0.655, 0.92),
        legend.spacing.x = unit(0.3, 'cm'),
        strip.background = element_blank(), strip.text = element_text(size = 22),
        panel.background = element_blank(), panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.margin = unit(c(5, 7.5, 2, 2), "mm")) +
  facet_wrap(~source, scales = "free_y")

ggsave("figures/Sorrentino_et_al_Fig2_legend_grayscale.tiff", dpi = 300, width = 12, height = 7)


## Location:Season Interaction

ltlipids.lp.summary <- ltlipids.filt %>% group_by(source, season, location, source) %>% 
  summarize(mean.lp = mean(avg_perc_lipid),
            sd.lp = sd(avg_perc_lipid),
            se.lp = sd.lp/sqrt(n())) %>% ungroup() %>% 
  group_by(season) %>% 
  mutate(width = 0.15 * n())

ggplot(ltlipids.lp.summary, aes(x = season, y = mean.lp, color = location, shape = location, group = location, width = width)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = mean.lp - se.lp, ymax = mean.lp + se.lp), 
                position = position_dodge(width = 0.3), size = 1) +
  geom_line(size = 1, position = position_dodge(width = 0.2)) +
  scale_y_continuous(limits = c(10, 30), breaks = seq(10, 30, 5), expand = c(0, 0)) +
  #scale_colour_manual("", values = c("black", "gray70", "gray40")) + 
  scale_colour_manual("", values = c("#2b83ba", "#d7191c", "#fdae61")) + 
  scale_shape_manual("", values = c(15, 16, 17)) +
  labs(y = 'Mean % Lipid Content', x = 'Season', color = "Location") +
  theme(axis.text = element_text(size = 20), axis.line.x = element_line(), 
        axis.line.y = element_line(), 
        axis.title.x = element_text(size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(size = 22, margin = margin(0, 15, 0, 0)),
        axis.ticks.length = unit(2.5, 'mm'), 
        legend.text = element_text(size = 16), legend.title = element_blank(),
        legend.key.width = unit(1.2, 'cm'), legend.key.height = unit(0.8, 'cm'), 
        legend.key = element_blank(), legend.position = c(0.09, 0.90),
        legend.spacing.x = unit(0.3, 'cm'),
        strip.background = element_blank(), strip.text = element_text(size = 22),
        panel.background = element_blank(), panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.margin = unit(c(5, 7.5, 2, 2), "mm")) +
  facet_wrap(~source, scales = "free_y")

ggsave("figures/Sorrentino_et_al_Fig3.tiff", dpi = 300, width = 12, height = 7)
ggsave("figures/Sorrentino_et_al_Fig3_LowRes.tiff", dpi = 100, width = 12, height = 7)


# Complete analysis of pXRF data ----

library(tidyverse)
library(tidymodels)
library(tidytext)
library(patchwork)
`%notin%` <- Negate(`%in%`)

pXRF_light_raw <- read_csv("data/pXRF_light.csv")
pXRF_heavy_raw <- read_csv("data/pXRF_heavy.csv")

## Normalise ----

normalise <- function (value, against) {return(value/against)}

### Normalise against silicon for light elements.
pXRF_light <- pXRF_light_raw %>% 
  mutate_at(vars(-id:-place), 
            ~normalise(., pXRF_light_raw$Si_K12)) %>% 
  select(-Si_K12)

rm(pXRF_light_raw)
  
### Normalise against Rhodium for light elements.
pXRF_heavy <- pXRF_heavy_raw %>% 
  mutate_at(vars(-id:-place), 
            ~normalise(., pXRF_heavy_raw$Rh_K12)) %>% 
  select(-Rh_K12)

rm(pXRF_heavy_raw)

##Plotting setup ---- 
theme_set(theme_bw())
colours <- c('BAT'='red',
             'GT'='magenta',
             'ZT'='midnightblue',
             'ZW'='blue',
             'BEL' = 'lightsalmon1',
             'BN' = 'lightsalmon2',
             'CM' = 'lightsalmon3',
             'LJ' = 'lightsalmon4',
             'RP' = 'brown',
             'TR' = 'gold',
             'UNID' = 'mediumseagreen') #set plotting colours

shapes <- c("Colonial" = 15,
             "Dutch" = 16,
             "Trial" = 18,
             "Unknown" = 3
             ) #set plotting colours

## PCA ----

### Function to compute variance explained
var_ex <- function(pca, pc) {
  sdev <- pca$steps[[2]]$res$sdev
  var <- round(sdev^2 / sum(sdev^2),3)*100
  paste0('PC',pc,': ', var[pc], '%')
}

### Function to compute variance explained
id_vars <- c('id', 'regno', 'site', 'class', 'construction', 'year', 'place')

### Function to compute variance explained
pca_light <- pXRF_light %>% recipe(~.) %>%
  update_role(all_of(id_vars), new_role = "id") %>%
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors()) %>% prep()

### Function to compute variance explained
pca_heavy <- pXRF_heavy %>%
  filter(site != 'TR') %>% 
  recipe(~.) %>%
  update_role(all_of(id_vars), new_role = "id") %>%
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors()) %>% prep()

#Plot both PCAs
pca_plot_light <- 
  juice(pca_light) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(color = site, shape=class), alpha = 0.7, size = 2) +
  scale_color_manual(values = colours) +
  scale_shape_manual(values = shapes) +
  labs(title = "PCA Light Elements",
       colour = "Site",
       shape = "Period",
       x=var_ex(pca_light, 1),
       y=var_ex(pca_light, 2)) +
  theme(legend.position = "none")

pca_plot_heavy <- 
  juice(pca_heavy) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(color = site, shape=class), alpha = 0.7, size = 2) +
  scale_color_manual(values = colours) +
  scale_shape_manual(values = shapes) +
  labs(title = "PCA Heavy Elements",
       colour = "Site",
       shape = "Period",
       x=var_ex(pca_heavy, 1),
       y=var_ex(pca_heavy, 2))

pcas <- pca_plot_light + pca_plot_heavy +
  plot_layout(guides = 'collect')

pcas

ggsave("output/pcas.png",pcas, width=9, height=4.5)

## Individual elements ----

#Summaries of all elements
pivot_longer(pXRF_light, Al_K12:Zn_K12) %>% 
  ggplot(aes(site, value)) + geom_jitter(aes(colour=site), position=position_jitter(0.2), size=1, alpha=0.7) +
  scale_color_manual(values = colours) + facet_wrap(~name, scales = 'free')
ggsave("output/elements_light.pdf", width = 22, height = 16)

pivot_longer(pXRF_heavy, Ba_K12:Zr_L1) %>% 
  ggplot(aes(site, value)) + geom_jitter(aes(colour=site), position=position_jitter(0.2), size=1, alpha=0.7) +
  scale_color_manual(values = colours) + facet_wrap(~name, scales = 'free')
ggsave("output/elements_heavy.pdf", width = 22, height = 16)

# Elements are are highly correlated, so we should reduce dimensions with PCA

pca_heavy %>%
  prep() %>% 
  tidy(2) %>% 
  filter(component %in% paste0("PC", 1:2)) %>%
  group_by(component) %>%
  top_n(7, abs(value)) %>%
  ungroup() %>%
  mutate(terms = reorder_within(terms, abs(value), component)) %>%
  ggplot(aes(value, terms, fill = value > 0)) +
  geom_col() +
  facet_wrap(~component, scales = "free_y") +
  scale_y_reordered() +
  labs(
    x = "Contribution to PC",
    y = NULL)  +
  theme(legend.position = "none")

ggsave("output/pca_heavy_loadings.png", width=4, height=2)

# Dutch wrecks have more Ba, Rb, K, Mn, V, As, and less Y.

rbk <- pXRF_heavy %>% filter(site != "TR") %>% 
ggplot(aes(Rb_K12, K_K12)) + geom_jitter(aes(colour=site), position=position_jitter(0.2), size=3, alpha=0.7) +
  scale_color_manual(values = colours) + scale_y_log10() + scale_x_log10() +
  labs(title= "Rb and K", subtitle = paste0('Correlation: ', round(cor(pXRF_heavy$Rb_K12, pXRF_heavy$K_K12),2))) 

k_year <- pXRF_heavy %>% filter(site != "TR") %>% 
  ggplot(aes(reorder(site, year), K_K12)) + geom_jitter(aes(colour=site), position=position_jitter(0.2), size=3, alpha=0.7) +
  scale_color_manual(values = colours) + scale_y_log10() + labs(x = 'Site (ordered by wrecking year, youngest to the left)', title="K vs Year Wrecked", subtitle = "Similar trend in Ba, Rb, K, Mn, V & As ")

rbk + k_year + plot_layout(guides = 'collect')
ggsave("output/rbk_and_kyear.png", width=9, height=4.5)

#Periods can be reliably separated on PC1, potassium here as an example.

### Iron ----

pXRF_heavy %>% arrange(-Fe_K12) %>% select(regno, Fe_K12) %>% top_n(10)
#remove ZT3380, which is an order of a magnitude higher than the others, and we noted it was on an iron stain.

pXRF_heavy %>% filter(regno != "ZT3380") %>% filter(site != "TR") %>% 
ggplot(aes(construction, Fe_K12)) + geom_jitter(aes(colour=site), position=position_jitter(0.2), size=1.5, alpha=0.7) +
  scale_color_manual(values = colours) + scale_y_log10()
ggsave("output/iron_construction.png", width=5, height=3.5)

# It looks like site is a factor here, with iron and composite wrecks all having higher levels of iron. However, the ranger overlaps.
# Figure 4. Iron levels tend to be higher in iron and composite ships, but are within variation for wooden ones. One ZT removed with an Fe_K12 count of 30.9, because of visible iron staining. However, this is based on a single artefact from an iron wreck, and so needs much more work.
## Trial PCA ----

#Compute PCA for heavy elements with Trial, but without ZW, ZT or UNIDs
pca_trial <- pXRF_heavy %>%
  filter(site %notin% c('ZW', "UNID", "ZT")) %>% 
  recipe(~.) %>%
  update_role(all_of(id_vars), new_role = "id") %>%
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors()) %>% prep()

pca_plot_trial <- 
  juice(pca_trial) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(color = site, shape=class), alpha = 0.7, size = 3) +
  scale_color_manual(values = colours) +
  scale_shape_manual(values = shapes) +
  labs(title = "PCA - Heavy Elements with Trial",
       colour = "Site",
       shape = "Period",
       x=var_ex(pca_trial, 1),
       y=var_ex(pca_trial, 2))
pca_plot_trial

pc_year <- juice(pca_trial) %>% ggplot(aes(PC1, year)) +
  geom_point(aes(colour=site, shape=class), size=3, alpha=0.7) +
  labs(colour = "Site", shape = "Period") +
  scale_color_manual(values = colours) + scale_shape_manual(values = shapes)

pca_plot_trial + pc_year + plot_layout(guides = 'collect')

ggsave("output/pca_trial.png", width=9, height=4.5)
#Figure 5. PCA with Trial included, and ZW/ZT removed. Trial is closest in time to Batavia but is English. The fact it plots with the colonial wrecks indicates PC1 reflects manufacture/origin rather than period or time underwater

#Complete analysis of pXRF data

library(tidyverse)
library(tidymodels)
library(patchwork)

pXRF_light_raw <- read_csv("data/pXRF_light.csv")
pXRF_heavy_raw <- read_csv("data/pXRF_heavy.csv")

## Normalise

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

## Plotting setup
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

## PCA

#Function to compute variance explained
var_ex <- function(pca, pc) {
  sdev <- pca$steps[[2]]$res$sdev
  var <- round(sdev^2 / sum(sdev^2),3)*100
  paste0('PC',pc,': ', var[pc], '%')
}

#ID Variables
id_vars <- c('id', 'regno', 'site', 'class', 'construction', 'year', 'place')

#Compute PCA for light elements
pca_light <- pXRF_light %>% recipe(~.) %>%
  update_role(all_of(id_vars), new_role = "id") %>%
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors()) %>% prep()

#Compute PCA for heavy elements
pca_heavy <- pXRF_heavy %>% recipe(~.) %>%
  update_role(all_of(id_vars), new_role = "id") %>%
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors()) %>% prep()

#Plot both PCAs
pca_plot_light <- 
  juice(pca_light) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(color = site, shape=class), alpha = 0.7, size = 2) +
  scale_color_manual(values = colours) +
  labs(title = "Light Elements",
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
  labs(title = "Heavy Elements",
       colour = "Site",
       shape = "Period",
       x=var_ex(pca_heavy, 1),
       y=var_ex(pca_heavy, 2))

pcas <- pca_plot_light + pca_plot_heavy +
  plot_layout(guides = 'collect')

ggsave("output/pcas.png",pcas, width=9, height=4.5)

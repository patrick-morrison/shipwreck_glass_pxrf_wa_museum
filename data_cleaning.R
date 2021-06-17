#Data cleaning
library(tidyverse)
library(readxl)
`%notin%` <- Negate(`%in%`)

light_elements_raw <- read_excel("data/raw/tracer_elements_light_glass_july2020_final.xlsx", 
                           sheet = "Points")

light_meta <- read_excel("data/raw/metadata_light_glass_july2020_final.xlsx")

light_elements <- light_elements_raw %>%
  rename_all(~str_replace_all(., "\\s+", "_")) %>%
  mutate(id = as.integer(str_remove(`...1`, '-Spectrometer Mode@150720_114558'))) %>%
  right_join(light_meta) %>% filter(site %notin% c('Control', 'Exclude')) %>% 
  mutate(class = as.factor(class),
         site = as.factor(site),
         regno = str_remove(`regno`, stringr::fixed("(Blue)")),
         year = case_when(
           site == 'BAT' ~ 1629,
           site == 'GT' ~ 1656,
           site == 'CM' ~ 1830,
           site == 'LJ' ~ 1921,
           site == 'ZT' ~ 1712,
           site == 'ZW' ~ 1727,
           site == 'BEL' ~ 1824,
           site == 'BN' ~ 1888,
           site == 'RP' ~ 1811),
         
         construction = as.factor(
           case_when(
             site %in% c('BAT', 'GT', 'ZT', 'ZW', 'LJ', "RP") ~ "Wooden",
             site %in% c('CM', 'BEL') ~ "Composite",
             site %in% c("BN") ~ "Iron",
             site %in% c('UNID') ~ "Unknown"))) %>% 
  select(id, regno, site, class, construction, year, place, everything()) %>% distinct() %>% 
  filter(id %notin% c(120L, 121L, 115L, 47L, 118L))

write_csv(light_elements, "data/pXRF_light.csv")

heavy_elements_raw <- read_excel("data/raw/tracer_elements_heavy_glass_may2021_final_trial.xlsx", 
                           sheet = "Points")

heavy_elements <- heavy_elements_raw %>%
  rename_all(~str_replace_all(., "\\s+", "_")) %>%
  mutate(...1 = str_remove(...1, stringr::fixed("(Blue)")),
         ...1 = str_remove(...1, stringr::fixed("HeYellowFil1240_")),
         ...1 = str_remove(...1, stringr::fixed("@150720_101609")),
         ...1 = str_remove(...1, stringr::fixed("@081220")),
         ...1 = str_remove(...1, stringr::fixed("@290421")),
         ...1 = str_remove(...1, stringr::fixed("2020")),
         ...1 = str_replace(...1, stringr::fixed("__"), '_')
  ) %>% 
  separate(...1, into=c('regno', 'place','id'), convert = TRUE, sep=c("_")) %>% 
  
  mutate(
    site = as.factor(
      case_when(
        str_detect(regno, 'BAT') ~ "BAT",
        str_detect(regno, 'GT') ~ "GT",
        str_detect(regno, 'CM') ~ "CM",
        str_detect(regno, 'LJ') ~ "LJ",
        str_detect(regno, 'ZT') ~ "ZT",
        str_detect(regno, 'ZW') ~ "ZW",
        str_detect(regno, 'TR') ~ "TR",
        str_detect(regno, 'BEL') ~ "BEL",
        str_detect(regno, 'BN') ~ "BN",
        str_detect(regno, 'RP') ~ "RP",
        str_detect(regno, 'UNID') ~ "UNID",
      )), 
    class = as.factor(
      case_when(
        site %in% c('BAT', 'GT', 'ZT', 'ZW') ~ "Dutch",
        site %in% c('CM', 'LJ', 'BEL', "BN", "RP", "TR") ~ "Colonial",
        site %in% c('UNID') ~ "Unknown",
      )),
    year = case_when(
      site == 'BAT' ~ 1629,
      site == 'GT' ~ 1656,
      site == 'CM' ~ 1830,
      site == 'LJ' ~ 1921,
      site == 'ZT' ~ 1712,
      site == 'ZW' ~ 1727,
      site == 'BEL' ~ 1824,
      site == 'BN' ~ 1888,
      site == 'TR' ~ 1622,
      site == 'RP' ~ 1811),
    
    construction = as.factor(
      case_when(
        site %in% c('BAT', 'GT', 'ZT', 'ZW', 'LJ', "RP", "TR") ~ "Wooden",
        site %in% c('CM', 'BEL') ~ "Composite",
        site %in% c("BN") ~ "Iron",
        site %in% c('UNID') ~ "Unknown",
      ))) %>% 
  filter(regno %notin% c("trueblank07", "SiO2")) %>% 
  select(id, regno, site, class, construction, year, place, everything()) %>% distinct()

write_csv(heavy_elements, "data/pXRF_heavy.csv")
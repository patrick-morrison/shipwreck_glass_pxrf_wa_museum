library(geojsonsf)
library(ggspatial)
library(ggrepel)
library(ozmaps)
library(tidyverse)
theme_set(theme_minimal())

sites <- geojson_sf("data/shipwreck_sites.geojson")
wa <- ozmap_data(data = "states") %>% filter(NAME=='Western Australia')

ggplot(sites) +
  geom_sf(data=wa, alpha=1) +
  geom_sf() +
  geom_text_repel(stat = "sf_coordinates", aes(label=name, geometry = geometry), 
                   size=3.2, fontface="italic",
                   nudge_x=c(0, 3.2, 0, 4,   4,  4,   4, 3, 3),
                   nudge_y=c(-1.3,   0, 0.2, 0, -.3, .2, -.3, -.4, -.1)) +
  annotation_scale(location = "br", height = unit(0.06, "in")) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0, "in"), pad_y = unit(.1, "in"),
                         style = north_arrow_nautical) +
  labs(x=NULL, y=NULL)
ggsave("output/glass_sites_map.png", width = 3.2, height = 4.5, bg='white')  
#Source: WA Museum Shipwrecks Database. Basemap: Australian Bureau of Statistics via the ozmaps packages in R."
 
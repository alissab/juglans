
require(rgdal)
require(raster)
require(dplyr)
require(rgeos)
require(ggplot2)

# projection you want to use
proj_out <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5
+lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# WGS84 (lat/long), coordinate format of data downloads
proj_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
+towgs84=0,0,0"

# need North America shapefile
na_shp <- readOGR("shapefiles/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_out)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
# glacier shapefile (start with 1k YBP to present)
shp_1k <- readOGR("ice/1ka_0.91cal_Dalton_et_al_2020_QSR.shp", "1ka_0.91cal_Dalton_et_al_2020_QSR")
shp_1k <- sp::spTransform(shp_1k, proj_out)

# read in pollen data, select observations from 1k YBP - present
pollen <- readRDS("juglans_data.RData")
pollen_1k <- pollen %>% filter(age <= 1000) %>% select(-juglans.nigra, -cut)
pollen_1k <- pollen_1k %>% filter(juglans.cinerea >0, !is.na(juglans.cinerea))

# sum site counts together (should have 1 count per site)
site_1k <- pollen_1k[!duplicated(pollen_1k[, c("site.name")]), c("site.name",'x','y')]

count_1k <- pollen_1k %>% select(site.name, juglans.cinerea) %>% 
  group_by(site.name) %>% 
  summarise(juglans.cinerea = sum(juglans.cinerea))
count_1k <- ungroup(count_1k)

dat_1k <- left_join(count_1k, site_1k, by = "site.name")

# plot glaciers and pollen sites together (1k YBP - current)
ggplot() +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
            alpha = 0.5) +
  geom_path(data = shp_1k, aes(x = long, y = lat, group = group), 
            size = 2, color = 'blue') +
  geom_point(data = dat_1k, aes(x = x, y = y), 
             color = 'red', size = 4, alpha = 0.5) +
  scale_y_continuous(limits = c(-1500000, 4000000)) +
  scale_x_continuous(limits = c(-1000000, 3000000)) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()



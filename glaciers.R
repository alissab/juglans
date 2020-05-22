
require(rgdal)
require(raster)
require(dplyr)
require(rgeos)
require(ggplot2)
require(gridExtra)

# require(geosphere)
# require(sp)

# read in pollen data
pollen <- readRDS("juglans_data.RData")

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

# START WITH 1K YBP
# glacier shapefile (start with 1k YBP to present)
shp_1k <- readOGR("ice/1ka_0.91cal_Dalton_et_al_2020_QSR.shp", "1ka_0.91cal_Dalton_et_al_2020_QSR")
shp_1k <- sp::spTransform(shp_1k, proj_out)

# select pollen observations from 1k YBP - present
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



# PERFORM PREVIOUS STEPS ON ALL DATA IN 1000-YEAR TIME CHUNKS

# cut data into time chunks (if it doesn't already have a 'cut' column)
# dat.cut <- cut(pollen$age, breaks = c(min(pollen$age, na.rm = TRUE), seq(1000, 21000, by = 1000)),
#                labels = FALSE)
# pollen$cut <- dat.cut
# pollen <- pollen[!is.na(pollen$cut),]

# separate dataframe into list of dataframes, one for each time chunk
# only use j. cinerea
time <- split(pollen, f = pollen$cut)
cinerea_list <- lapply(time, "[", c(4, 7:9))

# remove sites that contain NA or 0 pollen counts for each time chunk
  for(i in 1:21){
    cinerea_list[[i]] <- subset(cinerea_list[[i]], 
                                 !is.na(cinerea_list[[i]][,"juglans.cinerea"]) & 
                                  cinerea_list[[i]][,"juglans.cinerea"] != 0)
  }

# sum counts when there are multiple counts at a site during each time period
# (need one count per site per time period)
agg_list <- list()

for(i in 1:21){
agg_list[[i]] <- cinerea_list[[i]] %>% group_by(x, y) %>% 
  summarise(juglans.cinerea = sum(juglans.cinerea))
agg_list[[i]] <- as.data.frame(agg_list[[i]])
}


# PLOTTING LOOPS OVER TIME INTERVALS, CIRCLES FILLED BY RELATIVE ABUNDANCE
plot_list <- list(list())

  for(i in 1:length(agg_list)){
    plot_list[[i]] <- ggplot(data = agg_list[[i]]) +
      geom_point(aes_string(x = 'x', y = 'y', fill = 'juglans.cinerea'), 
                 size = 4, pch = 21, alpha = 0.5, color = 'red') +
      geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
                alpha = 0.5) +
      labs(subtitle = paste0(i, "000")) +
      scale_y_continuous(limits = c(-1500000, 4000000)) +
      scale_x_continuous(limits = c(-1000000, 3000000)) +
      theme_classic() +
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            line = element_blank(),
            legend.position = 'none') +
      coord_fixed()
  }


# READ IN GLACIER SHAPEFILES AND ADD THEM TO PLOTS BY TIME CHUNK
for(i in 1:21){
ice <- readOGR("ice/1ka_0.91cal_Dalton_et_al_2020_QSR.shp", 
                  "1ka_0.91cal_Dalton_et_al_2020_QSR")
ice <- sp::spTransform(shp_1k, proj_out)
}


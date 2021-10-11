
setwd('C:/Users/abrow/Documents/juglans')
require(rgdal)
require(raster)
require(dplyr)
require(rgeos)
require(geosphere)
require(ggplot2)
require(gridExtra)
require(ggspatial)
require(cowplot)

# READ IN POLLEN DATA
pollen <- readRDS("butternut_pollen_fixed_dates.RDS")
pollen <- pollen %>% filter(long > -110)

# READ IN AND ORGANIZE GLACIER SHAPEFILES
shp_21k <- readOGR("ice/17.5ka_21.1cal_Dalton_et_al_2020_QSR.shp", "17.5ka_21.1cal_Dalton_et_al_2020_QSR")
shp_20k <- readOGR("ice/17ka_20.5cal_Dalton_et_al_2020_QSR.shp", "17ka_20.5cal_Dalton_et_al_2020_QSR")
shp_19k <- readOGR("ice/16ka_19.3cal_Dalton_et_al_2020_QSR.shp", "16ka_19.3cal_Dalton_et_al_2020_QSR")
shp_18k <- readOGR("ice/15ka_18.0cal_Dalton_et_al_2020_QSR.shp", "15ka_18.0cal_Dalton_et_al_2020_QSR")
shp_17k <- readOGR("ice/14ka_16.8cal_Dalton_et_al_2020_QSR.shp", "14ka_16.8cal_Dalton_et_al_2020_QSR")
shp_16k <- readOGR("ice/13.5ka_16.1cal_Dalton_et_al_2020_QSR.shp", "13.5ka_16.1cal_Dalton_et_al_2020_QSR")
shp_15k <- readOGR("ice/12.5ka_14.9cal_Dalton_et_al_2020_QSR.shp", "12.5ka_14.9cal_Dalton_et_al_2020_QSR")
shp_14k <- readOGR("ice/11.5ka_13.5cal_Dalton_et_al_2020_QSR.shp", "11.5ka_13.5cal_Dalton_et_al_2020_QSR")
shp_13k <- readOGR("ice/11ka_12.8cal_Dalton_et_al_2020_QSR.shp", "11ka_12.8cal_Dalton_et_al_2020_QSR")
shp_12k <- readOGR("ice/10.5ka_12.1cal_Dalton_et_al_2020_QSR.shp", "10.5ka_12.1cal_Dalton_et_al_2020_QSR")
shp_11k <- readOGR("ice/9.6ka_11cal_Dalton_et_al_2020_QSR.shp", "9.6ka_11cal_Dalton_et_al_2020_QSR")
shp_10k <- readOGR("ice/8.5ka_9.6cal_Dalton_et_al_2020_QSR.shp", "8.5ka_9.6cal_Dalton_et_al_2020_QSR")
shp_9k <- readOGR("ice/8ka_9cal_Dalton_et_al_2020_QSR.shp", "8ka_9cal_Dalton_et_al_2020_QSR")
shp_8k <- readOGR("ice/7ka_7.9cal_Dalton_et_al_2020_QSR.shp", "7ka_7.9cal_Dalton_et_al_2020_QSR")
shp_7k <- readOGR("ice/6ka_6.8cal_Dalton_et_al_2020_QSR.shp", "6ka_6.8cal_Dalton_et_al_2020_QSR")
shp_6k <- readOGR("ice/5.5ka_6.3cal_Dalton_et_al_2020_QSR.shp", "5.5ka_6.3cal_Dalton_et_al_2020_QSR")
shp_5k <- readOGR("ice/5ka_5.71cal_Dalton_et_al_2020_QSR.shp", "5ka_5.71cal_Dalton_et_al_2020_QSR")
shp_4k <- readOGR("ice/4ka_4.5cal_Dalton_et_al_2020_QSR.shp", "4ka_4.5cal_Dalton_et_al_2020_QSR")
shp_3k <- readOGR("ice/3ka_3.2cal_Dalton_et_al_2020_QSR.shp", "3ka_3.2cal_Dalton_et_al_2020_QSR")
shp_2k <- readOGR("ice/2ka_2cal_Dalton_et_al_2020_QSR.shp", "2ka_2cal_Dalton_et_al_2020_QSR")
shp_1k <- readOGR("ice/1ka_0.91cal_Dalton_et_al_2020_QSR.shp", "1ka_0.91cal_Dalton_et_al_2020_QSR")

# fix projections as needed (should be unprojected lat/long)
shp_2k <- spTransform(shp_2k, proj4string(shp_1k))
shp_3k <- spTransform(shp_3k, proj4string(shp_1k))
shp_4k <- spTransform(shp_4k, proj4string(shp_1k))

# in case you need proj4string for lat/long:
# spTransform(XXX, "+init=epsg:4326")

nchunks <- 21

# place all glacier shapefiles in a list
glacier_list <- list()
for(i in 1:nchunks){
  glacier_list[[i]] <- get(x = paste0('shp_',i,'k'))
}
saveRDS(glacier_list, 'glacier_shapefiles_21-1k.RDS')

# FIND AVERAGE DISTANCE FROM SITES TO GLACIER POLYGON FOR EACH TIME CHUNK
avg <- as.vector(rep(NA, nchunks))
for(i in 1:nchunks){
  hi.age <- i*1000 + 500
  low.age <- i*1000 - 500
  pol <- pollen %>% filter(age <= hi.age, age >= low.age)  # select pollen data within time chunk
  pol <- pol[!duplicated(pol[,c('long','lat')]), c('long','lat')]
  dist <- dist2Line(line = glacier_list[[i]], p = pol)  # distance to glacier for each site
  avg[i] <- mean(dist[,1])  # calculate average distance to glacier for each time chunk
  cat(paste0("iteration ", i, "\n"))
}

# CALCULATE APPROXIMATE BIOTIC VELOCITY
bv <- data.frame(dist = avg, ybp = seq(1000, 21000, by = 1000))
bv$bv <- NA
for(i in 1:nrow(bv)){
  if(i<21){
    bv$bv[i] <- ( bv[i+1, 1] - bv[i, 1] ) / 1000}
  else{
    bv$bv[i] <- NA
  }
}

# convert meters to km
bv$dist <- bv$dist/1000
# saveRDS(bv, 'biot_vel.RDS')
# bv <- readRDS('biot_vel.RDS')

# plot bv over time
png('figs/dist_to_glacier.png', width = 12, height = 10, units = 'cm', res = 300)
par(mar = c(5,5,2,2))
plot(bv$ybp, bv$dist, type = 'l', xlab = 'Years before present', 
     ylab = 'Distance from glacier (km)')
abline(h = mean(bv$dist), lty = 'dotted')
dev.off()

png('figs/biot_vel.png', width = 12, height = 10, units = 'cm', res = 300)
par(mar = c(5,5,2,2))
plot(bv$ybp, bv$bv, type = 'l', xlab = 'Years before present', 
     ylab = 'Biotic velocity (m/yr)')
abline(h = mean(bv$bv, na.rm = TRUE), lty = 'dotted')
dev.off()


# FIND CENTROID OF POLLEN OCCURRENCES (TAKE MEAN OF X, Y COORDS FOR EACH TIME CHUNK)
# WHERE 1 TIME CHUNK = 1000 YEARS
projection <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
nchunks <- 22
avg <- data.frame(x = rep(NA, nchunks), y = rep(NA, nchunks), 
                  dist = rep(NA, nchunks),
                  chunk = rep(NA, nchunks))
for(i in 1:nchunks){
  hi.age <- (i-1)*1000 + 500
  low.age <- (i-1)*1000 - 500
  pol <- pollen %>% filter(age < hi.age, age >= low.age)  
  avg$x[i] <- mean(pol$x)
  avg$y[i] <- mean(pol$y)
  avg$chunk[i] <- as.character(i-1)  
}

# calculate distance between centroids i+1 and i
coords <- avg[, c('x','y')]
sp::coordinates(coords) <- ~x + y
sp::proj4string(coords) <- sp::CRS(projection)

avg$dist <- NA
for(i in 1:(nchunks-1)){
  avg$dist[i+1] <- pointDistance(coords[i+1], coords[i], lonlat = FALSE)
}

# calculate BV by dividing dist by 1000years and plot
avg$bv <- avg$dist/1000
avg$chunk <- as.integer(avg$chunk)

ggplot(data = avg, aes(x = chunk, y = bv)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(limits = c(1,21), breaks = seq(1,21, by = 1)) +
  scale_y_continuous(limits = c(0,500)) +
  labs(x = '1,000 years before present', y = 'Centroid velocity (m/yr)') +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 3))




# SAME THING AS ABOVE, BUT WITH ORDONEZ TIME CHUNKS
# fast climate change: 16-14, 14-12, 12-10
# slow climate change: 10-7, 7-4, 4-1
nchunks <- 7
avg_ord <- data.frame(x = rep(NA, nchunks), y = rep(NA, nchunks), 
                  dist = rep(NA, nchunks), chunk = rep(NA, nchunks),
                  age = c(1000, 4000, 7000, 10000, 12000, 14000, 16000),
                  bv = rep(NA, nchunks))

for(i in 1:nchunks){
  hi_age <- avg_ord$age[i] + 500
  lo_age <- avg_ord$age[i] - 500
  pol <- pollen %>% filter(age < hi_age, age >= lo_age)  
  avg_ord$x[i] <- mean(pol$x)
  avg_ord$y[i] <- mean(pol$y)
}

# calculate distance between centroids i+1 and i
coords <- avg_ord[ , c('x','y')]
sp::coordinates(coords) <- ~x + y
sp::proj4string(coords) <- sp::CRS(projection)

for(i in 1:(nchunks-1)){
  avg_ord$dist[i] <- pointDistance(coords[i+1], coords[i], lonlat = FALSE)
  avg_ord$chunk[i] <- paste0(avg_ord$age[i], "-", avg_ord$age[i+1], " ybp")
}

# calculate BV by dividing dist by years and plot
avg_ord$bv[1:3] <- avg_ord$dist/3000
avg_ord$bv[4:6] <- avg_ord$dist/2000
print(avg_ord[,c('chunk','bv')])
#         chunk        bv
# 1   1000-4000 ybp  13.17870
# 2   4000-7000 ybp  12.78266
# 3  7000-10000 ybp 116.49299
# 4 10000-12000 ybp  19.76805
# 5 12000-14000 ybp  19.17398
# 6 14000-16000 ybp 174.73949




# USE ADAM SMITH'S BIOTIC VELOCITY FUNCTION
# first need to convert HSM maps to raster stacks
# NOTE THAT BELOW CODE DOES NOT MASK OUT GLACIERS BEFORE CALCULATING BV
library(raster)
library(enmSdm)

# read in habitat suitability rasters
lig <- raster('shapefiles/HSM_rasters_dif_timepoints/lig_hsm.tif')
lgm <- raster('shapefiles/HSM_rasters_dif_timepoints/lgm_hsm.tif')
hs <- raster('shapefiles/HSM_rasters_dif_timepoints/hs_hsm.tif')
ba <- raster('shapefiles/HSM_rasters_dif_timepoints/ba_hsm.tif')
yds <- raster('shapefiles/HSM_rasters_dif_timepoints/yds_hsm.tif')
e_hol <- raster('shapefiles/HSM_rasters_dif_timepoints/eh_hsm.tif')
mid_hol <- raster('shapefiles/HSM_rasters_dif_timepoints/mid_hol_hsm.tif')
late_hol <- raster('shapefiles/HSM_rasters_dif_timepoints/lh_hsm.tif')
hsm_stack <- stack(list(lig, lgm, hs, ba, yds, e_hol, mid_hol, late_hol))
times <- c(-130, -22, mean(c(-17, -14.7)), mean(c(-14.7, -12.9)), mean(c(-12.9, -11.7)), 
           mean(c(-11.7, -8.326)), mean(c(-8.326, -4.2)), mean(c(-4.2, -0.0003)))
bv <- bioticVelocity(hsm_stack, times = times)
bv$centroidVelocity <- bv$centroidVelocity / 1000
bv$period <- c('LGM','HS','BA','YDS','EH','mid-Hol','late-Hol')
bv$timeTo_plot <- -1 * bv$timeTo

ggplot(bv, aes(x = timeTo_plot, y = centroidVelocity)) + 
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0, 225), name = 'Centroid velocity (m/yr)') +
  scale_x_continuous(name = '1,000 years before present', breaks = seq(0, 22, by = 1)) +
  geom_segment(aes(x = 22.5, y = 25, xend = 21.5, yend = 25), size = 1, col = 'red') +
  geom_segment(aes(x = 17, y = 25, xend = 14.7, yend = 25), size = 1, col = 'orange') +
  geom_segment(aes(x = 14.7, y = 25, xend = 12.9, yend = 25), size = 1, col = 'yellow') +
  geom_segment(aes(x = 12.9, y = 25, xend = 11.7, yend = 25), size = 1, col = 'green') +
  geom_segment(aes(x = 11.7, y = 25, xend = 8.326, yend = 25), size = 1, col = 'blue') +
  geom_segment(aes(x = 8.326, y = 25, xend = 4.2, yend = 25), size = 1, col = 'blueviolet') +
  geom_segment(aes(x = 4.2, y = 25, xend = 0.0003, yend = 25), size = 1, col = 'violet') +
  geom_text(aes(label = period), nudge_y = -bv$centroidVelocity + 30) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 3))
  


# COMBINE HSM- AND POLLEN-CALCULATED BVS INTO ONE PLOT
hsm_pollen <- data.frame(time = c(avg[2:nrow(avg), 'chunk'],
                                  bv$timeTo_plot),
                         bv = c(avg[2:nrow(avg), 'bv'],
                                bv$centroidVelocity),
                         type = c(rep('pollen', nrow(avg)-1), 
                                  rep('hsm', nrow(bv))),
                         period = c(rep(NA, nrow(pollen)), bv$period)
                         )
ggplot(data = hsm_pollen, aes(x = time, y = bv, color = type)) + 
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_y_continuous(name = 'Centroid velocity (m/yr)', limits = c(-30, 500)) +
  scale_x_continuous(name = '1,000 years before present', breaks = seq(1, 22, by = 1)) +
  geom_segment(aes(x = 22.5, y = -5, xend = 21.5, yend = -5), size = 1, col = 'red') +
  geom_segment(aes(x = 17, y = -5, xend = 14.7, yend = -5), size = 1, col = 'orange') +
  geom_segment(aes(x = 14.7, y = -5, xend = 12.9, yend = -5), size = 1, col = 'yellow') +
  geom_segment(aes(x = 12.9, y = -5, xend = 11.7, yend = -5), size = 1, col = 'green') +
  geom_segment(aes(x = 11.7, y = -5, xend = 8.326, yend = -5), size = 1, col = 'blue') +
  geom_segment(aes(x = 8.326, y = -5, xend = 4.2, yend = -5), size = 1, col = 'blueviolet') +
  geom_segment(aes(x = 4.2, y = -5, xend = 0.0003, yend = -5), size = 1, col = 'violet') +
  
  geom_segment(aes(x = 12, xend = 10, y = 89.5, yend = 89.5), size = 1, col = 'black', alpha = 0.5) +
  geom_segment(aes(x = 11, xend = 11, y = 176, yend = 3), size = 1, col = 'black') +
  
  geom_text(aes(label = period), y = -20, color = 'black') +
  annotate('text', label = 'Ordonez', x = 11, y = 190) +
  
  labs(color = 'Data source') +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 3))





# PROJECT COORDINATES AND PLOT POLLEN OCC WITH GLACIERS
# projection you want to use
proj_out <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5
+lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# WGS84 (lat/long), coordinate format of pollen data and glaciers
# proj_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
# +towgs84=0,0,0"

# read in North America shapefile
na_shp <- readOGR("shapefiles/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_out)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("shapefiles/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj_out)

# read pollen data
pollen <- readRDS("butternut_pollen_fixed_dates.RDS")
pollen <- pollen %>% filter(long > -110)

# project glacier polygons
glacier_list <- readRDS('glacier_shapefiles_21-1k.RDS')
nchunks <- length(glacier_list)

glacier_list_pr <- list()
for(i in 1:nchunks){
  glacier_list_pr[[i]] <- sp::spTransform(glacier_list[[i]], proj_out)
}

# POLLEN SITES ALREADY PROJECTED; CREATE LIST OF POLLEN DATA WITHIN TIME CHUNKS
pollen_list <- list()
for(i in 1:nchunks){
  hi.age <- i*1000 + 500
  low.age <- i*1000 - 500
  pollen_list[[i]] <- pollen %>% 
    filter(age <= hi.age, age >= low.age) %>% 
    select(x, y, Juglans.cinerea, site.name) %>% 
    group_by(x, y, site.name) %>% 
    summarise(count = sum(Juglans.cinerea, na.rm = TRUE))
}

# which min/max x/y values should you use to make maps
miny <- which.min(lapply(pollen_list, function(x) min(x$y, na.rm=TRUE)))
minx <- which.min(lapply(pollen_list, function(x) min(x$x, na.rm=TRUE)))
maxy <- which.max(lapply(pollen_list, function(x) max(x$y, na.rm=TRUE)))
maxx <- which.max(lapply(pollen_list, function(x) max(x$x, na.rm=TRUE)))


# PLOT GLACIERS AND POLLEN SITES TOGETHER FOR EACH TIME CHUNK
# (for selected time periodS)
times <- c(1, 11:21)

plot_list <- list()
for(i in times){
plot_list[[i]] <- ggplot() +
  geom_polygon(data = glacier_list_pr[[i]], aes(x = long, y = lat, group = group),
               size = 0.5, fill = 'slategray1', color = 'slategray') +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
            alpha = 0.5) +
  geom_point(data = pollen_list[[i]], aes(x = x, y = y), 
             pch = 21, color = 'black', fill = 'red', alpha = 0.9, size = 2) +
  
  # geom_text_repel(data = pollen_list[[i]], 
  #           aes(x = x, y = y, label = site.name)) +
  # scale_y_continuous(limits = c(-1500000, glacier_list_pr[[21]]@bbox[2,2])) +
  # scale_x_continuous(limits = c(-1000000, 3000000)) +
  ggtitle(paste0(i, ",000 YBP")) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_blank(),
        legend.text = element_blank(),
        # legend.title = element_text(size = 16),
        # legend.text = element_text(size = 14),
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.margin = unit(rep(0.25, 4), 'cm')) +
  coord_fixed(xlim = c(min(pollen_list[[11]]$x) - 500, max(pollen_list[[maxx]]$x) + 500),
              ylim = c(min(pollen_list[[11]]$y) - 500, 1700000))
}

plot_list2 <- plot_list[c(1,11:21)]

# grid.arrange - prints on screen
# do.call(grid.arrange, c(plot_list2, ncol = 4))

# arrangeGrob returns grob without drawing it
# allows you to use ggsave on the entire series of plots (not just
# the most recent plot drawn, which is what grid.arrange will do I think)
# marrangeGrob returns grob across multiple pages, if you need more space

layout <- arrangeGrob(grobs = plot_list2, nrow = 4, ncol = 3)

ggsave('figs/zoomed_pollen_with_glaciers_for_pub_reversed_times.png', plot = layout, 
       width = 6.5, height = 7, units = 'in', dpi = 300)

saveRDS(layout, 'plot_layout_for_publication.RDS')

# REVERSE ORDER OF PLOT LAYOUT
layout <- readRDS('plot_layout_for_publication.RDS')
grobs <- layout$grobs
grobs <- rev(grobs)
layout <- arrangeGrob(grobs = grobs, nrow = 4, ncol = 3)




#### PLOT HABITAT SUITABILITY WITH GLACIAL COVERAGE
require(rgdal)
require(raster)
require(dplyr)
require(ggplot2)
require(cowplot)
require(ggspatial)

#### START WITH MAPS THAT INCLUDE SHORELINES AND GLACIER COVERAGE (MAKE LOOP)
#### THEN MAKE LIG, PRESENT DAY, AND ZOOMED IN LGM

# LGM THROUGH LATE HOLOCENE MAPS
# read in habitat suitability rasters
lgm <- raster('shapefiles/HSM_rasters_dif_timepoints/lgm_hsm.tif')
hs <- raster('shapefiles/HSM_rasters_dif_timepoints/hs_hsm.tif')
ba <- raster('shapefiles/HSM_rasters_dif_timepoints/ba_hsm.tif')
yds <- raster('shapefiles/HSM_rasters_dif_timepoints/yds_hsm.tif')
e_hol <- raster('shapefiles/HSM_rasters_dif_timepoints/eh_hsm.tif')
mid_hol <- raster('shapefiles/HSM_rasters_dif_timepoints/mh_hsm.tif')
late_hol <- raster('shapefiles/HSM_rasters_dif_timepoints/lh_hsm.tif')
hs_list <- list(lgm, hs, ba, yds, e_hol, mid_hol, late_hol)

# project glacier polygons
glacier_list <- readRDS('glacier_shapefiles_21-1k.RDS')
nchunks <- length(glacier_list)
proj_out <- proj4string(lgm)

glacier_list_pr <- list()
for(i in 1:nchunks){
  glacier_list_pr[[i]] <- sp::spTransform(glacier_list[[i]], proj_out)
}
glacier_list <- glacier_list_pr[c(21, 16, 13, 12, 10, 6, 2)]

# read in North America shapefiles
na_shp <- readOGR("shapefiles/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_out)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("shapefiles/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj_out)

# read shoreline shapefiles
sh_lgm <- raster('shapefiles/shorelines/bio_2_lgm.tif')
sh_hs <- raster('shapefiles/shorelines/bio_2_hs.tif')
sh_ba <- raster('shapefiles/shorelines/bio_2_ba.tif')
sh_yds <- raster('shapefiles/shorelines/bio_2_yds.tif')
sh_eh <- raster('shapefiles/shorelines/bio_2_eh.tif')
sh_mh <- raster('shapefiles/shorelines/bio_2_mh.tif')
sh_lh <- raster('shapefiles/shorelines/bio_2_lh.tif')

shore_list <- list(sh_lgm, sh_hs, sh_ba, sh_yds, sh_eh, sh_mh, sh_lh)
n_times <- length(shore_list)

# specify extent you want to plot
plot_ext <- c(0, 3200000, -1200000, 1600000)
plot_ext <- extent(plot_ext)

# specify time period for each iteration
times <- c('21,000 YBP', '17,000 - 14,700 YBP', 
           '14,700 - 12,900 YBP', '12,900 - 11,700 YBP', '11,700 - 8,326 YBP', 
           '8,326 - 4,200 YBP', '4,200 - 300 YBP')

plot_list <- list()
for(i in 1:n_times){
  df <- as.data.frame(hs_list[[i]], xy = TRUE)
  names(df) <- c('x','y','hsm')
  
  shore <- shore_list[[i]]
  shore_pr <- projectRaster(from = shore, crs = CRS(proj_out))  # this takes a minute
  shore_crop <- crop(shore_pr, y = plot_ext)
  shore_df <- as.data.frame(shore_crop, xy = TRUE)
  shore_pol <- shore_df
  names(shore_pol) <- c('x','y','shore')
  shore_pol$shore <- ifelse(!is.na(shore_pol$shore), 1, 0)
  coordinates(shore_pol) <- ~x + y
  gridded(shore_pol) <- TRUE
  shore_pol <- raster(shore_pol)
  shore_cont <- rasterToContour(shore_pol)
  # simplify contour 
  shore_cont <- fortify(shore_cont)
  shore_cont <- shore_cont[shore_cont$id == 'C_1',]
  
  plot_list[[i]] <- ggplot() +
    geom_tile(data = df[!is.na(df$hsm),], aes(x = x, y = y, fill = hsm)) +
    scale_fill_gradientn(colours = terrain.colors(10, rev = TRUE)) +
    geom_path(data = lake_shp, aes(x = long, y = lat, group = group),
              size = 0.5, alpha = 0.4) +
    geom_polygon(data = glacier_list[[i]], aes(x = long, y = lat, group = group),
                 size = 0.5, fill = 'slategray1', color = 'slategray') +
    geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
              size = 0.5, alpha = 0.4) +
    geom_path(data = shore_cont, aes(x = long, y = lat, group = group), 
              size = 1, alpha = 0.8) +
    ggtitle(times[i]) +
    # labs(fill = 'Habitat \nsuitability') +
    theme(panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          legend.position = 'none',
          # legend.title = element_text(size = 8),
          # legend.text = element_text(size = 8),
          # legend.key.height = unit(1, 'cm'),
          plot.title = element_text(size = 10, hjust = 0.5)) +
    coord_fixed(xlim = c(100000, 2900000),
                ylim = c(-1000000, 1300000))
}
plot_grid(plotlist = plot_list)


#### LGM ZOOMED IN MAP
# convert suitability raster to dataframe
df <- as.data.frame(lgm, xy = TRUE)
names(df) <- c('x','y','hsm')

glacier_lgm <- glacier_list_pr[[21]]

# read shoreline shapefiles and convert to contour line
sh_lgm <- raster('shapefiles/shorelines/bio_2_lgm.tif')
shore_pr <- projectRaster(from = sh_lgm, crs = CRS(proj_out))  # this takes a minute
shore_crop <- crop(shore_pr, y = plot_ext)
shore_df <- as.data.frame(shore_crop, xy = TRUE)
shore_pol <- shore_df
names(shore_pol) <- c('x','y','shore')
shore_pol$shore <- ifelse(!is.na(shore_pol$shore), 1, 0)
coordinates(shore_pol) <- ~x + y
gridded(shore_pol) <- TRUE
shore_pol <- raster(shore_pol)
shore_cont <- rasterToContour(shore_pol)

# simplify contour 
shore_cont <- fortify(shore_cont)
shore_cont <- shore_cont[shore_cont$id == 'C_1',]  # should be rows 1:726


#### ZOOMED IN LGM MAP
plot_lgm_zoom <- ggplot() +
  geom_tile(data = df[!is.na(df$hsm),], aes(x = x, y = y, fill = hsm)) +
  scale_fill_gradientn(colours = terrain.colors(10, rev = TRUE)) +
  # geom_path(data = lake_shp, aes(x = long, y = lat, group = group),
  #           size = 0.5, alpha = 0.4) +
  geom_polygon(data = glacier_lgm, aes(x = long, y = lat, group = group),
               size = 0.5, fill = 'slategray1', color = 'slategray') +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
            size = 0.5, alpha = 0.4) +
  geom_path(data = shore_cont, aes(x = long, y = lat, group = group),
            size = 1, alpha = 0.8) +
  # labs(fill = 'Habitat \nsuitability') +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.position = 'none') +
        # plot.margin = unit(c(3,3,3,3), 'cm')) +
        # legend.title = element_text(size = 8),
        # legend.text = element_text(size = 8),
        # legend.key.height = unit(1, 'cm'),
        # plot.title = element_text(size = 12, hjust = 0.5)) +
  coord_fixed(xlim = c(1900000, 2900000),
              ylim = c(450000, 1300000))



#### PLOT LIG (130K YBP)
lig <- raster('shapefiles/HSM_rasters_dif_timepoints/lig_hsm.tif')
df <- as.data.frame(lig, xy = TRUE)
names(df) <- c('x','y','hsm')

sh_lig <- raster('shapefiles/shorelines/bio_2_lig.tif')
shore_pr <- projectRaster(from = sh_lig, crs = CRS(proj_out))  # this takes a minute
shore_crop <- crop(shore_pr, y = plot_ext)
shore_df <- as.data.frame(shore_crop, xy = TRUE)
shore_pol <- shore_df
names(shore_pol) <- c('x','y','shore')
shore_pol$shore <- ifelse(!is.na(shore_pol$shore), 1, 0)
coordinates(shore_pol) <- ~x + y
gridded(shore_pol) <- TRUE
shore_pol <- raster(shore_pol)
shore_cont <- rasterToContour(shore_pol)

# simplify contour 
shore_cont <- fortify(shore_cont)
shore_cont <- shore_cont[shore_cont$id == 'C_1',]  # should be rows 1:781

plot_lig <- ggplot() +
  geom_tile(data = df[!is.na(df$hsm),], aes(x = x, y = y, fill = hsm)) +
  scale_fill_gradientn(colours = terrain.colors(10, rev = TRUE)) +
  # geom_path(data = lake_shp, aes(x = long, y = lat, group = group),
  #           size = 0.5, alpha = 0.4) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
            size = 0.5, alpha = 0.4) +
  geom_path(data = shore_cont, aes(x = long, y = lat, group = group),
            size = 1, alpha = 0.8) +
  ggtitle('130,000 YBP') +
  # labs(fill = 'Habitat \nsuitability') +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.position = 'none',
        # legend.title = element_text(size = 8),
        # legend.text = element_text(size = 8),
        # legend.key.height = unit(1, 'cm'),
        plot.title = element_text(size = 10, hjust = 0.5)) +
  coord_fixed(xlim = c(100000, 2900000),
              ylim = c(-1000000, 1300000))


#### PLOT PRESENT DAY
present <- raster('shapefiles/HSM_rasters_dif_timepoints/present_hsm.tif')
df <- as.data.frame(present, xy = TRUE)
names(df) <- c('x','y','hsm')

plot_present <- ggplot() +
  geom_tile(data = df[!is.na(df$hsm),], aes(x = x, y = y, fill = hsm)) +
  scale_fill_gradientn(colours = terrain.colors(10, rev = TRUE)) +
  geom_path(data = lake_shp, aes(x = long, y = lat, group = group),
            size = 0.5, alpha = 0.4) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
            size = 0.5, alpha = 0.4) +
  # geom_path(data = shore_cont, aes(x = long, y = lat, group = group), 
  #           size = 0.5, alpha = 0.4) +
  ggtitle('Present') +
  # labs(fill = 'Habitat \nsuitability') +
  annotation_scale(location = "br", text_cex = 0.5) +
  annotation_north_arrow(location = "br", height = unit(0.75, 'cm'), width = unit(0.75, 'cm'),
                         pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in")) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.position = 'none',
        # legend.direction = 'horizontal',
        # legend.title = element_text(size = 12),
        # legend.text = element_text(size = 12),
        # legend.key.height = unit(1, 'cm'),
        # legend.key.width = unit(1.4, 'cm'),
        plot.title = element_text(size = 10, hjust = 0.5)) +
  coord_fixed(xlim = c(100000, 2900000),
              ylim = c(-1000000, 1300000))

plot_list[[8]] <- plot_present
plot_list[[9]] <- plot_lig
plot_list[[10]] <- plot_lgm_zoom
plot_list <- plot_list[c(9, 1, 10, 2:7, 8)]
plot_grid(plotlist = plot_list, nrow = 4, ncol= 3)
saveRDS(plot_list, 'hsm_shore_ice_figs_list.RDS')
ggsave('figs/hsm_shore_ice_maps.png', 
       height = 10, width = 9.4, units = 'in')



#### PLOT OF GENETIC SAMPLE LOCATIONS WITH INSET
require(dplyr)
require(ggplot2)
require(sp)
require(rgdal)
require(ggspatial)

dat <- read.csv('butternut_coord_df.csv', stringsAsFactors = FALSE)
names(dat) <- c('id', 'pop', 'color', 'x', 'y')

proj_from <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
proj_to <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# project coordinates
xy <- dat[,c('x','y')]
spdf <- SpatialPointsDataFrame(coords = xy, data = dat,
            proj4string = CRS(proj_from))
spdf <- spTransform(spdf, CRSobj = CRS(proj_to))
dat <- as.data.frame(spdf)

# create vector of colors for plotting
cols <- dat %>% select(color) %>% distinct()

# read in North America shapefiles
na_shp <- readOGR("shapefiles/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_to)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("shapefiles/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj_to)

# plot
ggplot() + 
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
            size = 0.5) +
  geom_polygon(data = lake_shp, aes(x = long, y = lat, group = group),
               size = 0.5, color = 'black', fill = 'white') +
  geom_point(data = dat, aes(x = x.1, y = y.1, fill = factor(pop)),
             size = 6, pch = 21, alpha = 0.6) +
  scale_fill_manual(name = 'Population ID', values = cols$color) +
  annotation_scale(location = "br", text_cex = 1.5) +
  annotation_north_arrow(location = "br",
                         pad_x = unit(0.8, "in"), pad_y = unit(0.4, "in")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        # legend.position = 'none') +
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key = element_rect(fill = "white")) +
        # legend.key.height = unit(1, 'cm'),
        # plot.title = element_text(size = 12, hjust = 0.5)) +
  coord_fixed(xlim = c(300000, 2500000),
              ylim = c(-700000, 1300000))
ggsave('figs/genetic_pops_map.png')


#### NORTH AMERICA INSET
ggplot() + 
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
            size = 0.5) +
  geom_polygon(data = lake_shp, aes(x = long, y = lat, group = group),
               size = 0.5, color = 'black', fill = 'white') +
  annotation_scale(location = "bl", text_cex = 2) +
  annotation_north_arrow(location = "bl",
                         pad_x = unit(0.8, "in"), pad_y = unit(0.4, "in")) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.position = 'none')
ggsave('figs/NA_inset_map.png')




#### RE-PLOTTING POPULATIONS SO COLOR SCHEME COINCIDES WITH REGION, NOT POP ID
require(dplyr)
require(ggplot2)
require(sp)
require(rgdal)
require(ggspatial)

dat <- read.csv('butternut_coord_df.csv', stringsAsFactors = FALSE)
names(dat) <- c('id', 'pop', 'color', 'x', 'y')
proj_from <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
proj_to <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# read in North America shapefiles
na_shp <- readOGR("shapefiles/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_to)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("shapefiles/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj_to)

# project coordinates
xy <- dat[,c('x','y')]
spdf <- SpatialPointsDataFrame(coords = xy, data = dat,
                               proj4string = CRS(proj_from))
spdf <- spTransform(spdf, CRSobj = CRS(proj_to))
dat <- as.data.frame(spdf)

# create color schemes by region
ont <- c('125147', '126147', '151')
queb <- c('171188', '170')
nb <- c('568', '9101113b', '9101113a', '1014', '7917', '31')
dat$region <- ifelse(dat$pop %in% ont, 'Ontario', 
                     ifelse(dat$pop %in% queb, 'Quebec', 
                            ifelse(dat$pop %in% nb, 'New Brunswick', 'United States')))
dat$reg_color <- ifelse(dat$region == 'Ontario', 'red4', 
                        ifelse(dat$region == 'Quebec', 'lightsalmon1', 
                               ifelse(dat$region == 'New Brunswick', 'firebrick1', 'dodgerblue2')))
# create vector of colors for plotting
cols <- dat %>% dplyr::select(reg_color) %>% distinct()

# plot
ggplot() + 
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
            size = 0.5) +
  geom_polygon(data = lake_shp, aes(x = long, y = lat, group = group),
               size = 0.5, color = 'black', fill = 'white') +
  geom_point(data = dat, aes(x = x.1, y = y.1, fill = factor(region)),
             size = 6, pch = 21, alpha = 0.6) +
  scale_fill_manual(name = 'Region', values = cols$reg_color) +
  
  annotation_scale(location = "br", text_cex = 1.5) +
  annotation_north_arrow(location = "br",
                         pad_x = unit(0.8, "in"), pad_y = unit(0.4, "in")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        # legend.position = 'none') +
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key = element_rect(fill = "white")) +
  # legend.key.height = unit(1, 'cm'),
  # plot.title = element_text(size = 12, hjust = 0.5)) +
  coord_fixed(xlim = c(300000, 2500000),
              ylim = c(-700000, 1300000))
ggsave('figs/genetic_pops_map_color_by_region.png')



#### RE-PLOTTING POPULATIONS USING MEAN POP LAT/LONG
require(dplyr)
require(ggplot2)
require(sp)
require(rgdal)
require(ggspatial)

dat <- read.csv('data/butternut_stat_df.csv', stringsAsFactors = FALSE)
dat <- dat %>% 
  select(X, Mean.Longitude, Mean.Latitude) %>% 
  rename(id = X, x = Mean.Longitude, y = Mean.Latitude)

proj_from <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
proj_to <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# read in North America shapefiles
na_shp <- readOGR("shapefiles/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_to)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("shapefiles/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj_to)

# project coordinates
xy <- dat[,c('x','y')]
spdf <- SpatialPointsDataFrame(coords = xy, data = data.frame(dat$id),
                               proj4string = CRS(proj_from))
spdf <- spTransform(spdf, CRSobj = CRS(proj_to))
dat <- as.data.frame(spdf)

# create color schemes by region
ont <- c('125147', '126147', '151', '171188', '170')
nb <- c('568', '9101113b', '9101113a', '1014', '7917', '31')
dat$region <- ifelse(dat$dat.id %in% ont, 'Ontario', 
                            ifelse(dat$dat.id %in% nb, 'New Brunswick', 'United States'))
dat$reg_color <- ifelse(dat$region == 'Ontario', 'red4', 
                               ifelse(dat$region == 'New Brunswick', 'firebrick1', 'dodgerblue2'))
# create vector of colors for plotting
cols <- dat %>% dplyr::select(reg_color) %>% distinct()

# plot
ggplot() + 
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
            size = 0.5) +
  geom_polygon(data = lake_shp, aes(x = long, y = lat, group = group),
               size = 0.5, color = 'black', fill = 'white') +
  geom_point(data = dat, aes(x = x, y = y, fill = factor(region)),
             size = 6, pch = 21, alpha = 0.9) +
  scale_fill_manual(name = 'Region', values = cols$reg_color) +
  
  annotation_scale(location = "br", text_cex = 1.5) +
  annotation_north_arrow(location = "br",
                         pad_x = unit(0.8, "in"), pad_y = unit(0.4, "in")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        # legend.position = 'none') +
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key = element_rect(fill = "white")) +
  # legend.key.height = unit(1, 'cm'),
  # plot.title = element_text(size = 12, hjust = 0.5)) +
  coord_fixed(xlim = c(300000, 2500000),
              ylim = c(-700000, 1300000))
ggsave('figs/genetic_pops_map_color_by_region_V3.png')

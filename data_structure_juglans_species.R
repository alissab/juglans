
require(neotoma)
require(dplyr)
require(tidyr)
require(analogue)
require(Bchron)
require(raster)
require(geosphere)
require(sp)
require(rgdal)
require(rdist)
require(rgeos)
require(ggplot2)
require(gridExtra)

# READ IN JUGLANS SITES AND DATA
dat.dl <- readRDS("juglans_data_download.RData")
data <- compile_downloads(dat.dl)

# remove NA rows
data <- data[!is.na(rowSums(data[,c('Juglans.nigra','Juglans.cinerea')], na.rm = TRUE)), ]

# remove rows with 0 sums
data <- data[rowSums(data[,c('Juglans.nigra','Juglans.cinerea')], na.rm = TRUE) > 0, ]

# RESCALE POLLEN COUNTS at each site before pulling out juglans data
# first, what is typical pollen count per site?
test <- data[,11:ncol(data)] %>% mutate(sum = rowSums(., na.rm = TRUE))
summary(test$sum)  # about 500

dat.perc <- data
dat.perc[ ,11:ncol(dat.perc)] <- t(apply(dat.perc[ ,11:ncol(dat.perc)],
                                         1,
                                         function(x) x/sum(x, na.rm=TRUE)*500))
# check your work
all(round(rowSums(dat.perc[, 11:ncol(dat.perc)], na.rm=TRUE)) == 500) 

# EXTRACT J. nigra and J. cinerea DATA SEPARATELY
dat_jug <- dat.perc %>% dplyr::select(site.name, age, date.type, 
                               lat, long, Juglans.nigra, Juglans.cinerea)

# calibrate age - leave this alone for now

# what proportion of counts contain each juglans species?


# for now, take sum of pollen counts per location if there are more than one 
# record/location
coords <- dat_jug[,c("site.name","age","lat","long")]
coords <- coords[!duplicated(coords[, c("site.name", "age")]), ]
dat_jug <- dat_jug %>% select(-date.type) %>% 
  group_by(site.name, age) %>% 
  summarise(juglans.nigra = sum(Juglans.nigra), 
            juglans.cinerea = sum(Juglans.cinerea))
dat_jug <- ungroup(dat_jug)
dat_jug <- left_join(dat_jug, coords, by = c("site.name", "age"))



# MAPS OF JUGLANS LOCATION BY TIME IN 1000-YEAR CHUNKS
# cut data into time chunks
dat.cut <- cut(dat_jug$age, breaks = c(min(dat_jug$age, na.rm = TRUE), seq(1000, 21000, by = 1000)),
               labels = FALSE)
dat_jug$cut <- dat.cut
dat_jug <- dat_jug[!is.na(dat_jug$cut),]

# prepare coordinates
dat_coords <- dat_jug[, c("long","lat")]
names(dat_coords) <- c("x", "y")
coordinates(dat_coords) <- dat_coords
sp::proj4string(dat_coords) <- proj_WGS84
dat_coords_t <- sp::spTransform(dat_coords, proj_out)
coords = coordinates(dat_coords_t)
coords <- data.frame(coords)
dat_jug <- cbind(dat_jug, coords)

dat_jug_sp <- SpatialPointsDataFrame(coords = dat_jug[,c('x','y')],
                                  data = dat_jug[,c('juglans.nigra','juglans.cinerea','cut')])

# crop out water
cropped <- raster::intersect(dat_jug_sp, shp)
cropped <- data.frame(cropped)

# separate dataframe into list of dataframes, one for each time chunk
time <- split(cropped, f = cropped$cut)

# PLOTTING LOOP OVER TIME INTERVALS (plot site locations of all Juglans spp.)
nchunks <- nrow(cropped[unique(cropped$cut),])
plot_list <- list()

for(i in 1:nchunks){
  plot_list[[i]] <- ggplot(data = time[[i]]) +
    geom_point(aes(x = x, y = y), size = 4, alpha = 0.5, color = "blue") +
    geom_path(data = shp, aes(x = long, y = lat, group = group)) +
    labs(title = paste0("-", i, "000")) +
    scale_y_continuous(limits = c(100000, 1900000)) +
    scale_x_continuous(limits = c(-800000, 2760000)) +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          legend.title = element_blank()) +
    coord_fixed()
}
do.call(grid.arrange, plot_list)




# PLOTTING LOOPS OVER TIME INTERVALS FOR TWO SPECIES SEPARATELY, CIRCLES
# FILLED BY RELATIVE ABUNDANCE

# separate dataframe into list of dataframes, one for each time chunk
time <- split(cropped, f = cropped$cut)
nigra_list <- lapply(time, "[", c(1, 3:5))
cinerea_list <- lapply(time, "[", c(2:5))
all_list <- list(nigra_list, cinerea_list)

# remove sites that contain NA or 0 pollen counts for each species/time chunk
for(j in 1:2){
  for(i in 1:21){
    all_list[[j]][[i]] <- subset(all_list[[j]][[i]], 
                                  !is.na(all_list[[j]][[i]][,1]) & 
                                    all_list[[j]][[i]][,1] != 0)
  }
}

plot_list <- rep(list(list()), 2)
for(j in 1:length(all_list)){
  for(i in 1:length(all_list[[1]])){
      plot_list[[j]][[i]] <- ggplot(data = all_list[[j]][[i]]) +
        geom_point(aes_string(x = 'x', y = 'y', fill = names(all_list[[j]][[i]][1])), 
                   size = 4, pch = 21, alpha = 0.3) +
        # geom_point(aes(x = x, y = y), size = 4, pch = 21, color = "black") +
        geom_path(data = shp, aes(x = long, y = lat, group = group)) +
        labs(subtitle = paste0("-", i, "000")) +
        scale_y_continuous(limits = c(-1600000, 3500000)) +
        scale_x_continuous(limits = c(-800000, 2760000)) +
        theme_classic() +
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              line = element_blank(),
              legend.position = 'none') +
        coord_fixed()
  }
  do.call(grid.arrange, plot_list[[j]])
}


# NEXT STEPS
# 1. try using bacon package to convert radiocarbon years to calendar years
# 2. apply cutoff (5%) to abundance
# 3. first, need to calculate proportion of grains that are juglans


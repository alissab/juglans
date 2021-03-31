
setwd('C:/Users/abrow/Documents/juglans')
require(neotoma)
require(Bchron)

# find datasets
all.datasets <- get_dataset(datasettype = "pollen")
all.downloads <- get_download(all.datasets, verbose = FALSE)
compiled.cores <- compile_downloads(all.downloads)
# saveRDS(compiled.cores, 'all_raw_pollen_data.RDS')
# compiled.cores <- readRDS('all_raw_pollen_data.RDS')

# standardize pollen counts so every site has same total # pollen grains across taxa
stand.cores <- compiled.cores
stand.cores[ ,11:ncol(stand.cores)] <- t(apply(stand.cores[ ,11:ncol(stand.cores)],
                                               1,
                                               function(x) if (sum(x, na.rm=TRUE)==0) {x} else {x/sum(x, na.rm=TRUE)*500}))

# check your work
summary(rowSums(stand.cores[,11:ncol(stand.cores)], na.rm=TRUE))
table(rowSums(stand.cores[,11:ncol(stand.cores)], na.rm=TRUE))
# 43 rows with 0 sums, 111581 with sum of 500

# get rid of rows with sum of 0 (ie, sites with no pollen counts)
stand.cores2 <- stand.cores[ rowSums(stand.cores[,11:ncol(stand.cores)], na.rm=TRUE) != 0, ]

# remove excess taxa; remove rows with butternut counts of 0 or NA
keep <- c('.id','site.name','depth','age','age.old','age.young','date.type','lat','long','dataset','Juglans.cinerea')
but <- compiled.cores[,colnames(compiled.cores) %in% keep]
but <- but[!is.na(but$Juglans.cinerea) & but$Juglans.cinerea != 0, ]

# translate any dates in radiocarbon years to calendar years
radio.years <- (but$date.type %in% "Radiocarbon years BP") &
  (but$age > 71 ) &
  (but$age < 46401)
sryears <- sum(radio.years, na.rm = TRUE)

# BChronCalibrate is in the BChron package:
calibrated <- BchronCalibrate(but$age[radio.years],
                              ageSds = rep(100, sryears),
                              calCurves = rep("intcal13", sryears))

#  we want the weighted means from "calibrated"
wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))
but$age[radio.years] <- sapply(calibrated, wmean.date)
hist(but$age)

# quick visual of all butternut plot locations
map <- map_data("world")
ggplot(data = data.frame(map), aes(long, lat)) + 
  geom_polygon(aes(group=group), color = "steelblue", alpha = 0.2) +
  geom_point(data = but,
             aes(x = long, y = lat), color = 2) +
  xlab("Longitude West") + 
  ylab("Latitude North")

# convert projection to USA Contiguous albers equal area
proj_out <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 
+ellps=GRS80 +datum=NAD83 +units=m +no_defs'

sp::coordinates(but) <- ~long+lat
sp::proj4string(but) <- sp::CRS('+init=epsg:4326')
albers <- sp::spTransform(but, proj_out)
xy = data.frame(coordinates(albers))
colnames(xy) = c('x', 'y')

# construct data frame with re-projected coordinates
but = data.frame(xy, but)
but = but[,which(colnames(but)!= 'optional')]

saveRDS(but, 'butternut_pollen_fixed_dates.RDS')
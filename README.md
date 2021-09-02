# Biogeographic history of *Juglans cinerea* (butternut)
This project focuses on the migration history of a tree of conservation concern in eastern North America over the last 21,000 years (i.e., since the Last Glacial Maximum).
The larger project (in prep) looks for evidence of postglacial migration in butternut using genetic, species occurrence/niche modeling, and fossil pollen data. This repo contains R code for the extraction, analysis, and visualization of fossil pollen data. For more information and code for the genetic analyses and species distribution modeling, see [ekschumacher's repo](https://github.com/ekschumacher/butternut).

### R
- **pull_butternut_data.R**
This code pulls fossil pollen records from the [Neotoma database](https://www.neotomadb.org/); performs data cleaning; and projects site coordinates in an Albers equal-area projection.
- **biotic_velocity.R** 
This code analyzes butternut pollen data (pollen record distance to glacier over time; migration speed/biotic velocity over time) and visualizes pollen records over time in relation to glacial coverage and habitat suitability estimates. This code uses [Dalton et al.'s (2020)](https://www.sciencedirect.com/science/article/abs/pii/S0277379119307619) ice sheet shapefiles. 

### Shapefiles
- North American shapefiles
- Great Lakes (North America) shapefiles

### Data
- **butternut_pollen_fixed_dates.RDS**
Contains output of pull_butternut_data.R code (ie, butternut fossil pollen data)
- **butternut_stat_df_csv**
Contains lat/long for each population used in genetic analyses. Used to create map of population locations (Fig 1)

### Figs
This folder contains figures used in main text and supplement of manuscript.

library(tidyverse); library(lubridate)

# Function reads and formats GP5WSHELL csv export files

read_gp5 <- function(x)
{
  require(tidyverse)
  read.csv(x, skip = 2, header = TRUE) %>% # Rows to skip depend on GP5WShell version!!
    dplyr::rename(Temperature = X.1.oC) %>%
    mutate(Time = as.POSIXct(Time, 
                             format =" %d.%m.%Y %H:%M:%S", 
                             tz = "UTC"),
           Logger = gsub("\\_.*", "" , x)) %>%
    dplyr::select(Logger, Time, Temperature)
}

# Read logger header data

read.table("data/Loggers.txt", sep = "\t", header = TRUE) %>%
  mutate(Installed = as.POSIXct(Installed, 
                                format ="%d/%m/%Y", 
                                tz = "UTC")) %>%
  select(-c(Name, ED50Zone:ED50Y, Datum)) -> header

header %>%
  group_by(Site) %>%
  summarise(Elevation = round(mean(Elevation), 0)) -> average.elevation

header %>% select(-Elevation) %>% 
  merge(average.elevation) -> header # Only one elevation value per site

# Read and bind all loggers, merge with logger header data

wd <- getwd()
setwd(paste(wd, "/data", sep = ""))
do.call(rbind, lapply(list.files(pattern = "*.csv"), read_gp5)) %>% 
  merge(header, by = "Logger") %>%
  na.omit -> logs
setwd(wd)

# Keep only records with two working loggers per site

logs %>% 
  group_by(Site, Groundwater) %>%
  summarise(Last.record = max(Time)) %>%
  group_by(Site) %>%
  summarise(End = min(Last.record)) -> ends

logs %>% 
  group_by(Site, Groundwater) %>%
  summarise(First.record = min(Time)) %>%
  group_by(Site) %>%
  summarise(Start = max(First.record)) -> starts

merge(starts, ends) -> period

merge(logs, period, by = "Site") %>%
  group_by(Logger) %>%
  filter(Time > (Installed + 60*60*24*7)) %>% # Remove records up to one week after installation
  filter(Time >= Start & Time <= End) %>% # Remove logs with only one working logger per site
  dplyr::select(Logger, Time, Temperature) ->
  temperatures

# Get CHELSA bioclimatic variables (they are in a external folder, they need to be downloaded https://chelsa-climate.org/bioclim/)

bio2 <- raster::raster("../#wwfmap/CHELSA/CHELSA_bio10_02.tif") # Diurnal range from CHELSA
bio5 <- raster::raster("../#wwfmap/CHELSA/CHELSA_bio10_05.tif") # Max of the warmest from CHELSA
bio6 <- raster::raster("../#wwfmap/CHELSA/CHELSA_bio10_06.tif") # Min of the coldest from CHELSA

pts <- header
sp::coordinates(pts) <- ~Longitude+Latitude # Make points a spatial object
sp::proj4string(pts) <- sp::CRS("+proj=longlat + ellps=WGS84") # Asign projection
Chelsa <- data.frame(sp::coordinates(pts)) # Create data.frame with the Enscobase points
Chelsa$`Diurnal range` <- raster::extract(bio2, pts)/10 # Extract diurnal range from CHELSA layer, correct units
Chelsa$`Summer max` <- raster::extract(bio5, pts)/10 # Extract max warmest from CHELSA layer, correct units
Chelsa$`Winter min` <- raster::extract(bio6, pts)/10 # Extract min coldest from CHELSA layer, correct units
Chelsa$`Annual range` <- Chelsa$`Summer max` - Chelsa$`Winter min` # Calculate bio7, annual range
header <- merge(header, Chelsa, by = c("Longitude", "Latitude"))

# Save clean data

save(temperatures, header, file = "results/logs.RData")

library(tidyverse); library(sp); library(lubridate)

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

# # Functions convert coordinates to decimal degrees
# 
# decdeg30T <- function(x) { x -> o1 
#   o1 %>% dplyr::select(ED50X, ED50Y) -> o2
#   SpatialPoints(o2, proj4string = CRS("+proj=utm +zone=30T +datum=WGS84")) -> o3
#   spTransform(o3, CRS("+proj=longlat +datum=WGS84")) %>% data.frame -> o4
#   colnames(o4) <- c("Longitude", "Latitude")
#   o1 %>% dplyr::select(-c(ED50X, ED50Y)) -> o5
#   cbind(o5, o4)
# }
# 
# decdeg29T <- function(x) { x -> o1 
#   o1 %>% dplyr::select(ED50X, ED50Y) -> o2
#   SpatialPoints(o2, proj4string = CRS("+proj=utm +zone=29T +datum=WGS84")) -> o3
#   spTransform(o3, CRS("+proj=longlat +datum=WGS84")) %>% data.frame -> o4
#   colnames(o4) <- c("Longitude", "Latitude")
#   o1 %>% dplyr::select(-c(ED50X, ED50Y)) -> o5
#   cbind(o5, o4)
# }

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

# Map of the loggers

header %>%
  ggplot(aes(x = Longitude, y = Latitude)) + 
  labs(x = "Longitude (º)", y = "Latitude (º)") +
  geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), 
               color = "black", fill = "grey") + 
  geom_point(aes(color = "indianred"), alpha = 0.5, size = 3, show.legend = FALSE) +
  geom_text(aes(x = -3, y = 42, label = "Spain"), size = 4) +
  geom_text(aes(x = 1, y = 44.8, label = "France"), size = 4) +
  geom_text(aes(x = -7.7, y = 40.8, label = "Portugal"), size = 4) +
  geom_text(aes(x = -5.5, y = 45.5, label = "Atlantic"), size = 3) +
  coord_cartesian(xlim = c(-10, 2), ylim = c(40, 46)) +
  theme(legend.position = "top", legend.justification = "center",
        panel.background = element_rect(fill = "darkgrey", colour = "darkgrey"),
        legend.text = element_text(size = 8), legend.title = element_blank())


# Read and bind all loggers, merge with logger header data
# Could not use filepath, 
# " Error in file(file, "rt") : cannot open the connection"

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
  filter(Time >= Start & Time <= End) -> # Remove logs with only one working logger per site
  logs.clean

# Visual check of the temperatures

logs.clean %>%
  ggplot(aes(x = Time, y = Temperature, color = Groundwater)) +
  geom_line() +
  facet_wrap(~ Site)

# Waterlogged logger at Los Cándanos looks like it sunk into the bog

logs.clean %>%
  filter(Time < "2014-10-25") %>%
  filter(Site == "Los Cándanos" & Groundwater == "Waterlogged") %>%
  group_by(day = floor_date(Time, "day")) %>%
  summarise(dif = max(Temperature) - min(Temperature)) %>%
  ggplot(aes(x = day, y = dif)) + 
  geom_line() # Problem seems to be around 2014-10-20

# Because of the few records I'll remove this site

logs.clean %>% 
  filter(! Site %in% c("Los Cándanos")) -> logs.clean2

# Visual check again

logs.clean2 %>%
  ggplot(aes(x = Time, y = Temperature, color = Groundwater)) +
  geom_line() +
  facet_wrap(~ Site)

# Final clean file

logs.clean2 %>%
  dplyr::select(Logger, Time, Temperature) ->
  temperatures

header %>%
  filter(Logger %in% temperatures$Logger) ->
  header

# Get CHELSA bioclimatic variables

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

write.csv(temperatures, "results/temperatures.csv", row.names = FALSE)
write.csv(header, "results/header.csv", row.names = FALSE)
save(temperatures, header, file = "results/logs.RData")

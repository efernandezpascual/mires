library(tidyverse); library(ggthemes); library(lubridate)

# Calculations

## Bioclimatic variables US Geological Survey, WorlClim

### Bio 2 Annual Mean Diurnal Range

temperatures %>%
  group_by(Logger, Day = floor_date(Time, "day")) %>%
  summarise(Tmax = max(Temperature), 
            Tmin = min(Temperature),
            Tdif = Tmax - Tmin) %>%
  group_by(Logger) %>%
  summarise(Mean.Diurnal.Range = mean(Tdif)) -> bio2

### Bio 5 Max Temperature of Warmest Month

temperatures %>%
  group_by(Logger, Day = floor_date(Time, "day")) %>%
  summarise(Tmax = max(Temperature)) %>%
  mutate(Month = month(Day)) %>%
  group_by(Logger, Month) %>%
  summarise(Max.Temperature.of.Warmest.Month = mean(Tmax)) %>%
  group_by(Logger) %>%
  filter(Max.Temperature.of.Warmest.Month == max(Max.Temperature.of.Warmest.Month)) %>%
  select(-Month) -> bio5

### Bio 6 Min Temperature of the Coldest Month

temperatures %>%
  group_by(Logger, Day = floor_date(Time, "day")) %>%
  summarise(Tmin = min(Temperature)) %>%
  mutate(Month = month(Day)) %>%
  group_by(Logger, Month) %>%
  summarise(Min.Temperature.of.Coldest.Month = mean(Tmin)) %>%
  group_by(Logger) %>%
  filter(Min.Temperature.of.Coldest.Month == min(Min.Temperature.of.Coldest.Month)) %>%
  select(-Month) -> bio6

### Merge and calculate Bio 7 Annual Temperature Range

bio2 %>%
  merge(bio5) %>%
  merge(bio6) %>%
  mutate(Temperature.Annual.Range = 
           Max.Temperature.of.Warmest.Month -
           Min.Temperature.of.Coldest.Month) -> bioclimatics

## Paired T-tests

t.testEFP <- function(df) 
{
  t <- t.test(df$Waterlogged, df$Dry, alternative = "less", paired = TRUE)
  data.frame(t$statistic,
             t$p.value,
             t$estimate)
}

t.testEFPgreater <- function(df) 
{
  t <- t.test(df$Waterlogged, df$Dry, alternative = "greater", paired = TRUE)
  data.frame(t$statistic,
             t$p.value,
             t$estimate) # In the case of the minimum temperature is the opposite hypothesis
}

rbind(
  merge(header, bioclimatics, by = "Logger") %>% 
    gather(Trait, Value, Mean.Diurnal.Range:Temperature.Annual.Range) %>%
    select(Site, Groundwater, Value, Trait, Value) %>%
    spread(Groundwater, Value) %>%
    filter(! Trait %in% "Min.Temperature.of.Coldest.Month") %>%
    group_by(Trait) %>%
    do(t.testEFP(.)),
  merge(header, bioclimatics, by = "Logger") %>% 
    gather(Trait, Value, Mean.Diurnal.Range:Temperature.Annual.Range) %>%
    select(Site, Groundwater, Value, Trait, Value) %>%
    spread(Groundwater, Value) %>%
    filter(Trait %in% "Min.Temperature.of.Coldest.Month") %>%
    group_by(Trait) %>%
    do(t.testEFPgreater(.))) -> ttests

#----------------------------------------------------------------------------

# Tables

## Table 1 - Study sites

temperatures %>%
  group_by(Logger) %>%
  summarise(Start = min(Time),
            End = max(Time),
            Length = End - Start,
            Length = round(as.numeric(Length), 0)) -> 
  recording.period

merge(header, recording.period) %>%
  group_by(Site, Mire, Elevation, Length) %>%
  summarise(Latitude = mean(Latitude), 
            Longitude = mean (Longitude)) %>%
  select(Site, Mire, Elevation, Latitude, Longitude, Length) %>%
  mutate(Latitude = round(Latitude, 4),
         Longitude = round(Longitude, 4)) %>%
  rename(Habitat = Mire,
         `Elevation (m)` = Elevation,
         `Latitude (ºC)` = Latitude,
         `Longitude (ºC)` = Longitude,
         `Records (days)` = Length) -> table1

## Table 2 - Bioclimatic variables

merge(header, bioclimatics) %>%
  select(Site, Groundwater, Mean.Diurnal.Range:Temperature.Annual.Range) %>%
  gather(Trait, Value, Mean.Diurnal.Range:Temperature.Annual.Range) %>%
  mutate(Trait = gsub("\\.", " ", Trait), 
         Value = round(Value, 2)) %>%
  mutate(Column = paste(Trait, Groundwater, sep = "\n")) %>%
  select(-c(Trait, Groundwater)) %>%
  spread(Column, Value) -> table2 # add header with https://rdrr.io/cran/kableExtra/man/add_header_above.html

# Figures

## Figure 1 - Datalogger records

merge(temperatures, header) %>%
  arrange(Elevation, Site) %>%
  mutate(Facet = paste(Elevation, Site, sep = " m a.s.l. - ")) -> recordsdf

recordsdf %>%
  select(Elevation, Facet) %>%
  unique %>%
  arrange(-Elevation) %>%
  pull(Facet) -> Order

recordsdf$Facet <- factor(recordsdf$Facet, levels = Order)

recordsdf %>%
  ggplot(aes(x = Time, y = Temperature, color = Groundwater)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(alpha = 0.85) +
  facet_wrap(~ Facet, ncol = 1) +
  xlab("Time (hourly records)") + ylab("Temperature (ºC)") +
  #scale_color_manual(values = c("darkorchid", "gold")) +
  theme_tufte() +
  theme(text = element_text(size = 10),
        legend.position = "top", 
        legend.justification = "center",
        legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha=1))) -> plot1

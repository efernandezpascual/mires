library(tidyverse); library(ggthemes); library(lubridate); library(here)

load(here::here("results", "logs.RData"))

# Calculations

## Bioclimatic variables US Geological Survey, WorlClim

### Bio 2 Annual Mean Diurnal Range

temperatures %>%
  group_by(Logger, Day = floor_date(Time, "day")) %>%
  summarise(Tmax = max(Temperature), 
            Tmin = min(Temperature),
            Tdif = Tmax - Tmin) %>%
  group_by(Logger) %>%
  summarise(`Diurnal range` = mean(Tdif)) -> bio2

### Bio 5 Max Temperature of Warmest Month

temperatures %>%
  group_by(Logger, Day = floor_date(Time, "day")) %>%
  summarise(Tmax = max(Temperature)) %>%
  mutate(Month = month(Day)) %>%
  group_by(Logger, Month) %>%
  summarise(`Summer max` = mean(Tmax)) %>%
  group_by(Logger) %>%
  filter(`Summer max` == max(`Summer max`)) %>%
  select(-Month) -> bio5

### Bio 6 Min Temperature of the Coldest Month

temperatures %>%
  group_by(Logger, Day = floor_date(Time, "day")) %>%
  summarise(Tmin = min(Temperature)) %>%
  mutate(Month = month(Day)) %>%
  group_by(Logger, Month) %>%
  summarise(`Winter min` = mean(Tmin)) %>%
  group_by(Logger) %>%
  filter(`Winter min` == min(`Winter min`)) %>%
  select(-Month) -> bio6

### Merge and calculate Bio 7 Annual Temperature Range

bio2 %>%
  merge(bio5) %>%
  merge(bio6) %>%
  mutate(`Annual range` = 
           `Summer max` -
           `Winter min`) -> bioclimatics

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
  merge(header[, 1:8], bioclimatics, by = "Logger") %>% 
    gather(Trait, Value, `Diurnal range`:`Annual range`) %>%
    select(Site, Groundwater, Value, Trait, Value) %>%
    spread(Groundwater, Value) %>%
    filter(! Trait %in% "Winter min") %>%
    group_by(Trait) %>%
    do(t.testEFP(.)),
  merge(header[, 1:8], bioclimatics, by = "Logger") %>% 
    gather(Trait, Value, `Diurnal range`:`Annual range`) %>%
    select(Site, Groundwater, Value, Trait, Value) %>%
    spread(Groundwater, Value) %>%
    filter(Trait %in% "Winter min") %>%
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
         `Latitude` = Latitude,
         `Longitude` = Longitude,
         `Records (days)` = Length) %>%
  arrange(`Elevation (m)`) -> table1

## Table 2 - Bioclimatic variables

merge(header[, 1:8], bioclimatics) %>%
  select(Site, Groundwater, `Diurnal range`:`Annual range`) %>%
  gather(Trait, Value, `Diurnal range`:`Annual range`) %>%
  spread(Groundwater, Value) %>%
  mutate(Buffer = Waterlogged - Dry,
         Trait = gsub("\\.", " ", Trait), 
         Buffer = round(Buffer, 2)) %>%
  select(-c(Dry, Waterlogged)) %>%
  spread(Trait, Buffer) %>%
  merge(header[, c(3, 8)], by = "Site") %>%
  arrange(Elevation) %>%
  select(Site, `Annual range`:`Winter min`) %>%
  unique -> table2 

## Table 3 - Chelsa

bioclimatics %>%
  gather(Variable, Value, `Diurnal range`:`Annual range`) %>%
  mutate(Source = "Loggers") -> bioLoggers

header %>%
  select(Logger, `Diurnal range`:`Annual range`) %>%
  gather(Variable, Value, `Diurnal range`:`Annual range`) %>%
  mutate(Source = "CHELSA") -> bioCHELSA

header %>%
  select(Site, Logger, Groundwater, Mire) %>%
  merge(rbind(bioLoggers, bioCHELSA), by = "Logger") %>%
  spread(Source, Value) -> bio

flm <- function(df) {
  m <- lm(Loggers ~ CHELSA, data = df)
  s <- summary(m)
  data.frame(t = s$coefficients[2, 3],
             p = s$coefficients[2, 4],
             R2 = s$adj.r.squared)}

bio %>%
  group_by(Variable, Groundwater) %>%
  do(flm(.)) %>%
  mutate(t =  round(t, 3),
         p = round(p, 3),
         R2 = round(R2, 2)) -> table3

# Figures

## Figure 1 - Datalogger records

merge(temperatures, header) %>%
  arrange(Elevation, Site) %>%
  mutate(Facet = paste(Elevation, Site, sep = " m a.s.l. - ")) -> recordsdf

recordsdf %>%
  select(Elevation, Facet) %>%
  unique %>%
  arrange(Elevation) %>%
  pull(Facet) -> Order

recordsdf$Facet <- factor(recordsdf$Facet, levels = Order)

recordsdf %>%
  ggplot(aes(x = Time, y = Temperature, color = Groundwater)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(alpha = 0.85) +
  facet_wrap(~ Facet, ncol = 4) +
  xlab("Time (hourly records)") + ylab("Temperature (ºC)") +
  ggthemes::theme_tufte() +
  theme(panel.background = element_rect(color = "grey96", fill = "grey96"),
        text = element_text(size = 10),
        legend.position = "top", 
        legend.justification = "center",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_color_manual(values = c("red3", "turquoise3")) -> plot1

## Fig 2 - Effect sizes

merge(header[, 1:8], bioclimatics) %>%
  select(Site, Groundwater, `Diurnal range`:`Annual range`) %>%
  gather(Trait, Value, `Diurnal range`:`Annual range`) %>%
  mutate(Trait = gsub("\\.", " ", Trait), 
         Value = round(Value, 2)) %>%
  group_by(Trait, Groundwater) %>%
  summarise(y = mean(Value), SD = sd(Value), n = length(Value),
            SE = SD/sqrt(n)) %>%
  ggplot(aes(x = Groundwater, y = y, ymin = y - SE, ymax = y + SE,
             fill = Groundwater)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~ Trait, scales = "free", strip.position = "left", ncol = 4) +
  geom_errorbar(width = .2, position = position_dodge(.9), 
                color = "grey25", show.legend = FALSE) +
  ggthemes::theme_tufte() + 
  theme(panel.background = element_rect(color = "grey96", fill = "grey96"), 
        legend.position = "top", 
        legend.justification = "center",
        legend.title = element_blank(),
        strip.placement = "outside",
        axis.title = element_blank()) +
  scale_fill_manual(values = c("red3", "turquoise3")) -> plot2

## Fig 3 Chelsa

bio %>%
  ggplot(aes(CHELSA, Loggers, color = Groundwater)) + 
  geom_point(size = 4) + 
  geom_smooth(method = "lm", se = F, linetype = "dashed") +
  facet_wrap( ~Variable, ncol = 2, scales = "free") +
  ggthemes::theme_tufte() + 
  xlab("CHELSA air temperatures (°C)") + ylab("Logger soil temperatures (°C)") +
  # coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme(panel.background = element_rect(color = "grey96", fill = "grey96"), 
        legend.position = "top", 
        legend.justification = "center",
        legend.title = element_blank(),
        strip.placement = "outside") +
  scale_color_manual(values = c("red3", "turquoise3")) -> plot3

# Save

save(ttests, table1, table2, table3, plot1, plot2, plot3, file = here::here("results", "MSoutput.RData"))

ggsave(plot1, file = here::here("doc", "manuscript", "Fig1.png"), 
       path = NULL, scale = 1, width = 170, height = 120, units = "mm", dpi = 600)

ggsave(plot2, file = here::here("doc", "manuscript", "Fig2.png"), 
       path = NULL, scale = 1, width = 170, height = 70, units = "mm", dpi = 600)

ggsave(plot3, file = here::here("doc", "manuscript", "Fig3.png"), 
       path = NULL, scale = 1, width = 170, height = 150, units = "mm", dpi = 600)
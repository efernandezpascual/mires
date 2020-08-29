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

#----------------------------------------------------------------------------

# Figures

## Figure 1 - Map

library(raster)

### Get map of Spain

spaingeo <- getData("GADM", country = "spain", level = 1)
spaindf <- fortify(spaingeo)

spainalt <- getData("alt", country = "spain", level = 1)
altdf  <- data.frame(rasterToPoints(spainalt))
colnames(altdf) = c("lon", "lat", "alt")
altdf %>% filter(alt < 2650) -> altdf

### Draw map of Spain

header %>%
  group_by(Site, Mire) %>%
  summarise(long = mean(Longitude), lat = mean(Latitude)) %>%
  ggplot(aes(long, lat)) +
  geom_tile(data = altdf, aes(x = lon, y = lat, fill = alt)) +
  scale_fill_gradient(name = "Elevation (m asl)", low = "grey96", high = "black",
                      breaks = c(0, 2600)) +
  geom_blank(data = spaindf, aes(x = long, y = lat)) +
  geom_map(data = spaindf, map = spaindf,
           aes(group = group, map_id = id),
           fill = NA, color = "black", size = .3) +
  geom_point(size = 2) +
  ggrepel::geom_label_repel(aes(label = Site), box.padding = .6, size = 2.5) +
  labs(x = "Longitude (º)", y = "Latitude (º)") +
  geom_text(aes(x = -6.1, y = 43.75, label = "Cantabrian Sea"), size = 3.5) +
  geom_text(aes(x = -6.1, y = 43.4, label = "Asturies"), size = 5) +
  geom_text(aes(x = -5.1, y = 42.9, label = "León"), size = 5) +
  scale_x_continuous(limits = c(-6.5, -4.7), expand = c(0, 0)) +
  scale_y_continuous(limits = c(42.8, 43.9), expand = c(0, 0)) +
  ggthemes::theme_map() +
  theme(panel.background = element_rect(color = "black", fill = "grey96"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, 0, 0),
        plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
        axis.title = element_blank()) -> g1

### Draw map of Europe

ggplotGrob(
  ggplot() +
    geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), 
                 color = "black", fill = "grey", size = 0.25) +
    geom_path(data = data.frame(long = c(-6.5, -6.5, -4.7, -4.7, -6.5),
                                lat = c(42.8, 43.9, 43.9, 42.8, 42.8)), 
              aes(x = long, y = lat), size = 0.25) +
    labs(x = "Longitude (º)", y = "Latitude (º)") +
    coord_cartesian(xlim = c(-8.8, 25), ylim = c(36.5, 58.5)) +
    ggthemes::theme_map() +
    theme(panel.background = element_rect(color = "grey96", fill = "grey96"),
          legend.position = "top", legend.title = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))) -> g2
 
### Put together

g1 +
  annotation_custom(grob = g2, xmin = -5.54, xmax = -4.7,
                    ymin = 43.445, ymax = 43.9) -> plot1a

### Add picture

img1 <- png::readPNG("data/riotuertu.png")
plot1b <- grid::rasterGrob(img1, width = unit(85, "mm"), height = unit(85, "mm"), just = "centre") 

cowplot::plot_grid(plot1a, plot1b, ncol = 2, labels = c("a", "b")) -> plot1

## Figure 2 - Datalogger records

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
  scale_color_manual(values = c("red3", "turquoise3")) -> plot2

## Fig 3 - Effect sizes

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
  scale_fill_manual(values = c("red3", "turquoise3")) -> plot3

## Fig 4 Chelsa

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
  scale_color_manual(values = c("red3", "turquoise3")) -> plot4

# Save

save(ttests, table1, table2, table3, plot1, plot2, plot3, plot4, file = here::here("results", "MSoutput.RData"))

ggsave(plot1, file = here::here("doc", "manuscript", "Fig1.png"), 
       path = NULL, scale = 1, width = 170, height = 85, units = "mm", dpi = 600)

ggsave(plot2, file = here::here("doc", "manuscript", "Fig2.png"), 
       path = NULL, scale = 1, width = 170, height = 120, units = "mm", dpi = 600)

ggsave(plot3, file = here::here("doc", "manuscript", "Fig3.png"), 
       path = NULL, scale = 1, width = 170, height = 70, units = "mm", dpi = 600)

ggsave(plot4, file = here::here("doc", "manuscript", "Fig4.png"), 
       path = NULL, scale = 1, width = 170, height = 150, units = "mm", dpi = 600)
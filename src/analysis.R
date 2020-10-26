library(tidyverse); library(lubridate)

load("results/logs.RData")

# Calculations

## Bioclimatic variables as defined by US Geological Survey, WorlClim

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

## Check for normality http://www.sthda.com/english/wiki/paired-samples-t-test-in-r#preleminary-test-to-check-paired-t-test-assumptions

merge(header[, 1:8], bioclimatics, by = "Logger") %>% 
  gather(Trait, Value, `Diurnal range`:`Annual range`) %>%
  dplyr::select(Site, Groundwater, Value, Trait, Value) %>%
  spread(Groundwater, Value) %>%
  mutate(Buffer = Waterlogged - Dry) -> o1

o1 %>%
  filter(Trait == "Annual range") %>%
  pull(Buffer) %>%
  shapiro.test()

o1 %>%
  filter(Trait == "Diurnal range") %>%
  pull(Buffer) %>%
  shapiro.test()

o1 %>%
  filter(Trait == "Winter min") %>%
  pull(Buffer) %>%
  shapiro.test()

o1 %>%
  filter(Trait == "Summer max") %>%
  pull(Buffer) %>%
  shapiro.test()

## Paired T-tests

t.testEFP <- function(df) 
{
  t <- t.test(df$Waterlogged, df$Dry, alternative = "less", paired = TRUE)
  data.frame(t$statistic,
             t$p.value,
             t$estimate,
             t$conf.int[1],
             t$conf.int[2]
  )
}

t.testEFPgreater <- function(df) 
{
  t <- t.test(df$Waterlogged, df$Dry, alternative = "greater", paired = TRUE)
  data.frame(t$statistic,
             t$p.value,
             t$estimate,
             t$conf.int[1],
             t$conf.int[2]
  ) # In the case of the minimum temperature is the opposite hypothesis
}

rbind(
  merge(header[, 1:8], bioclimatics, by = "Logger") %>% 
    gather(Trait, Value, `Diurnal range`:`Annual range`) %>%
    dplyr::select(Site, Groundwater, Value, Trait, Value) %>%
    spread(Groundwater, Value) %>%
    filter(! Trait %in% "Winter min") %>%
    group_by(Trait) %>%
    do(t.testEFP(.)),
  merge(header[, 1:8], bioclimatics, by = "Logger") %>% 
    gather(Trait, Value, `Diurnal range`:`Annual range`) %>%
    dplyr::select(Site, Groundwater, Value, Trait, Value) %>%
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
  group_by(Site, Mire, pH, Conductivity, Elevation, Length) %>%
  summarise(Latitude = mean(Latitude), 
            Longitude = mean (Longitude)) %>%
  select(Site, Mire, pH, Conductivity, Elevation, Latitude, Longitude, Length) %>%
  mutate(Latitude = round(Latitude, 4),
         Longitude = round(Longitude, 4),
         pH = round(pH, 1),
         Conductivity = round(Conductivity, 1)) %>%
  rename(Habitat = Mire,
         `Conductivity (μS/cm)` = Conductivity,
         `Elevation (m)` = Elevation,
         `Latitude` = Latitude,
         `Longitude` = Longitude,
         `Records (days)` = Length) %>%
  arrange(`Elevation (m)`) -> table1

## Table 2 - Bioclimatic variables

merge(header[, 1:10], bioclimatics) %>%
  dplyr::select(Site, Groundwater, `Diurnal range`:`Annual range`) %>%
  gather(Trait, Value, `Diurnal range`:`Annual range`) %>%
  spread(Groundwater, Value) %>%
  mutate(Buffer = Waterlogged - Dry,
         Trait = gsub("\\.", " ", Trait), 
         Buffer = round(Buffer, 2)) %>%
  select(-c(Dry, Waterlogged)) %>%
  spread(Trait, Buffer) %>%
  merge(header[, c(3, 10)], by = "Site") %>%
  arrange(Elevation) %>%
  select(Site, `Annual range`:`Winter min`) %>%
  unique %>%
  remove_rownames() -> table2 

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

RMSE <- function(df) {
  require(Metrics)
  data.frame(RMSE = rmse(df$Loggers, df$CHELSA))}

CORR <- function(df) {
  c <- cor.test(df$Loggers, df$CHELSA)
  data.frame(t = c$statistic,
             df = c$parameter,
             p = c$p.value,
             r = c$estimate)}

bio %>%
  group_by(Variable, Groundwater) %>%
  do(CORR(.)) -> corrtable # Pearson's correlation

bio %>%
  group_by(Variable, Groundwater) %>%
  do(RMSE(.)) -> rmsetable # Root-mean-square error (RMSE)

merge(corrtable, rmsetable) %>%
  mutate(RMSE = round(RMSE, 2),
         r = round(r, 2),
         t = round(t, 3),
         p = round(p, 3)) %>%
  rename(`RMSE (ºC)` = RMSE,
         `Pearson's r` = r) -> table3

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
        legend.title = element_text(size = 16, face = "bold"),
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

img1 <- png::readPNG("data/Riotuertu.png")
plot1b <- grid::rasterGrob(img1, width = unit(85, "mm"), height = unit(85, "mm"), just = "centre") 

cowplot::plot_grid(plot1a, plot1b, ncol = 2, labels = c("(a)", "(b)")) -> plot1

## Figure 2 - Datalogger records

merge(temperatures, header) %>%
  arrange(Elevation, Site) %>%
  mutate(Facet = paste(Elevation, Mire, sep = " m a.s.l. - ")) -> recordsdf

recordsdf %>%
  dplyr::select(Elevation, Facet) %>%
  unique %>%
  arrange(-Elevation) %>%
  pull(Facet) -> Order

recordsdf$Facet <- factor(recordsdf$Facet, levels = Order)

mean(recordsdf$Time) -> MeanTime

recordsdf %>%
  group_by(Facet, Site) %>%
  summarise(Time = MeanTime, Temperature = 40) -> labelsdf

recordsdf %>%
  ggplot(aes(x = Time, y = Temperature)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(aes(color = Groundwater), alpha = 0.65) +
  facet_wrap(~ Facet, ncol = 2) +
  xlab("Time (hourly records)") + ylab("Temperature (ºC)") +
  geom_text(aes(label = Facet, y = Temperature - 5), data = labelsdf, fontface = "bold") + 
  geom_text(aes(label = Site), data = labelsdf, fontface = "bold") + 
  scale_x_datetime(breaks = as.POSIXct(c("2015-01-01 12:00:00 UTC",
                                         "2016-01-01 12:00:00 UTC",
                                         "2017-01-01 12:00:00 UTC",
                                         "2018-01-01 12:00:00 UTC", 
                                         "2019-01-01 12:00:00 UTC")),
                   date_labels = "%Y") +
  ggthemes::theme_tufte() +
  theme(panel.background = element_rect(color = "grey96", fill = "grey96"),
        legend.position = "top", 
        legend.justification = "center",
        legend.title = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, -7.5, 0),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  scale_color_manual(values = c("red3", "turquoise3")) -> plot2

## Fig 3 - Effect sizes

merge(header[, 1:8], bioclimatics) %>%
  dplyr::select(Site, Groundwater, `Diurnal range`:`Annual range`) %>%
  gather(Trait, Value, `Diurnal range`:`Annual range`) %>%
  mutate(Trait = gsub("\\.", " ", Trait), 
         Value = round(Value, 2)) %>%
  group_by(Trait, Groundwater) %>%
  summarise(y = mean(Value), SD = sd(Value), n = length(Value),
            SE = SD/sqrt(n)) %>%
  ggplot(aes(x = Groundwater, y = y, ymin = y - SE, ymax = y + SE,
             fill = Groundwater)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.65) +
  facet_wrap(~ Trait, scales = "free", strip.position = "left", ncol = 4) +
  geom_errorbar(width = .2, position = position_dodge(.9), 
                color = "grey25", show.legend = FALSE) +
  ggthemes::theme_tufte() + 
  theme(panel.background = element_rect(color = "grey96", fill = "grey96"), 
        legend.position = "top", 
        legend.justification = "center",
        legend.title = element_blank(),
        strip.placement = "outside",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 12, face = "bold"),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 16, face = "bold"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, -5, 0),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text = element_text(size = 14, face = "bold")) +
  scale_fill_manual(values = c("red3", "turquoise3")) -> plot3

## Fig 4 Chelsa

bio %>%
  filter(Variable == "Annual range") %>%
  ggplot(aes(CHELSA, Loggers)) + 
  geom_point(size = 4, aes(color = Groundwater), alpha = 0.65) + 
  geom_smooth(method = "lm", se = F, linetype = "dashed", aes(color = Groundwater), show.legend = FALSE, alpha = 0.65) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  ggthemes::theme_tufte() + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
  xlab("CHELSA air temperatures (°C)") + ylab("Logger soil temperatures (°C)") +
  coord_fixed(ratio = 1, xlim = c(10, 30), ylim = c(10, 30), expand = TRUE, clip = "off") +
  theme(panel.background = element_rect(color = "grey96", fill = "grey96"), 
        legend.position = "none", 
        legend.justification = "center",
        legend.title = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(0, -0.9, 0, -1), "cm"),
        plot.title = element_text(margin = margin(t = 10, b = -20), face = "bold")) +
  labs(title = "  Annual range", x = NULL, y = NULL) +
  scale_color_manual(values = c("red3", "turquoise3")) -> plot4a

bio %>%
  filter(Variable == "Diurnal range") %>%
  ggplot(aes(CHELSA, Loggers)) + 
  geom_point(size = 4, aes(color = Groundwater), alpha = 0.65) + 
  geom_smooth(method = "lm", se = F, linetype = "dashed", aes(color = Groundwater), show.legend = FALSE, alpha = 0.65) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  ggthemes::theme_tufte() + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
  xlab("CHELSA air temperatures (°C)") + ylab("Logger soil temperatures (°C)") +
  coord_fixed(ratio = 1, xlim = c(0, 10), ylim = c(0, 10), expand = TRUE, clip = "off") +
  theme(panel.background = element_rect(color = "grey96", fill = "grey96"), 
        legend.position = "none", 
        legend.justification = "center",
        legend.title = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(0, 0, 0, -0.9), "cm"),
        plot.title = element_text(margin = margin(t = 10, b = -20), face = "bold")) +
  labs(title = "  Diurnal range", x = NULL, y = NULL) +
  scale_color_manual(values = c("red3", "turquoise3")) -> plot4b


bio %>%
  filter(Variable == "Summer max") %>%
  ggplot(aes(CHELSA, Loggers)) + 
  geom_point(size = 4, aes(color = Groundwater), alpha = 0.65) + 
  geom_smooth(method = "lm", se = F, linetype = "dashed", aes(color = Groundwater), show.legend = FALSE, alpha = 0.65) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  ggthemes::theme_tufte() + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
  xlab("CHELSA air temperatures (°C)") + ylab("Logger soil temperatures (°C)") +
  coord_fixed(ratio = 1, xlim = c(10, 35), ylim = c(10, 35), expand = TRUE, clip = "off") +
  theme(panel.background = element_rect(color = "grey96", fill = "grey96"), 
        legend.position = "none", 
        legend.justification = "center",
        legend.title = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(0, -0.9, 0, -1), "cm"),
        plot.title = element_text(margin = margin(t = 10, b = -20), face = "bold")) +
  labs(title = "  Summer max", x = NULL, y = NULL) +
  scale_color_manual(values = c("red3", "turquoise3")) -> plot4c

bio %>%
  filter(Variable == "Winter min") %>%
  ggplot(aes(CHELSA, Loggers)) + 
  geom_point(size = 4, aes(color = Groundwater), alpha = 0.65) + 
  geom_smooth(method = "lm", se = F, linetype = "dashed", aes(color = Groundwater), show.legend = FALSE, alpha = 0.65) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  ggthemes::theme_tufte() + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
  xlab("CHELSA air temperatures (°C)") + ylab("Logger soil temperatures (°C)") +
  coord_fixed(ratio = 1, xlim = c(-6.5, 6.5), ylim = c(-6.5, 6.5), expand = TRUE, clip = "off") +
  theme(panel.background = element_rect(color = "grey96", fill = "grey96"), 
        legend.position = "none", 
        legend.justification = "center",
        legend.title = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(0, 0, 0, -0.9), "cm"),
        plot.title = element_text(margin = margin(t = 10, b = -20), face = "bold")) +
  labs(title = "  Winter min", x = NULL, y = NULL) +
  scale_color_manual(values = c("red3", "turquoise3")) -> plot4d

y.grob <- grid::textGrob("Logger soil temperature (ºC)", 
                         gp = grid::gpar(fontface="bold", fontsize = 14), rot = 90)

x.grob <- grid::textGrob("CHELSA logger temperature (ºC)", 
                         gp = grid::gpar(fontface = "bold", fontsize = 14))

gridExtra::grid.arrange(gridExtra::arrangeGrob(plot4a, plot4b, plot4c, plot4d, left = y.grob, bottom = x.grob)) -> plot4z

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(plot4a + theme(legend.position = "top"))

gridExtra::grid.arrange(mylegend, plot4z, nrow = 2, heights = c(1, 10)) -> plot4

# Plot5 peak analysis

temperatures %>%
  group_by(Logger, Day = floor_date(Time, "day")) %>%
  summarise(Tmax = max(Temperature)) %>%
  merge(header) %>%
  filter(Groundwater == "Dry") %>%
  dplyr::select(Site, Day, Tmax) %>%
  group_by(Site) %>%
  filter(Tmax > quantile(Tmax, 0.95)) %>%
  dplyr::select(Site, Day) -> hotdays

temperatures %>%
  group_by(Logger, Day = floor_date(Time, "day")) %>%
  summarise(Tmax = max(Temperature)) %>%
  merge(header) %>%
  dplyr::select(Site, Groundwater, Day, Tmax) %>%
  merge(hotdays) %>%
  dplyr::select(Site, Groundwater, Tmax, Day) %>%
  spread(Groundwater, Tmax) %>%
  mutate(Buffer = Dry - Waterlogged) -> peakdf

peakdf %>% group_by(Site) %>%
  summarise(Dry = mean(Dry), Waterlogged = mean(Waterlogged)) -> meanpeaks

peakdf %>%
  ggplot(aes(Dry, Waterlogged, color = Site)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.65) +
  ggrepel::geom_label_repel(aes(label = Site), data = meanpeaks, fill = alpha(c("white"), 0.5), 
                            fontface = "bold", min.segment.length = 100, size = 2.5,
                            label.padding = 0.1) +
  coord_fixed(ratio = 1, xlim = c(12.5, 42.5), ylim = c(12.5, 42.5), expand = TRUE, clip = "off") +
  xlab("Dry point daily max (°C)") + ylab("Watelogged point daily max (°C)") +
  scale_color_manual(values = c("darkcyan", "darkorange", "dodgerblue3", "purple", "forestgreen", "gold3", "darkred", "grey15")) + 
  ggthemes::theme_tufte() + 
  theme(panel.background = element_rect(color = "grey96", fill = "grey96"), 
        legend.position = "none", 
        legend.justification = "center",
        legend.title = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) -> plot5

# Save

save(ttests, table1, table2, table3, file = "results/MSoutput.RData") # Save numerical results and tables to build manuscript

ggsave(plot1, file = "results/Fig1.png", 
       path = NULL, scale = 1, width = 170, height = 85, units = "mm", dpi = 300) # Save figures to build manuscript

ggsave(plot2, file = "results/Fig2.png", 
       path = NULL, scale = 1, width = 170, height = 220, units = "mm", dpi = 300)

ggsave(plot3, file = "results/Fig3.png", 
       path = NULL, scale = 1, width = 170, height = 85, units = "mm", dpi = 300)

ggsave(plot4, file = "results/Fig4.png", 
       path = NULL, scale = 1, width = 170, height = 170, units = "mm", dpi = 300)

ggsave(plot5, file = "results/Fig5.png", 
       path = NULL, scale = 1, width = 85, height = 85, units = "mm", dpi = 300)


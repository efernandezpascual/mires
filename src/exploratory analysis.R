library(tidyverse); library(lubridate)
load("results/clean_data/temperatures.RData")

# Days of recordings

temperatures %>%
  group_by(Logger) %>%
  summarise(Start = min(Time),
            End = max(Time),
            Length = End - Start) -> recording.period

# Bioclimatic variables US Geological Survey, WorlClim

# Bio1 Annual Mean Temperature

temperatures %>%
  group_by(Logger) %>%
  summarise(Annual.Mean.Temperature = mean(Temperature)) -> bio1

# Bio 2 Annual Mean Diurnal Range

temperatures %>%
  group_by(Logger, Day = floor_date(Time, "day")) %>%
  summarise(Tmax = max(Temperature), 
            Tmin = min(Temperature),
            Tdif = Tmax - Tmin) %>%
  group_by(Logger) %>%
  summarise(Mean.Diurnal.Range = mean(Tdif)) -> bio2

# Bio 4 Temperature Seasonality (Standard Deviation)

temperatures %>%
  group_by(Logger) %>%
  summarise(Temperature.Seasonality = sd(Temperature)) -> bio4

# Bio 5 Max Temperature of Warmest Month

temperatures %>%
  group_by(Logger, Day = floor_date(Time, "day")) %>%
  summarise(Tmax = max(Temperature)) %>%
  mutate(Month = month(Day)) %>%
  group_by(Logger, Month) %>%
  summarise(Max.Temperature.of.Warmest.Month = mean(Tmax)) %>%
  group_by(Logger) %>%
  filter(Max.Temperature.of.Warmest.Month == max(Max.Temperature.of.Warmest.Month)) %>%
  select(-Month) -> bio5

# Bio 6 Min Temperature of the Coldest Month

temperatures %>%
  group_by(Logger, Day = floor_date(Time, "day")) %>%
  summarise(Tmin = min(Temperature)) %>%
  mutate(Month = month(Day)) %>%
  group_by(Logger, Month) %>%
  summarise(Min.Temperature.of.Coldest.Month = mean(Tmin)) %>%
  group_by(Logger) %>%
  filter(Min.Temperature.of.Coldest.Month == min(Min.Temperature.of.Coldest.Month)) %>%
  select(-Month) -> bio6

# Merge and calculate Bio 7 Annual Temperature Range

bio1 %>%
  merge(bio2) %>%
  merge(bio4) %>%
  merge(bio5) %>%
  merge(bio6) %>%
  mutate(Temperature.Annual.Range = 
           Max.Temperature.of.Warmest.Month -
           Min.Temperature.of.Coldest.Month) -> bioclimatics

# Paired T-tests

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
  gather(Trait, Value, Annual.Mean.Temperature.y:Temperature.Annual.Range.y) %>%
  select(Site, Groundwater, Value, Trait, Value) %>%
  spread(Groundwater, Value) %>%
  filter(! Trait %in% "Min.Temperature.of.Coldest.Month.y") %>%
  group_by(Trait) %>%
  do(t.testEFP(.)),
merge(header, bioclimatics, by = "Logger") %>% 
  gather(Trait, Value, Annual.Mean.Temperature.y:Temperature.Annual.Range.y) %>%
  select(Site, Groundwater, Value, Trait, Value) %>%
  spread(Groundwater, Value) %>%
  filter(Trait %in% "Min.Temperature.of.Coldest.Month.y") %>%
  group_by(Trait) %>%
  filter(!Site %in% "Los Cándanos") %>%
  do(t.testEFPgreater(.)))

# Plot vs elevation

merge(header, bioclimatics) %>% 
  gather(Trait, Value, Annual.Mean.Temperature:Annual.Temperature.Range) %>%
  select(Elevation, Site, Groundwater, Value, Trait, Value) %>%
  spread(Groundwater, Value) %>%
  ggplot(aes(x = Elevation, y = Dry - Waterlogged)) +
  geom_text(aes(label = Site)) + 
  facet_wrap(~ Trait) + 
  geom_smooth(method = "lm")


# PCA

pcaEFP <- function(df, id = 1, header = NULL)
{
  z <- list()
  z$pca <- PCA(df[, -id], scale.unit = TRUE, graph = FALSE)
  if(is.null(header)) z$inds <- data.frame(Logger = df[, id], z$pca$ind$coord[, 1:2])
  if(! is.null(header)) z$inds <- merge(header, data.frame(Logger = df[, id], z$pca$ind$coord[, 1:2]))
  z$vars <- rownames_to_column(data.frame(z$pca$var$coord[, 1:2]), var = "Variable")
  z$vars$Variable <- gsub(".", " ", z$vars$Variable, fixed = TRUE)
  class(z) <- "mvEFP"
  z
}  

plot.mvEFP <- function(pca, colorby = NULL)
{
  colorby <- enquo(colorby)
  p <- ggplot(pca$inds, aes(x = Dim.1, y = Dim.2)) +
    geom_point(aes(color = !! colorby), size = 3, show.legend = T, alpha = 0.6) +
    #theme_classic2() + 
    theme(legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(size = 7.5),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
          plot.margin = unit(x = c(0, 0, 0, 0), units = "cm"),
          text = element_text(size = 8)) +
    
    scale_color_manual(values = c("forestgreen", "gold", "purple"))
  
  if(! is.null(pca$pca)) {p <- p + 
    scale_x_continuous(name = paste("Axis 1 (", round(pca$pca$eig[1, 2], 0),
                                    "% variance explained)", sep = "")) + 
    scale_y_continuous(name = paste("Axis 2 (", round(pca$pca$eig[2, 2], 0), 
                                    "% variance explained)", sep = "")) +
    #    geom_segment(data = pca$vars, aes(x = 0, y = 0, xend = 4*Dim.1, yend = 4*Dim.2), alpha = 0.4) +
    geom_text(data = pca$vars, aes(x = 4*Dim.1, y = 4*Dim.2, label = Variable), 
              alpha = 0.8, size = 3)} else {
                p <- p + 
                  scale_x_continuous(name = "Axis 1") + 
                  scale_y_continuous(name = "Axis 2") +
                  #    geom_segment(data = pca$vars, aes(x = 0, y = 0, xend = 4*Dim.1, yend = 4*Dim.2), alpha = 0.4) +
                  geom_text(data = pca$vars, 
                            aes(x = Dim.1, y = Dim.2, label = Variable), 
                            alpha = 0.8, size = 3)}
  p
}

summary.mvEFP <- function(pca)
{
  z <- data.frame(pca$pca$var$contrib[, 1:2], pca$pca$var$cor[, 1:2])
  is.num <- sapply(z, is.numeric)
  z[is.num] <- lapply(z[is.num], round, 2)
  colnames(z) <- c("cont.Dim.1", "cont.Dim.2", "corr.Dim.1", "corr.Dim.2")
  z
}

library(FactoMineR)

bioclimatics %>%
  pcaEFP(header = header) -> pcaBio
plot(pcaBio, colorby = Groundwater) + geom_text(aes(label = Site, color = Groundwater))

# Wodlclim comparison


bioclimatics %>%
  gather(Trait, Value, Annual.Mean.Temperature:Temperature.Annual.Range) %>%
  merge(header) %>%
  select(Site, Trait, Groundwater, Value) %>%
  spread(Groundwater, Value) %>%
  mutate(buffer = Dry - Waterlogged) %>%
  select(-c(Dry, Waterlogged)) -> buffer1

header %>%
  gather(Trait, Value, Annual.Mean.Temperature:Temperature.Annual.Range) %>%
  group_by(Site, Trait) %>%
  filter(Trait %in% buffer1$Trait) %>%
  filter(! Trait %in% c("Annual.Mean.Temperature", "Temperature.Seasonality")) %>% 
  summarise(worldclim = mean(Value)) %>% 
  merge(buffer1) %>%
  ggplot(aes(x = worldclim, y = buffer)) + 
  geom_point() +
  geom_text(aes(label = Site)) + 
  geom_smooth(method = "lm") + 
  facet_wrap(~ Trait, scales = "free")


header %>%
  gather(Trait, Value, Annual.Mean.Temperature:Temperature.Annual.Range) %>%
  group_by(Site, Trait) %>%
  filter(Trait %in% buffer1$Trait) %>%
  filter(! Trait %in% c("Annual.Mean.Temperature", "Temperature.Seasonality")) %>% 
  summarise(worldclim = mean(Value)) %>% 
  merge(buffer1) %>%
  filter(Trait == "Temperature.Annual.Range") -> df1

lm(buffer ~ worldclim, df1) -> m1

summary(m1)

header %>%
  gather(Trait, Value, Annual.Mean.Temperature:Temperature.Annual.Range) %>%
  group_by(Site, Trait) %>%
  filter(Trait %in% buffer1$Trait) %>%
  filter(! Trait %in% c("Annual.Mean.Temperature", "Temperature.Seasonality")) %>% 
  summarise(worldclim = mean(Value)) %>% 
  merge(buffer1) %>%
  filter(Trait == "Mean.Diurnal.Range") -> df1

lm(buffer ~ worldclim, df1) -> m1

summary(m1)

header %>%
  gather(Trait, Value, Annual.Mean.Temperature:Temperature.Annual.Range) %>%
  group_by(Site, Trait) %>%
  filter(Trait %in% buffer1$Trait) %>%
  filter(! Trait %in% c("Annual.Mean.Temperature", "Temperature.Seasonality")) %>% 
  summarise(worldclim = mean(Value)) %>% 
  merge(buffer1) %>%
  filter(Trait == "Min.Temperature.of.Coldest.Month") -> df1

lm(buffer ~ worldclim, df1) -> m1

summary(m1)

header %>%
  gather(Trait, Value, Annual.Mean.Temperature:Temperature.Annual.Range) %>%
  group_by(Site, Trait) %>%
  filter(Trait %in% buffer1$Trait) %>%
  filter(! Trait %in% c("Annual.Mean.Temperature", "Temperature.Seasonality")) %>% 
  summarise(worldclim = mean(Value)) %>% 
  merge(buffer1) %>%
  filter(Trait == "Max.Temperature.of.Warmest.Month") -> df1

lm(buffer ~ worldclim, df1) -> m1

summary(m1)

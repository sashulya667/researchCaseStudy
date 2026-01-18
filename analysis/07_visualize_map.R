# ============================================================
# COMPLETE, CLEAN, RUNNABLE SCRIPT
# ============================================================

# ----------------------------
# Libraries
# ----------------------------
library(sf)
library(dplyr)
library(ggplot2)
library(tidyr)

# ----------------------------
# Load spatial data
# ----------------------------
gpkg_path <- "/Users/sasha/Desktop/University/3rd semester/RCS/map/final_map.gpkg"

map_sf <- st_read(gpkg_path, layer = "output", quiet = TRUE)

# Tell sf which column is geometry
st_geometry(map_sf) <- "geom"

# ----------------------------
# Create PROV inside each REG
# ----------------------------
map_sf <- map_sf %>%
  group_by(REG) %>%
  mutate(PROV = row_number()) %>%
  ungroup()

# ----------------------------
# Plot 1: REGIONS
# ----------------------------
p_regions <- ggplot(map_sf) +
  geom_sf(
    aes(fill = factor(REG)),
    color = "white",
    linewidth = 0.5
  ) +
  theme_minimal() +
  labs(
    fill = "Region",
    title = "Regions"
  )

print(p_regions)

# ----------------------------
# Merge provinces in Region 2
# ----------------------------
reg2 <- map_sf %>% filter(REG == 2)
others <- map_sf %>% filter(REG != 2)

reg2_merged <- reg2 %>%
  mutate(
    PROV = case_when(
      PROV %in% c(1, 2) ~ 1,  # merge first two provinces
      TRUE ~ 2
    )
  ) %>%
  group_by(REG, PROV) %>%
  summarise(
    geom = st_union(geom),
    .groups = "drop"
  )

# Restore geometry column
st_geometry(reg2_merged) <- "geom"

# Recombine map
map_sf_fixed <- bind_rows(others, reg2_merged)

# ----------------------------
# Plot 2: PROVINCES
# ----------------------------
p_provinces <- ggplot(map_sf_fixed) +
  geom_sf(
    aes(fill = interaction(REG, PROV)),
    color = "white",
    linewidth = 0.4
  ) +
  theme_minimal() +
  labs(
    fill = "Region–Province",
    title = "Provinces (Region 2 merged)"
  )

print(p_provinces)

# ----------------------------
# Load Dataset_B
# ----------------------------
devtools::load_all()
library(here)
load_data <- function(filename) read.csv(file.path(here("data"), filename))

df_B <- load_data("Dataset_B.csv")

# ============================================================
# AGGREGATE DATA BY REGION
# ============================================================

# Mean EDI per region
edi_by_region <- df_B %>%
  group_by(REG) %>%
  summarise(
    mean_EDI = mean(EDI, na.rm = TRUE),
    .groups = "drop"
  )

# Number of distinct DOU per region (IMPORTANT: numeric)
dou_by_region <- df_B %>%
  group_by(REG) %>%
  summarise(
    n_DOU = n_distinct(DOU),
    .groups = "drop"
  )

# Household count per region
hh_by_region <- df_B %>%
  count(REG, name = "n_households")

# Combine all region indicators
region_stats <- edi_by_region %>%
  left_join(dou_by_region, by = "REG") %>%
  left_join(hh_by_region, by = "REG")

# ============================================================
# CREATE REGION GEOMETRY (FROM PROVINCES)
# ============================================================

regions_geom <- map_sf_fixed %>%
  group_by(REG) %>%
  summarise(
    geom = st_union(geom),
    .groups = "drop"
  )

st_geometry(regions_geom) <- "geom"

# Join statistics to geometry
regions_data <- regions_geom %>%
  left_join(region_stats, by = "REG")

# ============================================================
# PLOT 3: MEAN EDI BY REGION (MAP)
# ============================================================

p_edi_region <- ggplot(regions_data) +
  geom_sf(
    aes(fill = mean_EDI),
    color = "white",
    linewidth = 0.6
  ) +
  scale_fill_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(
    fill = "Mean EDI",
    title = "Mean EDI by Region"
  )

print(p_edi_region)

# ============================================================
# PLOT 4: NUMBER OF DOU BY REGION (MAP)
# ============================================================

p_dou_region <- ggplot(regions_data) +
  geom_sf(
    aes(fill = n_DOU),
    color = "white",
    linewidth = 0.6
  ) +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal() +
  labs(
    fill = "Number of DOU",
    title = "Administrative Units (DOU) by Region"
  )

print(p_dou_region)

# ============================================================
# BAR CHART: MEAN EDI BY REGION
# ============================================================

p_edi_bar <- ggplot(region_stats, aes(x = factor(REG), y = mean_EDI)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(
    x = "Region",
    y = "Mean EDI",
    title = "Mean EDI Comparison Across Regions"
  )

print(p_edi_bar)

# ============================================================
# AGGREGATE DATA BY PROVINCE (WITHIN REGION)
# ============================================================

province_stats <- df_B %>%
  group_by(REG, PROV) %>%
  summarise(
    mean_EDI = mean(EDI, na.rm = TRUE),
    n_DOU = n_distinct(DOU),
    n_households = n(),
    .groups = "drop"
  )

# Join province statistics to spatial data
province_data <- map_sf_fixed %>%
  left_join(province_stats, by = c("REG", "PROV"))

# ============================================================
# PLOT 5: MEAN EDI BY PROVINCE
# ============================================================

p_edi_province <- ggplot(province_data) +
  geom_sf(
    aes(fill = mean_EDI),
    color = "white",
    linewidth = 0.4
  ) +
  scale_fill_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(
    fill = "Mean EDI",
    title = "Mean EDI by Province"
  )

print(p_edi_province)



library("sp")
d = load("/Users/sasha/Downloads/AMAP_v0.2.1.RData")

plot(AMAP)

plot(AMAP, col = AMAP$REG, border = "white")

library(sf)
library(ggplot2)

AMAP_sf <- st_as_sf(AMAP)

names(AMAP_sf)

ggplot(AMAP_sf) +
  geom_sf(aes(fill = factor(REG)), color = "white", linewidth = 0.2) +
  theme_minimal() +
  labs(fill = "Region", title = "AMAP regions")


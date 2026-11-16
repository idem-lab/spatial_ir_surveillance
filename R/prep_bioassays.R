# Prepare the bioassay data: WHO discriminating-concentration pyrethroid
# bioassays on the An. gambiae complex from the Vector Atlas, as cleaned by the
# ir_cube repo: https://github.com/idem-lab/ir_cube, with a pixel ID and an
# admin-2 district for each.
# 
# Inputs from ir_cube.
# The data files are not in the repo, so we copy them here (gitignored) and
# record the ir_cube commit that made them, and their md5 sums. This script
# stops if a file's checksum differs from the stated checksums, so that everyone
# runs the analysis on the same data. If a change to ir_cube changes the data
# and it's copied in, it can't slip in unnoticed. To update the data
# deliberately, update the commit and checksums here:
#
#   data/ir_cube/all_gambiae_complex_data.RDS: the cleaned bioassays, made by
#     ir_cube's R/prep_bioassays.R (not this script, with the same name)
#   data/ir_cube/raster_mask.tif: the 0.0417 degree (about 5 km) grid of
#     pixels used by ir_cube, made by ir_cube's R/prep_rasters.R
#
# If these files are not in data/ir_cube/, the script copies them from a local
# clone of ir_cube at ir_cube_dir on Nick's machine. If you need these files,
# ask Nick to share them and put them in data/ir_cube/ yourself. Admin
# boundaries are downloaded from GADM 4.1 to data/gadm/.
#
# Output: data/bioassays.RDS

source("R/packages.R")

# ir_cube commit whose data prep made the copied files (idem-lab/ir_cube)
ir_cube_commit <- "2e7e4a62a35bc8905b9213f5e7384849b251cb9f"
ir_cube_dir <- "~/Dropbox/github/ir_cube"
ir_cube_files <- c(
  all_gambiae_complex_data.RDS = "44772bd1da61751c30d7e6e1f020eaf7",
  raster_mask.tif = "54e89d34918b62c0296cfec92c890585"
)

# WHO discriminating concentrations (%) of the pyrethroids ir_cube models.
# These are also the most common concentration of each in the data, which is
# how ir_cube subsets them.
discriminating_concentrations <- c(
  "Alpha-cypermethrin" = 0.05,
  "Deltamethrin" = 0.05,
  "Lambda-cyhalothrin" = 0.05,
  "Permethrin" = 0.75
)

# ir_cube drops bioassays before this year
baseline_year <- 1995

# copy the ir_cube files, if needed, and check they match the recorded ones
dir.create("data/ir_cube", showWarnings = FALSE, recursive = TRUE)
for (file in names(ir_cube_files)) {
  local_file <- file.path("data/ir_cube", file)
  if (!file.exists(local_file)) {
    source_file <- file.path(ir_cube_dir, "data/clean", file)
    if (!file.exists(source_file)) {
      stop(local_file, " is missing and there is no ir_cube clone at ",
           ir_cube_dir, ". Ask Nick for the file and put it in data/ir_cube/")
    }
    file.copy(source_file, local_file)
  }
  if (tools::md5sum(local_file) != ir_cube_files[[file]]) {
    stop(local_file, " does not match the file from ir_cube commit ",
         ir_cube_commit)
  }
}

mask <- terra::rast("data/ir_cube/raster_mask.tif")

bioassays <- readRDS("data/ir_cube/all_gambiae_complex_data.RDS") |>
  dplyr::filter(
    insecticide_type %in% names(discriminating_concentrations),
    concentration == discriminating_concentrations[insecticide_type],
    year >= baseline_year,
    # a handful of records of minor members of the complex
    !species %in% c("Anopheles merus", "Anopheles quadriannulatus")
  ) |>
  dplyr::mutate(
    pixel_id = terra::cellFromXY(mask, cbind(longitude, latitude))
  ) |>
  # drop records outside the mask, as ir_cube does
  dplyr::filter(
    !is.na(terra::extract(mask, pixel_id)[, 1])
  )

# ISO3 codes for the countries in the data (ir_cube renames a few)
country_iso3 <- geodata::country_codes() |>
  dplyr::select(country_name = NAME, iso3 = ISO3) |>
  dplyr::mutate(
    country_name = dplyr::case_when(
      country_name == "Central African Republic" ~ "CAR",
      country_name == "Côte d'Ivoire" ~ "Côte d’Ivoire",
      country_name == "Democratic Republic of the Congo" ~ "DR Congo",
      country_name == "São Tomé and Príncipe" ~ "Sao Tome & Principe",
      .default = country_name
    )
  ) |>
  dplyr::filter(country_name %in% bioassays$country_name)

missing_countries <- setdiff(bioassays$country_name, country_iso3$country_name)
if (length(missing_countries) > 0) {
  stop("no ISO3 code for: ", paste(missing_countries, collapse = ", "))
}

# admin-2 districts for these countries, or admin-1 where GADM has no admin-2
get_districts <- function(iso3) {
  for (level in 2:1) {
    # geodata saves to a gadm/ folder inside path
    polys <- geodata::gadm(iso3, level = level, path = "data",
                           version = "4.1")
    if (!is.null(polys)) {
      gid <- paste0("GID_", level)
      name <- paste0("NAME_", level)
      polys$district_id <- polys[[gid]][, 1]
      polys$district_name <- polys[[name]][, 1]
      polys$admin_level <- level
      return(polys[, c("district_id", "district_name", "admin_level")])
    }
  }
  stop("no admin-1 or admin-2 boundaries for ", iso3)
}

# base rbind() dispatches to terra's method for SpatVectors
districts <- do.call(rbind, lapply(country_iso3$iso3, get_districts)) |>
  sf::st_as_sf() |>
  # some GADM polygons have self-intersecting edges, which sf rejects
  sf::st_make_valid()

n_bioassays <- nrow(bioassays)

bioassays <- bioassays |>
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326,
               remove = FALSE) |>
  # the district containing each bioassay, or the nearest district if it is
  # outside every district (e.g. on the coast). A bioassay exactly on a shared
  # border goes to the first of the districts it touches
  sf::st_join(districts, join = sf::st_nearest_feature) |>
  # how far each bioassay is from its district: 0 if it's inside
  dplyr::mutate(
    district_distance_km = as.numeric(sf::st_distance(
      geometry,
      districts$geometry[match(district_id, districts$district_id)],
      by_element = TRUE
    )) / 1000
  ) |>
  sf::st_drop_geometry() |>
  dplyr::select(
    record_id,
    study,
    country_name,
    district_id,
    district_name,
    admin_level,
    district_distance_km,
    pixel_id,
    latitude, longitude,
    year, month,
    species,
    insecticide_type, concentration,
    mosquito_number, died,
    n_imputed
  )

stopifnot(
  nrow(bioassays) == n_bioassays,
  !anyNA(bioassays$district_id)
)

saveRDS(bioassays, "data/bioassays.RDS")

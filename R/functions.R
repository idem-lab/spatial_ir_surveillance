# Functions used by the analysis scripts. Source this file after
# R/packages.R.

# ==== kdr genotypes (R/prep_kdr_samples.R, R/fit_kdr_predictor.R) ====

# genotypes at Vgsc codon 995 (L is the wild type, F and S the kdr alleles),
# the number of kdr alleles each carries, and its count column in each source
# (IR Mapper's survivor/dead file uses the same names as its other files)
genotypes <- data.frame(
  genotype     = c("L/L", "L/F", "F/F", "L/S", "S/S", "F/S"),
  kdr_alleles  = c(0, 1, 2, 1, 2, 2),
  vector_atlas = c("vgsc995l_vgsc995l_n", "vgsc995l_vgsc995f_n",
                   "vgsc995f_vgsc995f_n", "vgsc995l_vgsc995s_n",
                   "vgsc995s_vgsc995s_n", "vgsc995f_vgsc995s_n"),
  ir_mapper    = c("L/L (no.)", "L/F (no.)", "F/F (no.)", "L/S (no.)",
                   "S/S (no.)", "F/S (no.)")
)

# ==== pixels (all scripts) ====

# the grid of pixels: the grid of ir_cube's raster mask (its origin and
# resolution) extended to the whole globe, without values, so that every
# bioassay and kdr sample is in a pixel, whether or not the mask covers it
pixel_grid <- function() {
  terra::extend(terra::rast(terra::rast("data/ir_cube/raster_mask.tif")),
                terra::ext(-180, 180, -90, 90), snap = "in")
}

# the ID of the pixel containing each location: its cell number in
# pixel_grid()
pixel_id_at <- function(longitude, latitude) {
  terra::cellFromXY(pixel_grid(), cbind(longitude, latitude))
}

# the centre of each pixel, a two-column matrix of longitude and latitude
pixel_centres <- function(pixel_id) {
  terra::xyFromCell(pixel_grid(), pixel_id)
}

# half-normal prior on a standard deviation in greta: N+(0, sd^2)
half_normal <- function(sd = 1, dim = 1) {
  greta::normal(0, sd, truncation = c(0, Inf), dim = dim)
}

# ==== bioassays (R/fit_kdr_predictor.R, R/fit_bioassay_variance.R) ====

# the insecticides of both models: the pyrethroids with enough survivor/dead
# pairs to estimate the kdr effect (alpha-cypermethrin has no pairs, and
# cyfluthrin and lambda-cyhalothrin one each); the variance model uses the
# same bioassays
insecticides <- c("Deltamethrin", "Permethrin")

# bioassay intra-class correlation of mortality by insecticide: the posterior
# means of ir_cube's R/fig_illustrate_bioassay_variability.R (outputs/
# bioassay_rho_hierarchical.csv, untracked) on its branch rho-study-effect
# (idem-lab/ir_cube#72), which fits rho hierarchically to all replicate
# bioassays in a pixel, year and insecticide, with a study-year effect, so
# that rho is the noise between bioassays of one study
rho <- setNames(c(0.1363, 0.1280), insecticides)

# beta-binomial distribution in greta with mean p and intra-class correlation
# rho (ir_cube's, from its R/functions.R)
betabinomial_p_rho <- function(N, p, rho) {
  a <- p * (1 / rho - 1)
  b <- a * (1 - p) / p
  greta::beta_binomial(size = N, alpha = a, beta = b)
}

# the likelihood of the bioassays' deaths: beta-binomial around their
# expected mortality, 1 - survival (a greta array), with rho fixed for each
# insecticide
bioassay_likelihood <- function(bioassays, survival) {
  rho_i <- unname(rho[bioassays$insecticide_type])
  stopifnot(!anyNA(rho_i))
  died <- greta::as_data(bioassays$died)
  greta::distribution(died) <- betabinomial_p_rho(bioassays$mosquito_number,
                                                  1 - survival, rho_i)
}

# a study left out: Yadouleton et al. (2010, Malaria Journal 9:83) report
# 100% deltamethrin mortality at all 43 localities despite kdr frequencies up
# to 0.97, with no susceptible reference strain, which isn't credible as a
# measure of population mortality
excluded_study <- "yadouleton 2010, Benin"

# bioassays with the insecticides, other than of the excluded study,
# labelled by pixel-year, district-year and study-year
select_bioassays <- function(bioassays) {
  bioassays |>
    dplyr::filter(insecticide_type %in% insecticides,
                  study != excluded_study) |>
    dplyr::mutate(pixel_year = paste(pixel_id, year),
                  district_year = paste(district_id, year),
                  study_year = paste(study, year))
}

# ==== kdr predictor (R/fit_kdr_predictor.R) ====

# pyrethroids as spelled in the IR Mapper survivor/dead file, and as in the
# bioassay data
pair_pyrethroids <- setNames(insecticides, tolower(insecticides))

# ---- kdr predictor: data ----

# kdr samples in the same pixel as a bioassay pixel-year, collected in that
# year, one row per (pixel-year, sample) pair. A sample collected over
# several years that match bioassays in the same pixel is assigned to the one
# nearest the middle of its collection years (the earlier on a tie), so each
# sample is used once
find_same_pixel_samples <- function(bioassays, kdr_samples) {
  pixel_years <- dplyr::distinct(bioassays, pixel_year, pixel_id, year)
  kdr_samples |>
    dplyr::inner_join(pixel_years, by = dplyr::join_by(pixel_id),
                      relationship = "many-to-many") |>
    dplyr::filter(year_start <= year, year_end >= year) |>
    dplyr::group_by(record_id) |>
    dplyr::arrange(abs(year - (year_start + year_end) / 2), year) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()
}

# bioassays in pixel-years with kdr samples linked to them (links: one row per
# pixel-year and sample). selection is the pooled kdr frequency of those
# samples (kdr alleles over alleles genotyped) in three bands, for grouping
# posterior predictive checks
select_kdr_sampled <- function(bioassays, links) {
  sampled <- links |>
    dplyr::group_by(pixel_year) |>
    dplyr::summarise(kdr_frequency = sum(round(2 * n_tested * kdr_frequency)) /
                       sum(2 * n_tested)) |>
    dplyr::mutate(selection = cut(kdr_frequency, c(-Inf, 0.3, 0.7, Inf),
                                  c("kdr below 0.3", "kdr 0.3 to 0.7",
                                    "kdr 0.7 or above"),
                                  right = FALSE))
  dplyr::inner_join(bioassays, dplyr::select(sampled, pixel_year, selection),
                    by = dplyr::join_by(pixel_year))
}

# genotype counts of survivors and dead in IR Mapper's survivor/dead pairs
# tested with deltamethrin or permethrin, other than An. arabiensis (left
# out of the bioassays and kdr samples too, in R/prep_bioassays.R and
# R/prep_kdr_samples.R, as its kdr effect differs): one row
# per pair and number of kdr alleles a, with the numbers of survivors and dead
# genotyped. Where a row gives genotype percentages but no counts, counts are
# the percentages times the number genotyped
read_pair_genotypes <- function(file) {
  # missing values are blank, NA or NR
  rows <- read.csv(file, check.names = FALSE, fileEncoding = "latin1",
                   na.strings = c("", "NA", "NR"))
  # some column names end in (non-breaking) spaces
  names(rows) <- trimws(gsub("\u00a0", " ", names(rows)))
  counts <- as.matrix(rows[genotypes$ir_mapper])
  percents <- as.matrix(rows[paste(genotypes$genotype, "(%)")])
  n_tested <- rows[["No. mosquitoes tested"]]
  has_counts <- rowSums(!is.na(counts)) > 0
  counts[is.na(counts)] <- 0
  from_percents <- round(percents * n_tested / 100)
  from_percents[is.na(from_percents)] <- 0
  counts[!has_counts, ] <- from_percents[!has_counts, ]
  stopifnot(all(rows[["Subset of susceptibility tests"]] %in%
                  c("survivors", "dead")))

  # one row per subset (survivors or dead) and number of kdr alleles
  long <- data.frame(
    pair = rows[["Matched pair ID"]],
    subset = rows[["Subset of susceptibility tests"]],
    species = tolower(rows[["Anophelines tested"]]),
    insecticide = tolower(rows[["Insecticide tested"]]),
    a = rep(0:2, each = nrow(rows)),
    n = c(sapply(0:2, function(a) {
      rowSums(counts[, genotypes$kdr_alleles == a, drop = FALSE])
    }))
  ) |>
    dplyr::filter(insecticide %in% names(pair_pyrethroids),
                  !grepl("arabiensis", species)) |>
    dplyr::mutate(insecticide = unname(pair_pyrethroids[insecticide]))
  stopifnot(!anyDuplicated(long[, c("pair", "subset", "a")]))

  survivors <- long |>
    dplyr::filter(subset == "survivors") |>
    dplyr::select(pair, insecticide, a, survivors = n)
  dead <- long |>
    dplyr::filter(subset == "dead") |>
    dplyr::select(pair, a, dead = n)
  dplyr::inner_join(survivors, dead, by = dplyr::join_by(pair, a)) |>
    # pairs with mosquitoes of only one genotype, or with no survivors or no
    # dead, carry no information about the kdr effect
    dplyr::group_by(pair) |>
    dplyr::filter(sum(survivors + dead > 0) > 1,
                  sum(survivors) > 0, sum(dead) > 0) |>
    dplyr::ungroup() |>
    dplyr::filter(survivors + dead > 0)
}

# ---- kdr predictor: model ----

# log hazard of a fully susceptible population in a bioassay: WHO set the
# discriminating concentration at twice the LC99.9 of a susceptible colony,
# and under the hazard model doubling the concentration doubles the hazard,
# -log(0.001) at the LC99.9
alpha_star <- log(-2 * log(0.001))

# hazard of dying in a bioassay, H, for a population with kdr frequency f,
# kdr effect per allele kdr_effect, and the other drivers' frequency g,
# ilogit(eta), and effect gamma (greta arrays or numbers)
hazard <- function(f, kdr_effect, g, gamma) {
  exp(alpha_star - 2 * f * kdr_effect - 2 * g * gamma)
}

# bioassay survival, exp(-H), and mortality, 1 - exp(-H), from the hazard H
# (mortality computed so that it is precise when H is small)
survival <- function(H) {
  exp(-H)
}
mortality <- function(H) {
  -expm1(-H)
}

# greta model of the kdr predictor, from the selected bioassays, their kdr
# samples (from find_same_pixel_samples()) and the survivor/dead pairs (from
# read_pair_genotypes()). For bioassay i of insecticide j in pixel-year
# (l, t), the log hazard falls with the mean number of kdr alleles, 2f, and
# with the other drivers of resistance (other mechanisms and bioassay
# conditions), in the same form: an allele frequency, ilogit(eta), times a
# positive effect per allele, gamma_j:
#   log H_i = alpha* - 2 f_lt exp(lambda_j + b_lt) - 2 gamma_j ilogit(eta_lt)
#   1 - pi_i = exp(-H_i)
#   logit(f_lt) = mu_f + sigma_f z_f_lt,  eta_lt = sigma_eta z_eta_lt
# where z_f and z_eta are standard normal with correlation psi (co-selection:
# places with more kdr also have more of the other resistance). Returns the
# model; the greta array of each bioassay's expected survival, 1 - pi_i (for
# posterior predictive checks); the insecticides, in the order of lambda and
# gamma; the pixel-years, in the order of the pixel-years' latent variables;
# and the greta arrays of those latent variables, b, eta and f
build_kdr_model <- function(bioassays, samples, pairs) {
  # ---- indices ----
  pixel_years <- unique(bioassays$pixel_year)
  insecticides <- sort(unique(bioassays$insecticide_type))
  pair_ids <- unique(pairs$pair)
  lt_i <- match(bioassays$pixel_year, pixel_years)
  j_i <- match(bioassays$insecticide_type, insecticides)
  lt_k <- match(samples$pixel_year, pixel_years)
  h <- match(pairs$pair, pair_ids)
  j_h <- match(pairs$insecticide, insecticides)
  stopifnot(!anyNA(lt_k))
  # the pairs inform the kdr effects of the bioassays' insecticides
  stopifnot(all(pairs$insecticide %in% insecticides))
  n_lt <- length(pixel_years)

  # ---- priors ----
  # log of the median kdr effect per allele, and the effect of the other
  # drivers, for each insecticide
  lambda <- greta::normal(0, 1, dim = length(insecticides))
  gamma <- greta::lognormal(0, 1, dim = length(insecticides))
  sigma_b <- half_normal()
  mu_f <- greta::normal(0, 3)
  sigma_f <- half_normal(2)
  sigma_eta <- half_normal()
  psi <- greta::uniform(-1, 1)

  # ---- kdr frequency ----
  z_f <- greta::normal(0, 1, dim = n_lt)
  f <- greta::ilogit(mu_f + sigma_f * z_f)
  # kdr alleles A_k in each sample k of M_k mosquitoes: beta-binomial around
  # its pixel-year's frequency, with intra-class correlation rho_kdr (uniform
  # prior)
  rho_kdr <- greta::beta(1, 1)
  f_k <- f[lt_k]
  alleles <- 2 * samples$n_tested
  kdr_alleles <- greta::as_data(round(alleles * samples$kdr_frequency))
  greta::distribution(kdr_alleles) <- betabinomial_p_rho(alleles, f_k,
                                                         rho_kdr)

  # ---- bioassays ----
  b <- sigma_b * greta::normal(0, 1, dim = n_lt)
  eta <- sigma_eta * (psi * z_f + sqrt(1 - psi^2) *
                        greta::normal(0, 1, dim = n_lt))
  survival <- survival(hazard(f[lt_i], exp(lambda[j_i] + b[lt_i]),
                              greta::ilogit(eta[lt_i]), gamma[j_i]))
  bioassay_likelihood(bioassays, survival)

  # ---- survivor/dead pairs ----
  # genotyped survivors L of L + D with a kdr alleles in pair h: binomial with
  # logit(theta_ha) = kappa_h + logit(S_ha), where S_ha is the survival of
  # mosquitoes with a kdr alleles, exp(-exp(omega_h - exp(lambda_j + b'_h) a)),
  # kappa_h a nuisance intercept (for the shares of survivors and dead
  # genotyped), omega_h a nuisance log hazard without kdr, and b'_h the
  # pair's deviation in the kdr effect
  kappa_pair <- greta::normal(0, 3, dim = length(pair_ids))
  omega_pair <- greta::normal(0, 3, dim = length(pair_ids))
  b_pair <- sigma_b * greta::normal(0, 1, dim = length(pair_ids))
  # logit(S) = -x - log(1 - exp(-x)), with x = -log(S), computed stably
  x <- exp(omega_pair[h] - exp(lambda[j_h] + b_pair[h]) * pairs$a)
  theta <- greta::ilogit(kappa_pair[h] - x - log(-expm1(-x)))
  survivors <- greta::as_data(pairs$survivors)
  greta::distribution(survivors) <- greta::binomial(
    pairs$survivors + pairs$dead, theta
  )

  list(model = greta::model(lambda, gamma, sigma_b, mu_f, sigma_f, sigma_eta,
                            psi, rho_kdr),
       survival = survival, insecticides = insecticides,
       pixel_years = pixel_years, b = b, eta = eta, f = f)
}

# the median kdr effect per allele, exp(lambda), and bioassay mortality for
# each insecticide, for a population with the median kdr effect (b = 0) and
# the other drivers at their median (eta = 0): without kdr, and with kdr
# fixed (f = 1), summarised as summarise_matrix()
kdr_mortality <- function(draws, insecticides) {
  matrix_draws <- as.matrix(draws)
  lambda <- matrix_draws[, paste0("lambda[", seq_along(insecticides), ",1]"),
                         drop = FALSE]
  gamma <- matrix_draws[, paste0("gamma[", seq_along(insecticides), ",1]"),
                        drop = FALSE]
  estimates <- cbind(exp(lambda),
                     mortality(hazard(0, exp(lambda), plogis(0), gamma)),
                     mortality(hazard(1, exp(lambda), plogis(0), gamma)))
  colnames(estimates) <- c(paste0("exp_lambda_", insecticides),
                           paste0("mortality_kdr_absent_", insecticides),
                           paste0("mortality_kdr_fixed_", insecticides))
  summarise_matrix(estimates)
}

# ==== bioassay variance (R/fit_bioassay_variance.R) ====

# rate of the exponential prior on 1 / range of the pixel-year effects:
# P(range < 10 km) = 0.05
inv_range_prior_rate <- -10 * log(0.05)

# Wendland's taper psi(r) = (1 - r)^4 (1 + 4 r) of distances d, with r =
# d / taper_km, 0 beyond taper_km. With great-circle distances, the tapered
# exponential correlation is positive definite on the sphere (Gneiting 2013)
wendland_taper <- function(d, taper_km) {
  r <- pmin(d / taper_km, 1)
  (1 - r)^4 * (1 + 4 * r)
}

# ---- bioassay variance: model ----

# greta model of bioassay variance, from the selected bioassays. For
# bioassay i of insecticide j in pixel l, district d, province r, country c
# and year t:
#   logit(1 - pi_i) = W_jlt + B_jdt
# The within-district-year effect is
#   W_jlt = sigma_v_j v_lt + sigma_e_j e_lt
# with v_.t a Gaussian process over all pixels in year t, with unit variance
# and the correlation exp(-d / range) psi(d / taper_km) between pixels a
# great-circle distance d apart, and e_lt independent standard normal (a
# nugget). The between-district-year effect is
#   B_jdt = alpha_j + (beta_j + sigma_trend_j trend_c) (t - 2015) / 10 +
#           sigma_country_j b_c + sigma_country_year_j h_ct +
#           sigma_province_j p_r + sigma_province_year_j q_rt +
#           sigma_district_j c_d + sigma_district_year_j delta_dt
# The standardised effects (v, e, trend, b, h, p, q, c, delta) are standard
# normal and shared between the insecticides, whose intercepts, trends and SDs
# differ: each effect is the same up to its scale in both insecticides. The
# country effect is centred, as the data identify it well: the first
# insecticide's effect has SD sigma_country_1, and the others' are it scaled
# by their SDs; the rest are non-centred. Provinces are admin-1 units. The
# prior on 1 / range shrinks toward an infinite range (Fuglstad et al. 2019):
# exponential with rate inv_range_prior_rate. pixel_xy is a two-column matrix
# of each bioassay's pixel centre (longitude, latitude). Returns the model,
# the insecticides (in the order of the parameters), and what
# fresh_survival() needs
build_variance_model <- function(bioassays, pixel_xy, taper_km) {
  # ---- indices ----
  insecticides <- sort(unique(bioassays$insecticide_type))
  j_i <- match(bioassays$insecticide_type, insecticides)
  # the index of each bioassay's level of an effect
  level <- function(label) match(label, unique(label))
  # provinces: the admin-1 unit of each admin-2 district (GADM ids
  # "BFA.9.1_1" -> "BFA.9"); admin-1 districts are their own province
  province_id <- ifelse(bioassays$admin_level == 2,
                        sub("\\.[0-9]+_[0-9]+$", "", bioassays$district_id),
                        bioassays$district_id)
  country_i <- level(bioassays$country_name)
  lt_i <- level(bioassays$pixel_year)
  first_lt <- match(seq_len(max(lt_i)), lt_i)
  clusters <- pixel_year_clusters(pixel_xy[first_lt, , drop = FALSE],
                                  bioassays$year[first_lt], taper_km)
  decade <- (bioassays$year - 2015) / 10

  # ---- priors ----
  alpha <- greta::normal(0, 3, dim = length(insecticides))
  beta <- greta::normal(0, 1, dim = length(insecticides))
  sigma_trend <- half_normal(dim = length(insecticides))
  sigma_country <- half_normal(dim = length(insecticides))
  sigma_country_year <- half_normal(dim = length(insecticides))
  sigma_province <- half_normal(dim = length(insecticides))
  sigma_province_year <- half_normal(dim = length(insecticides))
  sigma_district <- half_normal(dim = length(insecticides))
  sigma_district_year <- half_normal(dim = length(insecticides))
  sigma_v <- half_normal(dim = length(insecticides))
  sigma_e <- half_normal(dim = length(insecticides))
  inv_range <- greta::exponential(inv_range_prior_rate)

  # ---- effects ----
  # a non-centred effect of the levels in index: a standardised effect
  # shared between the insecticides, scaled by each insecticide's SD
  effect <- function(sigma, index) {
    sigma[j_i] * greta::normal(0, 1, dim = max(index))[index]
  }
  country_1 <- greta::normal(0, sigma_country[1], dim = max(country_i))
  country <- (sigma_country[j_i] / sigma_country[1]) * country_1[country_i]
  between <- alpha[j_i] +
    (beta[j_i] + effect(sigma_trend, country_i)) * decade +
    country +
    effect(sigma_country_year,
           level(paste(bioassays$country_name, bioassays$year))) +
    effect(sigma_province, level(province_id)) +
    effect(sigma_province_year, level(paste(province_id, bioassays$year))) +
    effect(sigma_district, level(bioassays$district_id)) +
    effect(sigma_district_year, level(bioassays$district_year))
  within <- within_effect(sigma_v, sigma_e, inv_range, j_i, lt_i, clusters,
                          taper_km)

  # ---- bioassays ----
  survival <- greta::ilogit(between + within)
  bioassay_likelihood(bioassays, survival)
  model <- greta::model(alpha, beta, sigma_trend, sigma_country,
                        sigma_country_year, sigma_province,
                        sigma_province_year, sigma_district,
                        sigma_district_year, sigma_v, sigma_e, inv_range)

  list(model = model, insecticides = insecticides, between = between,
       sigma_v = sigma_v, sigma_e = sigma_e, inv_range = inv_range,
       j_i = j_i, lt_i = lt_i, clusters = clusters, taper_km = taper_km)
}

# the within-district-year effect of each bioassay: new pixel-year effects v
# (a GP over the points of clusters from pixel_year_clusters()) and e, shared
# between the insecticides and scaled by each insecticide's SDs. j_i and lt_i
# are each bioassay's insecticide and pixel-year
within_effect <- function(sigma_v, sigma_e, inv_range, j_i, lt_i, clusters,
                          taper_km) {
  v <- tapered_gp(clusters, inv_range, taper_km)
  e <- greta::normal(0, 1, dim = length(v))
  sigma_v[j_i] * v[lt_i] + sigma_e[j_i] * e[lt_i]
}

# expected survival of each bioassay with new pixel-year effects (v and e)
# drawn from their fitted SDs and range, and the between-district-year
# effects at their fitted values: new bioassays at the fitted pixels of the
# same district-years. fit is the list from build_variance_model() (or as
# saved by R/fit_bioassay_variance.R); greta::calculate() with the posterior
# draws as values and nsim draws the new effects for each posterior draw.
# Made outside build_variance_model() because greta::model() takes in every
# variable connected to its parameters, and would sample these as free
# parameters
fresh_survival <- function(fit) {
  greta::ilogit(fit$between + within_effect(fit$sigma_v, fit$sigma_e,
                                            fit$inv_range, fit$j_i, fit$lt_i,
                                            fit$clusters, fit$taper_km))
}

# clusters of points (rows of the two-column matrix xy, of longitude and
# latitude) within the same group, linked by chains of points under taper_km
# apart (all the points of a group without a taper). Points in different
# clusters are independent, so the correlation matrix is block-diagonal over
# clusters. Returns the clusters in classes by size, each class with: the
# points of each cluster (a matrix of clusters by slots, with NA in unused
# slots after the used ones), and the great-circle distances (km) between
# them (an array of clusters by slots by slots). Unused slots are given a huge
# distance from every other slot, so each cluster's correlation matrix is
# padded to the class size with an identity and stays positive definite; as
# they come after the used slots, they do not change the used slots'
# Cholesky factor
pixel_year_clusters <- function(xy, group, taper_km) {
  points <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]),
                         coords = c("x", "y"), crs = 4326)
  members <- unlist(lapply(split(seq_along(group), group), function(rows) {
    if (is.infinite(taper_km)) return(list(rows))
    near <- sf::st_is_within_distance(points[rows, ], dist = taper_km * 1000)
    graph <- igraph::graph_from_adj_list(near, mode = "all")
    split(rows, igraph::components(graph)$membership)
  }), recursive = FALSE)
  # clusters are padded to the largest in their size class, so the classes
  # keep the padding small while keeping the number of batches low
  size_class <- cut(lengths(members),
                    c(0, 10, 25, 50, 100, 150, 200, 300, Inf))
  lapply(split(members, size_class, drop = TRUE), function(class_members) {
    slots <- max(lengths(class_members))
    point <- t(vapply(class_members,
                      function(m) c(m, rep(NA, slots - length(m))),
                      numeric(slots)))
    padding <- matrix(1e6, slots, slots)
    diag(padding) <- 0
    distance <- array(rep(padding, each = length(class_members)),
                      c(length(class_members), slots, slots))
    for (k in seq_along(class_members)) {
      m <- class_members[[k]]
      distance[k, seq_along(m), seq_along(m)] <-
        as.numeric(sf::st_distance(points[m, ])) / 1000
    }
    list(point = point, distance = distance)
  })
}

# unit-variance pixel-year effects with the tapered exponential correlation
# of 1 / range (a scalar greta array) between the points of each cluster from
# pixel_year_clusters(): for each class, one batched Cholesky factor per
# cluster of its padded correlation matrix times new standard normals
# (tf_tapered_gp()). greta has no public interface for custom operations, so
# this uses its internal constructor (greta from GitHub, as in renv.lock)
tapered_gp <- function(clusters, inv_range, taper_km) {
  parts <- lapply(clusters, function(size_class) {
    used <- which(!is.na(size_class$point))
    # standard normals only for the used slots, which tf_tapered_gp() places
    # by slot (0-based positions in the normals, with the last, one past them,
    # a zero for the unused slots), so the unused slots add nothing to the
    # free state
    slot <- array(length(used), dim(size_class$point))
    slot[used] <- seq_along(used) - 1L
    effects <- greta::.internals$nodes$constructors$op(
      "tapered_gp", inv_range,
      greta::normal(0, 1, dim = length(used)),
      greta::as_data(size_class$distance),
      dim = dim(size_class$point),
      operation_args = list(taper_km = taper_km, slot = slot),
      tf_operation = "tf_tapered_gp", tf_function_env = environment()
    )
    list(effects = effects[used], point = size_class$point[used])
  })
  effects <- do.call(c, lapply(parts, `[[`, "effects"))
  effects[order(unlist(lapply(parts, `[[`, "point")))]
}

# tensorflow of tapered_gp(): inv_range (draws, 1, 1), z (draws, used slots,
# 1) and distance (1, clusters, slots, slots; broadcast over draws) tensors,
# and slot, the position in z of each cluster's slots. The taper is computed
# as in wendland_taper()
tf_tapered_gp <- function(inv_range, z, distance, taper_km, slot) {
  tf <- tensorflow::tf
  z <- tf$squeeze(tf$pad(z, list(c(0L, 0L), c(0L, 1L), c(0L, 0L))), -1L)
  z <- tf$gather(z, slot, axis = 1L)
  r <- tf$minimum(distance / taper_km, 1)
  taper <- (1 - r)^4 * (1 + 4 * r)
  correlation <- tf$exp(-distance * tf$expand_dims(inv_range, -1L)) * taper
  chol <- tf$linalg$cholesky(correlation)
  tf$squeeze(tf$matmul(chol, tf$expand_dims(z, -1L)), -1L)
}

# summaries of bioassay mortality for the checks of the variance model: those
# of mortality_summaries(), and the mean over district-years of the SD and of
# the variance of mortality between their pixel-years (the mean SD is smaller
# than the mean variance implies if the variance differs between
# district-years)
variance_summaries <- function(mortality, bioassays) {
  sds <- sapply(split(seq_along(mortality), bioassays$district_year),
                function(rows) {
                  sd(tapply(mortality[rows], bioassays$pixel_year[rows], mean))
                })
  c(mortality_summaries(mortality, bioassays),
    within_district_sd = mean(sds), within_district_var = mean(sds^2))
}

# ---- bioassay variance: simulated districts ----

# the archetypal district: a circle 100 km across (about the median diameter
# of districts with bioassay data), on a grid of 5 km pixels (as the
# bioassay data; centres in km)
district_grid <- function() {
  centres <- seq(-47.5, 47.5, by = 5)
  grid <- expand.grid(x = centres, y = centres)
  grid[sqrt(grid$x^2 + grid$y^2) <= 50, ]
}

# one realisation of a unit-variance GP on grid with the tapered exponential
# correlation of length scale range
simulate_gp <- function(grid, range, taper_km) {
  distance <- as.matrix(dist(grid))
  correlation <- exp(-distance / range) * wendland_taper(distance, taper_km)
  as.vector(t(chol(correlation + diag(1e-9, nrow(grid)))) %*%
              rnorm(nrow(grid)))
}

# mortality in each pixel of a district with pixel effects u (logit scale)
# whose mean mortality over its pixels is district_mean: the logit-scale
# district effect is solved for this realisation, so every simulated district
# has the stated mean
district_mortality <- function(u, district_mean) {
  shift <- uniroot(function(mu) mean(plogis(mu + u)) - district_mean,
                   qlogis(district_mean) + c(-1, 1) * (10 + max(abs(u))),
                   tol = 1e-10)$root
  plogis(shift + u)
}

# column labels of simulated districts: the parameters of each draw, one per
# line (markdown, for ggtext)
district_labels <- function(sigma_v, range, sigma_e = NULL) {
  label <- sprintf("\u03c3<sub>v</sub> = %.2f<br>\u2113 = %.0f km", sigma_v,
                   range)
  if (is.null(sigma_e)) return(label)
  paste0(label, sprintf("<br>\u03c3<sub>e</sub> = %.2f", sigma_e))
}

# mortality on the district grid for each of the pixel effects in effects (a
# list, one per draw, labelled by labels) and each district mean (labelled
# by row_labels), with the draws ordered by increasing variance of mortality
# between pixels, averaged over the district means
simulate_districts <- function(effects, labels, district_means,
                               row_labels = sprintf("Mean %.0f%%",
                                                    100 * district_means)) {
  variance <- sapply(effects, function(u) {
    mean(sapply(district_means, function(m) var(district_mortality(u, m))))
  })
  panels <- expand.grid(k = seq_along(effects), r = seq_along(district_means))
  do.call(rbind, Map(function(k, r) {
    data.frame(district_grid(),
               mortality = district_mortality(effects[[k]], district_means[r]),
               draw = factor(labels[k], labels[order(variance)]),
               district_mean = factor(row_labels[r], row_labels))
  }, panels$k, panels$r))
}

# the simulated districts plotted: mortality on the grid (fill), by district
# mean (rows) and parameter draw (columns), white to colour
plot_districts <- function(simulated, title, colour) {
  ggplot2::ggplot(simulated, ggplot2::aes(x, y, fill = mortality)) +
    ggplot2::geom_raster() +
    ggplot2::coord_equal() +
    ggplot2::facet_grid(district_mean ~ draw, switch = "both") +
    ggplot2::scale_fill_gradient(low = "white", high = colour,
                                 limits = c(0, 1),
                                 labels = scales::percent,
                                 name = "Mortality") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 18) +
    ggplot2::theme(strip.text.x = ggtext::element_markdown(size = 16,
                                                           lineheight = 1.3),
                   strip.text.y.left = ggplot2::element_text(size = 16,
                                                             angle = 90),
                   plot.title = ggplot2::element_text(size = 22),
                   plot.margin = ggplot2::margin(10, 10, 10, 10))
}

# simulated districts drawn from the priors of the variance model: n_draws
# draws of sigma_v, the range and sigma_e, with the pixel effects with and
# without e (the same realisation of v in both, and the columns in the same
# order: by increasing variance of mortality between pixels without e).
# Returns the two plots, white to a muted grey-brown, distinct from ir_cube's
# insecticide colours
prior_district_plots <- function(district_means, taper_km, n_draws) {
  grid <- district_grid()
  sigma_v <- abs(rnorm(n_draws))
  sigma_e <- abs(rnorm(n_draws))
  range <- 1 / rexp(n_draws, inv_range_prior_rate)
  v <- lapply(seq_len(n_draws), function(k) {
    sigma_v[k] * simulate_gp(grid, range[k], taper_km)
  })
  v_e <- lapply(seq_len(n_draws), function(k) {
    v[[k]] + sigma_e[k] * rnorm(nrow(grid))
  })
  labels_v <- district_labels(sigma_v, range)
  labels_v_e <- district_labels(sigma_v, range, sigma_e)
  without_e <- simulate_districts(v, labels_v, district_means)
  with_e <- simulate_districts(v_e, labels_v_e, district_means)
  with_e$draw <- factor(with_e$draw,
                        labels_v_e[match(levels(without_e$draw), labels_v)])
  title <- "Districts simulated from the model priors,"
  colour <- "#6E5A4B"
  list(with_e = plot_districts(with_e, paste(title, "with pixel effects"),
                               colour),
       without_e = plot_districts(without_e,
                                  paste(title, "without pixel effects"),
                                  colour))
}

# the pixel effects of posterior draw k (a row of the matrix of posterior
# draws) on the district grid, for each insecticide: one realisation of the
# standardised v and e, scaled by each insecticide's SDs
posterior_pixel_effects <- function(posterior, k, insecticides, taper_km) {
  grid <- district_grid()
  v <- simulate_gp(grid, 1 / posterior[k, "inv_range"], taper_km)
  e <- rnorm(nrow(grid))
  lapply(seq_along(insecticides), function(j) {
    posterior[k, sprintf("sigma_v[%d,1]", j)] * v +
      posterior[k, sprintf("sigma_e[%d,1]", j)] * e
  })
}

# the SD of mortality (percentage points) between the pixels of districts
# with each district mean, over n_draws posterior draws, each with its own
# realisation of the pixel effects: mean, median and 95% interval, by
# insecticide and district mean
district_sd_summary <- function(posterior, insecticides, district_means,
                                taper_km, n_draws) {
  sds <- do.call(rbind, lapply(sample(nrow(posterior), n_draws), function(k) {
    effects <- posterior_pixel_effects(posterior, k, insecticides, taper_km)
    do.call(rbind, lapply(seq_along(insecticides), function(j) {
      data.frame(insecticide = insecticides[j], district_mean = district_means,
                 sd = sapply(district_means, function(m) {
                   100 * sd(district_mortality(effects[[j]], m))
                 }))
    }))
  }))
  sds |>
    dplyr::group_by(insecticide, district_mean) |>
    dplyr::summarise(mean = mean(sd), median = median(sd),
                     lower = quantile(sd, 0.025), upper = quantile(sd, 0.975),
                     .groups = "drop")
}

# simulated districts drawn from the posterior of the variance model, for
# each insecticide: n_draws posterior draws, with rows labelled by the mean
# SD of mortality between pixels (from district_sd_summary()). Returns the
# plots, named by insecticide, white to ir_cube's pyrethroid blue
posterior_district_plots <- function(posterior, insecticides, district_means,
                                     district_sds, taper_km, n_draws) {
  draws <- sample(nrow(posterior), n_draws)
  effects <- lapply(draws, posterior_pixel_effects, posterior = posterior,
                    insecticides = insecticides, taper_km = taper_km)
  plots <- lapply(seq_along(insecticides), function(j) {
    sds <- district_sds[district_sds$insecticide == insecticides[j], ]
    sds <- sds[match(district_means, sds$district_mean), ]
    labels <- district_labels(posterior[draws, sprintf("sigma_v[%d,1]", j)],
                              1 / posterior[draws, "inv_range"],
                              posterior[draws, sprintf("sigma_e[%d,1]", j)])
    simulated <- simulate_districts(
      lapply(effects, `[[`, j), labels, district_means,
      sprintf("Mean %.0f%% (\u00b1 %.0f%%)", 100 * sds$district_mean,
              sds$mean)
    )
    plot_districts(simulated, paste("Districts simulated from the fitted",
                                    "model,", insecticides[j]), "#2166AC")
  })
  setNames(plots, insecticides)
}

# save a district figure, at a size that fits its labels
save_district_figure <- function(plot, file) {
  ggplot2::ggsave(file, plot, width = 14, height = 9.5, dpi = 200,
                  bg = "white")
}

# ==== kdr overdispersion by collection class (R/fit_kdr_rho_collection.R) ====

# runs of years at most one year apart: the number of the run each year is in
year_runs <- function(year) {
  years <- sort(unique(year))
  cumsum(c(1, diff(years) > 1))[match(year, years)]
}

# replicate groups of kdr samples (from data/kdr_samples.csv): samples of the
# same pixel and species group, from collection years linked when at most one
# year apart. Only groups of at least two samples are kept
replicate_groups <- function(samples) {
  samples |>
    dplyr::group_by(pixel_id, species_group) |>
    dplyr::mutate(group = paste(pixel_id, species_group,
                                year_runs(collection_year))) |>
    dplyr::group_by(group) |>
    dplyr::filter(dplyr::n() > 1) |>
    dplyr::ungroup()
}

# greta model of kdr overdispersion by collection class, ir_cube's
# hierarchical model of bioassay rho (R/fig_illustrate_bioassay_variability.R
# on its branch rho-study-effect) with kdr samples in place of bioassays and
# collection classes in place of insecticides, with one level of pooling of
# rho. The distribution of the replicate groups' kdr frequencies and the SD
# of the study-year effects are by species group (by insecticide in
# ir_cube). For kdr sample k with A_k kdr alleles of 2 M_k, in replicate
# group r of species group q, from study g in year t, of collection class c:
#   A_k ~ BetaBinomial(2 M_k, p_k, rho_c(k))
#   logit(p_k) = logit(f_r(k)) + s_g,t(k)
#   f_r ~ Beta(m_q(r) kappa_q(r), (1 - m_q(r)) kappa_q(r))
#   s_g,t ~ Normal(0, sigma_s,q(g,t))
#   logit(rho_c) ~ Normal(mu, sigma_rho)
# with ir_cube's priors. rho is in the order of classes
build_rho_collection_model <- function(samples, classes) {
  # ---- levels and indices ----
  species <- sort(unique(samples$species_group))
  groups <- unique(samples$group)
  study_years <- unique(samples$study_year)
  c_k <- match(samples$collection, classes)
  r_k <- match(samples$group, groups)
  gt_k <- match(samples$study_year, study_years)
  q_r <- match(samples$species_group[match(groups, samples$group)], species)
  q_gt <- match(samples$species_group[match(study_years, samples$study_year)],
                species)

  # ---- rho, by collection class ----
  # centred, as each class's rho is well informed by its samples
  mu <- greta::normal(qlogis(0.15), 1)
  sigma_rho <- half_normal(0.5)
  logit_rho <- greta::normal(mu, sigma_rho, dim = length(classes))
  rho <- greta::ilogit(logit_rho)

  # ---- kdr frequency of each replicate group ----
  # mean and concentration of the frequencies in each species group
  m <- greta::beta(1, 1, dim = length(species))
  kappa <- greta::lognormal(log(5), 2, dim = length(species))
  f <- greta::beta(m[q_r] * kappa[q_r], (1 - m[q_r]) * kappa[q_r],
                   dim = length(groups))

  # ---- study-year effect (non-centred), SD by species group ----
  sigma_s <- half_normal(dim = length(species))
  s <- sigma_s[q_gt] * greta::normal(0, 1, dim = length(study_years))

  # ---- kdr samples ----
  p <- greta::ilogit(log(f[r_k]) - log1p(-f[r_k]) + s[gt_k])
  alleles <- 2 * samples$n_tested
  kdr_alleles <- greta::as_data(round(alleles * samples$kdr_frequency))
  greta::distribution(kdr_alleles) <- betabinomial_p_rho(alleles, p,
                                                         rho[c_k])

  greta::model(rho, mu, sigma_rho, sigma_s)
}

# plot rho by collection class (from summarise_matrix(), with collection,
# groups and samples columns): the posterior median (point) and 95% interval
# (whisker) of each class, the first class at the top
plot_rho_collection <- function(estimates) {
  estimates <- estimates |>
    dplyr::mutate(
      label = sprintf("%s\n(%i groups, %i samples)", collection, groups,
                      samples),
      label = factor(label, levels = rev(label))
    )
  ggplot2::ggplot(estimates, ggplot2::aes(y = label)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper),
                            colour = "grey80", linewidth = 1.5) +
    ggplot2::geom_point(ggplot2::aes(x = median), colour = "grey20",
                        size = 2) +
    ggplot2::scale_x_continuous(limits = c(0, NA), n.breaks = 6) +
    ggplot2::labs(x = "kdr intra-class correlation, rho", y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey85")
    )
}

# ==== posterior summaries and checks (both models) ====

# posterior mean, SD and quantiles of each column of a matrix of draws
summarise_matrix <- function(matrix_draws) {
  data.frame(
    parameter = colnames(matrix_draws),
    mean = colMeans(matrix_draws),
    sd = apply(matrix_draws, 2, sd),
    lower = apply(matrix_draws, 2, quantile, 0.025),
    median = apply(matrix_draws, 2, quantile, 0.5),
    upper = apply(matrix_draws, 2, quantile, 0.975),
    row.names = NULL
  )
}

# posterior summaries of each parameter, with the potential scale reduction
# factor (R-hat) and effective sample size, naming the elements of parameters
# with one per insecticide (greta's alpha[1,1]...) by insecticide
summarise_draws <- function(draws, insecticides) {
  matrix_draws <- as.matrix(draws)
  pattern <- "^(.+)\\[(\\d+),1\\]$"
  vector_element <- grepl(pattern, colnames(matrix_draws))
  colnames(matrix_draws)[vector_element] <- paste0(
    sub(pattern, "\\1", colnames(matrix_draws)[vector_element]), "_",
    insecticides[as.integer(sub(pattern, "\\2",
                                colnames(matrix_draws)[vector_element]))]
  )
  summarise_matrix(matrix_draws) |>
    dplyr::mutate(
      r_hat = coda::gelman.diag(draws, multivariate = FALSE,
                                autoburnin = FALSE)$psrf[, 1],
      effective_samples = coda::effectiveSize(draws)
    )
}

# summaries of bioassay mortality, for posterior predictive checks
mortality_summaries <- function(mortality, bioassays) {
  pixel_year_mortality <- tapply(mortality, bioassays$pixel_year, mean)
  c(mean = mean(mortality),
    share_at_100 = mean(mortality >= 1),
    share_below_98 = mean(mortality < 0.98),
    share_below_90 = mean(mortality < 0.90),
    percentile_10 = unname(quantile(mortality, 0.10)),
    pixel_year_sd = sd(pixel_year_mortality))
}

# two-sided p-value of each value of x within x itself: the p-value each
# simulated dataset would get if it were the observed data
two_sided_p_within <- function(x) {
  n <- length(x)
  p_lower <- rank(x, ties.method = "max") / n
  p_upper <- (n - rank(x, ties.method = "min") + 1) / n
  pmin(1, 2 * pmin(p_lower, p_upper))
}

# posterior predictive checks: simulate each bioassay's deaths from posterior
# draws of its expected survival (survival_draws: draws by bioassays), and
# compare summaries of simulated and observed mortality within each group of
# bioassays (by row, a factor such as kdr band or study group, and
# insecticide), with the central 95% of the simulated summaries and a
# two-sided posterior predictive p-value. A check fails if its p-value is below
# min_p_threshold, the Westfall-Young min-p threshold for a 5% family-wise
# error rate across all the checks: the 5th percentile of the smallest
# p-value across the checks of each simulated dataset, which allows for the
# checks being correlated. summaries is a function of the mortality and the
# rows of bioassays in a group, returning a named vector
posterior_predictive_check <- function(survival_draws, bioassays, row,
                                       summaries = mortality_summaries) {
  # beta shapes a and b sum to 1 / rho - 1, as in betabinomial_p_rho()
  phi_i <- 1 / unname(rho[bioassays$insecticide_type]) - 1
  n_i <- bioassays$mosquito_number
  simulated <- apply(survival_draws, 1, function(survival) {
    p_dead <- rbeta(length(survival), (1 - survival) * phi_i,
                    survival * phi_i)
    rbinom(length(survival), n_i, p_dead) / n_i
  })
  observed <- bioassays$died / n_i
  groups <- split(seq_len(nrow(bioassays)),
                  interaction(row, bioassays$insecticide_type, drop = TRUE))
  first <- sapply(groups, `[`, 1)
  # summaries of the simulated data: checks by draws
  simulated_summaries <- lapply(groups, function(rows) {
    apply(simulated[rows, , drop = FALSE], 2, summaries,
          bioassays = bioassays[rows, ])
  })
  n_summaries <- sapply(simulated_summaries, nrow)
  simulated_summaries <- do.call(rbind, simulated_summaries)
  observed_summaries <- unlist(lapply(groups, function(rows) {
    summaries(observed[rows], bioassays[rows, ])
  }))
  stopifnot(length(observed_summaries) == nrow(simulated_summaries))

  p_upper <- rowMeans(simulated_summaries >= observed_summaries)
  p_lower <- rowMeans(simulated_summaries <= observed_summaries)
  p_value <- pmin(1, 2 * pmin(p_upper, p_lower))
  # smallest p-value across the checks of each simulated dataset
  min_p <- apply(apply(simulated_summaries, 1, two_sided_p_within), 1, min)
  min_p_threshold <- unname(quantile(min_p, 0.05))

  data.frame(
    row = rep(row[first], n_summaries),
    insecticide = rep(bioassays$insecticide_type[first], n_summaries),
    summary = rownames(simulated_summaries),
    observed = observed_summaries,
    simulated_median = apply(simulated_summaries, 1, median),
    simulated_lower = apply(simulated_summaries, 1, quantile, 0.025),
    simulated_upper = apply(simulated_summaries, 1, quantile, 0.975),
    p_value = p_value,
    min_p_threshold = min_p_threshold,
    fails = p_value < min_p_threshold,
    row.names = NULL
  )
}

# ---- plots ----

# plot posterior predictive checks (output of posterior_predictive_check()):
# for each summary of bioassay mortality (columns) and row (in the order of
# its factor levels), the observed summary for each insecticide (rows within
# a panel; point)
# against the central 95% of summaries of simulated data (whisker), with the
# two-sided posterior predictive p-value. Checks that fail (p-value below the
# min-p threshold for a 5% family-wise error rate) are in black
plot_checks <- function(checks) {
  # summaries, in the order and with the labels to plot them
  summaries <- c(
    mean = "mean mortality",
    share_at_100 = "share at 100%",
    share_below_98 = "share below 98%",
    share_below_90 = "share below 90%",
    percentile_10 = "10th percentile",
    pixel_year_sd = "SD between pixel-years",
    within_district_sd = "SD within district-years",
    within_district_var = "variance within district-years"
  )
  summaries <- summaries[names(summaries) %in% checks$summary]

  checks <- checks |>
    dplyr::mutate(
      # first insecticide (alphabetically) at the top
      insecticide = factor(insecticide,
                           levels = rev(sort(unique(insecticide)))),
      summary = factor(summaries[summary], levels = summaries),
      p_label = dplyr::case_when(
        p_value < 0.001 ~ "p < 0.001",
        p_value < 0.01 ~ sprintf("p = %.3f", p_value),
        .default = sprintf("p = %.2f", p_value)
      )
    )
  stopifnot(!anyNA(checks$summary))

  ggplot2::ggplot(checks, ggplot2::aes(y = insecticide)) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = simulated_lower, xmax = simulated_upper),
      colour = "grey80", linewidth = 1.5
    ) +
    ggplot2::geom_point(ggplot2::aes(x = observed, colour = fails), size = 2) +
    # p-value at the right edge of the panel, above the insecticide's row
    ggplot2::geom_text(
      ggplot2::aes(x = Inf, label = p_label, colour = fails, size = fails),
      hjust = 1.1, vjust = -1.1, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey45",
                                            `TRUE` = "black"),
                                 guide = "none") +
    ggplot2::scale_size_manual(values = c(`FALSE` = 2.3, `TRUE` = 3),
                               guide = "none") +
    # room above the top insecticide for its p-value
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.4, 0.9))) +
    # rows (kdr bands or study groups) labelled on the left, summaries on
    # top, each summary with its own horizontal scale
    ggplot2::facet_grid(row ~ summary, scales = "free_x", switch = "y") +
    # room on the right so the data stay clear of the p-values
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.35)), n.breaks = 4
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      subtitle = sprintf(
        "Checks fail if p < %.4f (min-p threshold, 5%% family-wise error rate)",
        checks$min_p_threshold[1]
      )
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.spacing.x = ggplot2::unit(1.2, "lines"),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey85"),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1)
    )
}

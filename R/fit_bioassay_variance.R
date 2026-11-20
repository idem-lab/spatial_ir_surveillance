# Fit the bioassay variance model: how much pyrethroid bioassay mortality
# varies between pixels within a district in one year, after the variation
# between districts and over time. For bioassay i of insecticide j in pixel
# l, district d, province r, country c and year t, in which D_i of N_i
# mosquitoes died (build_variance_model()):
#   D_i ~ BetaBinomial(N_i, pi_i, rho_j)
#   logit(1 - pi_i) = W_jlt + B_jdt
# The within-district-year effect W = sigma_v_j v_lt + sigma_e_j e_lt is a
# tapered exponential Gaussian process v over all pixels in a year (one range
# for both insecticides) plus an independent pixel-year effect e (a nugget).
# The between-district-year effect B is a linear trend per decade with
# country departures, and independent effects of country, country-year,
# province, province-year, district and district-year; the province effects
# take the variation shared by neighbouring districts. Each random effect is
# a standardised effect shared between the insecticides, scaled by each
# insecticide's SD; intercepts and trends differ. rho_j, the noise of
# individual bioassays, is fixed at ir_cube's estimates. Fitted to
# district-years with bioassays of that insecticide in several pixels: in
# the others the district-year effect and the pixel-year effects cannot be
# told apart, and their pixels would inform only the correlation between
# districts.
#
# The posterior predictive checks simulate new bioassays with new pixel-year
# effects (v and e, drawn from their fitted SDs and range) at the fitted
# pixels, with the between-district-year effects at their posterior draws, so
# the within-district-year summaries test whether the fitted parameters
# reproduce the spread between the pixels of a district-year. The other
# summaries span district-years, which the new v also shifts, as the GP spans
# district borders. The checks are split by district-years sampled by one
# study and by several, since differences between studies would add to the
# spread there. The figures show districts simulated from the priors and from
# the posterior (the archetypal district: a circle 100 km across, 5 km
# pixels, with mean mortality 50, 80 or 95%).
#
# Inputs: data/bioassays.RDS (R/prep_bioassays.R),
#   data/ir_cube/raster_mask.tif (R/get_data.R) and
#   metadata/bioassay_exclusions.csv
# Outputs: output/bioassay_variance_estimates.csv (posterior summaries);
#   output/bioassay_variance_fit.RDS (the posterior draws, with what
#   fresh_survival() needs); output/bioassay_variance_checks.csv and
#   figures/bioassay_variance_checks.png (posterior predictive checks);
#   output/bioassay_variance_district_sd.csv (the SD of mortality between
#   pixels of simulated districts); figures/prior_districts.png,
#   figures/prior_districts_no_pixel_effect.png and
#   figures/posterior_districts_<insecticide>.png (simulated districts)

source("R/packages.R")
source("R/functions.R")
set.seed(2026)

# MCMC settings: longer trajectories and warmup than greta's defaults (5-10
# leapfrog steps, 1000 warmup), which the SDs of the nested space-time
# effects need to mix
sampler <- greta::hmc(Lmin = 20, Lmax = 30)
n_samples <- 1000
warmup <- 2000
chains <- 4

# distance (km) beyond which pixel-year effects are independent. The taper
# makes the correlation matrix block-diagonal over clusters of nearby
# pixels, so each cluster gets its own Cholesky factor, saving time and
# memory; Inf gives the untapered model
taper_km <- 400

# simulated districts: their mean mortality, the number of prior and
# posterior draws to plot, and the number of posterior draws for the SD of
# mortality between their pixels
district_means <- c(0.5, 0.8, 0.95)
n_district_draws <- 6
n_sd_draws <- 1000

# ---- data ----

# bioassays excluded after checking against their sources (decisions and
# reasons, by rule, in metadata/bioassay_exclusions.csv). High spread: in
# district-years where the SD of the empirical logit of mortality between
# bioassays exceeded 3, records shown to be database errors, or from sources
# that cannot be checked. Duplicates: in district-years whose bioassays agree
# more closely than replicate tests of one population would 95% of the time
# (beta-binomial at rho, pooled mortality; district-years with all bioassays
# at 0% or all at 100% excluded from this screen), pairs of records from
# different sources with mortality within 0.5 percentage points; where two
# sources matched at least twice and in at least half of their pairs of
# records, the copy (from PMI or an institute database) of the data
# collector's result was dropped
excluded <- read.csv("metadata/bioassay_exclusions.csv") |>
  dplyr::filter(decision == "exclude")

selected <- select_bioassays(readRDS("data/bioassays.RDS"))
stopifnot(all(excluded$record_id %in% selected$record_id))

# only district-years with bioassays of that insecticide in several pixels
bioassays <- selected |>
  dplyr::filter(!record_id %in% excluded$record_id) |>
  dplyr::group_by(insecticide_type, district_year) |>
  dplyr::filter(dplyr::n_distinct(pixel_year) > 1) |>
  dplyr::ungroup()

# ---- fit ----

pixel_xy <- terra::xyFromCell(terra::rast("data/ir_cube/raster_mask.tif"),
                              bioassays$pixel_id)
built <- build_variance_model(bioassays, pixel_xy, taper_km)
draws <- greta::mcmc(built$model, sampler = sampler, n_samples = n_samples,
                     warmup = warmup, chains = chains, verbose = FALSE)

write.csv(summarise_draws(draws, built$insecticides),
          "output/bioassay_variance_estimates.csv", row.names = FALSE)
saveRDS(c(list(draws = draws), built[names(built) != "model"]),
        "output/bioassay_variance_fit.RDS")

# ---- checks ----

# within each insecticide, for district-years sampled by one study and by
# several, from all the posterior draws with new pixel-year effects
studies <- bioassays |>
  dplyr::group_by(insecticide_type, district_year) |>
  dplyr::mutate(studies = dplyr::if_else(dplyr::n_distinct(study_year) > 1,
                                         "several studies", "one study")) |>
  dplyr::pull(studies) |>
  factor(levels = c("one study", "several studies"))
survival_draws <- greta::calculate(
  fresh_survival(built), values = draws,
  nsim = coda::niter(draws) * coda::nchain(draws)
)[[1]][, , 1]
checks <- posterior_predictive_check(survival_draws, bioassays, studies,
                                     summaries = variance_summaries)
write.csv(checks, "output/bioassay_variance_checks.csv", row.names = FALSE)
ggplot2::ggsave("figures/bioassay_variance_checks.png",
                plot_checks(checks),
                width = 16, height = 4.5, dpi = 200, bg = "white")

# ---- simulated districts ----

prior_plots <- prior_district_plots(district_means, taper_km,
                                    n_district_draws)
save_district_figure(prior_plots$with_e, "figures/prior_districts.png")
save_district_figure(prior_plots$without_e,
                     "figures/prior_districts_no_pixel_effect.png")

posterior <- as.matrix(draws)
district_sds <- district_sd_summary(posterior, built$insecticides,
                                    district_means, taper_km, n_sd_draws)
write.csv(district_sds, "output/bioassay_variance_district_sd.csv",
          row.names = FALSE)
posterior_plots <- posterior_district_plots(posterior, built$insecticides,
                                            district_means, district_sds,
                                            taper_km, n_district_draws)
for (insecticide in names(posterior_plots)) {
  save_district_figure(posterior_plots[[insecticide]],
                       sprintf("figures/posterior_districts_%s.png",
                               insecticide))
}

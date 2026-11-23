# Fit the kdr predictor of pyrethroid bioassay mortality: survival without kdr
# (alpha), the effect of each kdr allele (beta), both for each insecticide,
# and how much each varies between populations (sigma_u, sigma_b). The model
# (build_kdr_model()) has three likelihoods:
#   - bioassays, whose logit survival is linear in the population's mean
#     number of kdr alleles, 2f, where f is its kdr frequency
#   - genotypes of bioassay survivors and dead (IR Mapper's matched pairs),
#     which inform beta within populations
#   - kdr samples from the same pixel and year as the bioassays, which
#     inform f
# There is no time trend: changes over time in survival without kdr are part
# of the variation between pixel-years.
#
# Inputs: data/irmapper/vgsc_survivors_dead.csv (R/get_data.R),
#   data/bioassays.RDS (R/prep_bioassays.R), data/kdr_samples.csv
#   (R/prep_kdr_samples.R)
# Outputs: output/kdr_predictor_estimates.csv (posterior summaries),
#   output/kdr_predictor_checks.csv and figures/kdr_predictor_checks.png
#   (posterior predictive checks),
#   output/kdr_predictor_fit.RDS (the posterior draws, with the greta arrays
#   of expected survival and the pixel-years' u, b and f, for
#   greta::calculate()), output/kdr_predictor_bioassays.RDS (the bioassays)

source("R/packages.R")
source("R/functions.R")
set.seed(2026)

# MCMC settings: twice greta's default number of leapfrog steps per HMC
# iteration (5 to 10)
sampler <- greta::hmc(Lmin = 10, Lmax = 20)
n_samples <- 5000
warmup <- 5000
chains <- 4
# number of posterior draws used for posterior predictive checks
n_check_draws <- 10000

# the kdr samples of the study whose bioassays select_bioassays() leaves out
# (Yadouleton et al. 2010), by their paper key from R/prep_kdr_samples.R
excluded_kdr_paper <- "10.1186/14752875983"

# ---- data ----

bioassays <- select_bioassays(readRDS("data/bioassays.RDS"))
kdr_samples <- read.csv("data/kdr_samples.csv") |>
  dplyr::filter(paper != excluded_kdr_paper)

# kdr samples from the same pixel and year as bioassays: one row per
# pixel-year and sample
links <- find_same_pixel_samples(bioassays, kdr_samples)

selected <- select_kdr_sampled(bioassays, links) |>
  dplyr::select(record_id, selection, pixel_year, insecticide_type,
                mosquito_number, died)
samples <- dplyr::filter(links, pixel_year %in% selected$pixel_year)
pairs <- read_pair_genotypes("data/irmapper/vgsc_survivors_dead.csv")

# ---- fit ----

built <- build_kdr_model(selected, samples, pairs)
draws <- greta::mcmc(built$model, sampler = sampler, n_samples = n_samples,
                     warmup = warmup, chains = chains, verbose = FALSE)

write.csv(dplyr::bind_rows(summarise_draws(draws, built$insecticides),
                           kdr_mortality(draws, built$insecticides)),
          "output/kdr_predictor_estimates.csv", row.names = FALSE)
saveRDS(c(list(draws = draws, survival = built$survival,
               pixel_years = built$pixel_years),
          built$latents),
        "output/kdr_predictor_fit.RDS")
saveRDS(selected, "output/kdr_predictor_bioassays.RDS")

# ---- checks ----

# within each band of kdr frequency and insecticide, from n_check_draws
# posterior draws of each bioassay's expected survival
survival_draws <- as.matrix(greta::calculate(built$survival, values = draws))
survival_draws <- survival_draws[sample.int(nrow(survival_draws),
                                            n_check_draws), ]
checks <- posterior_predictive_check(survival_draws, selected,
                                     selected$selection)
write.csv(checks, "output/kdr_predictor_checks.csv", row.names = FALSE)
ggplot2::ggsave("figures/kdr_predictor_checks.png", plot_checks(checks),
                width = 12, height = 6, dpi = 200, bg = "white")

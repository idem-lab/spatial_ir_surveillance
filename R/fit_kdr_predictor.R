# Fit the kdr predictor of pyrethroid bioassay mortality, a hazard model with
# co-selection. For a bioassay of insecticide j, the log hazard of dying is
#   log H = alpha* - 2 f exp(lambda_j + b) - 2 gamma_j ilogit(eta)
# with survival exp(-H), kdr frequency f, and eta the logit frequency of the
# other drivers of resistance, with correlation psi between eta and logit f
# (co-selection). alpha* is fixed from the WHO discriminating concentration
# (twice the LC99.9). The model (build_kdr_model()) has three likelihoods:
#   - bioassays, through the log hazard above
#   - genotypes of bioassay survivors and dead (IR Mapper's matched pairs),
#     which inform the kdr effect within populations
#   - kdr samples from the same pixel and year as the bioassays, which
#     inform f
# There is no time trend: changes over time in the other drivers are part of
# the variation between pixel-years.
#
# Inputs: data/irmapper/vgsc_survivors_dead.csv (R/get_data.R),
#   data/bioassays.RDS (R/prep_bioassays.R), data/kdr_samples.csv
#   (R/prep_kdr_samples.R)
# Outputs: output/kdr_predictor_estimates.csv (posterior summaries),
#   output/kdr_predictor_checks.csv and figures/kdr_predictor_checks.png
#   (posterior predictive checks),
#   output/kdr_predictor_fit.RDS (the posterior draws, with the greta arrays
#   of expected survival and the pixel-years' latent variables, for
#   greta::calculate()),
#   output/kdr_predictor_bioassays.RDS (the bioassays)

source("R/packages.R")
source("R/functions.R")
set.seed(2026)

# MCMC settings: 20 to 40 leapfrog steps (greta's default is 5 to 10), where
# effective draws per second level off; 4000 draws (greta's default is 2000)
# to reach the effective sample size target, rather than longer trajectories
sampler <- greta::hmc(Lmin = 20, Lmax = 40)
n_samples <- 4000
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
                     chains = chains, verbose = FALSE)

write.csv(dplyr::bind_rows(summarise_draws(draws, built$insecticides),
                           kdr_mortality(draws, built$insecticides)),
          "output/kdr_predictor_estimates.csv", row.names = FALSE)
saveRDS(c(list(draws = draws), built[names(built) != "model"]),
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

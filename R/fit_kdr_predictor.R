# Fit the kdr predictor of pyrethroid bioassay mortality, a hazard model with
# co-selection. For a bioassay of insecticide j, the log hazard of dying is
#   log H = alpha* - 2 f exp(lambda_j + b) - 2 gamma_j ilogit(eta)
# with survival exp(-H), kdr frequency f, and eta the logit frequency of the
# other drivers of resistance, correlated psi with logit f (co-selection).
# alpha* is fixed from the WHO discriminating concentration (twice the
# LC99.9). The model (build_kdr_model()) has three likelihoods:
#   - bioassays, whose log hazard falls with the population's mean number of
#     kdr alleles, 2f, where f is its kdr frequency, and with the other
#     drivers
#   - genotypes of bioassay survivors and dead (IR Mapper's matched pairs),
#     which inform the kdr effect within populations
#   - kdr samples from the same pixel and year as the bioassays, which
#     inform f
# There is no time trend: changes over time in the other drivers are part of
# the variation between pixel-years.
#
# The posterior predictive checks use each bioassay's fitted survival. The
# held-out checks leave out the bioassays (not the kdr samples) of each of
# 10 folds of pixel-years and predict them with new b and eta given the
# fitted logit f, giving the elpd and the PIT of pooled mortality.
#
# Inputs: data/irmapper/vgsc_survivors_dead.csv (R/get_data.R),
#   data/bioassays.RDS (R/prep_bioassays.R), data/kdr_samples.csv
#   (R/prep_kdr_samples.R)
# Outputs: output/kdr_predictor_estimates.csv (posterior summaries),
#   output/kdr_predictor_checks.csv and figures/kdr_predictor_checks.png
#   (posterior predictive checks), output/kdr_predictor_heldout_checks.csv
#   and figures/kdr_predictor_heldout_checks.png (PIT checks of held-out
#   bioassays), output/kdr_predictor_cv.csv (held-out elpd),
#   output/kdr_predictor_fit.RDS (the posterior draws, with the greta arrays
#   of expected survival, the pixel-years' latent variables and what
#   new_mortality() needs, for greta::calculate()),
#   output/kdr_predictor_bioassays.RDS (the bioassays)

source("R/packages.R")
source("R/functions.R")
set.seed(2026)

# MCMC settings: 20 to 40 leapfrog steps (greta's default is 5 to 10): the
# many pixel-year random effects force small steps, and sampler efficiency
# (effective draws per second) levels off from about 20 steps, so the final
# run uses more draws rather than longer trajectories
sampler <- greta::hmc(Lmin = 20, Lmax = 40)
n_samples <- 4000
chains <- 4
# sets of uniform PITs for the held-out calibration checks
n_uniform_draws <- 4000
# cross-validation folds, and the draws of the new b and eta for each
# posterior draw when integrating them out of the held-out likelihood
n_folds <- 10
n_fresh <- 5
# the fold fits use longer trajectories: with 20 to 40 steps, a chain of one
# fold stayed in a different region from the others (R-hat 1.5)
fold_sampler <- greta::hmc(Lmin = 50, Lmax = 100)

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

# within each band of kdr frequency and insecticide, from the posterior draws
# of each bioassay's expected survival
survival_draws <- as.matrix(greta::calculate(built$survival, values = draws))
checks <- posterior_predictive_check(survival_draws, selected,
                                     selected$selection)
write.csv(checks, "output/kdr_predictor_checks.csv", row.names = FALSE)
ggplot2::ggsave("figures/kdr_predictor_checks.png", plot_checks(checks),
                width = 12, height = 6, dpi = 200, bg = "white")

# ---- held-out checks ----

pixel_years <- unique(selected$pixel_year)
fold <- sample(rep(seq_len(n_folds), length.out = length(pixel_years)))
fold_i <- fold[match(selected$pixel_year, pixel_years)]
heldout <- lapply(seq_len(n_folds), function(k) {
  fold_built <- build_kdr_model(selected[fold_i != k, ], samples, pairs)
  fold_draws <- greta::mcmc(fold_built$model, sampler = fold_sampler,
                            n_samples = n_samples, chains = chains,
                            verbose = FALSE)
  c(heldout_scores(fold_built, fold_draws, selected[fold_i == k, ], n_fresh),
    max_r_hat = max(coda::gelman.diag(fold_draws, multivariate = FALSE,
                                      autoburnin = FALSE)$psrf[, 1]))
})

elpd <- loo::elpd(do.call(cbind, lapply(heldout, `[[`, "loglik")))
write.csv(data.frame(elpd = elpd$estimates["elpd", "Estimate"],
                     se = elpd$estimates["elpd", "SE"],
                     pixel_years = length(pixel_years),
                     max_fold_r_hat = max(sapply(heldout, `[[`,
                                                 "max_r_hat"))),
          "output/kdr_predictor_cv.csv", row.names = FALSE)

# PIT of the held-out pooled mortality in each band of kdr frequency and
# insecticide
heldout_checks <- pit_check(do.call(rbind, lapply(heldout, `[[`, "pit")),
                            n_uniform_draws)
write.csv(heldout_checks, "output/kdr_predictor_heldout_checks.csv",
          row.names = FALSE)
ggplot2::ggsave("figures/kdr_predictor_heldout_checks.png",
                plot_checks(heldout_checks),
                width = 6, height = 6, dpi = 200, bg = "white")

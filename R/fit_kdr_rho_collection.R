# Fit the overdispersion (intra-class correlation, rho) of kdr allele counts
# between replicate genotyped samples, separately by how the mosquitoes were
# collected (the collection class from R/prep_kdr_samples.R). Larval
# collections pool few breeding sites and contain many siblings, while
# host-seeking adults should be better mixed by flight, so rho may differ
# between them.
#
# The model (build_rho_collection_model()) is ir_cube's hierarchical model of
# bioassay rho: kdr alleles in each sample are beta-binomial around their
# replicate group's frequency, with a study-year effect, so that rho is the
# noise between samples of one study. rho is fitted for each collection
# class, pooled across classes. Replicate groups are samples of the same
# pixel and species group, from collection years linked when at most one
# year apart (pixels, species groups, collection years and study-years are
# defined in R/prep_kdr_samples.R). Samples pooled from several collection
# methods, or with an unknown method, are left out.
#
# Input: data/kdr_samples.csv (R/prep_kdr_samples.R)
# Outputs: output/kdr_rho_collection_estimates.csv (posterior summaries of
#   rho by collection class, with the probability that each class's rho
#   exceeds that of host-seeking adults), output/kdr_rho_collection_fit.RDS
#   (the joint posterior draws of rho, a matrix with a column for each
#   class) and figures/kdr_rho_collection.png

source("R/packages.R")
source("R/functions.R")
set.seed(2026)

# MCMC settings: greta's defaults, but with three times its leapfrog steps
# per HMC iteration (5 to 10), as ir_cube uses for this model. With the
# default steps, rho and the SDs of the study-year effects mixed poorly
# (R-hat up to 1.2), with rho centred (well informed) and the study-year
# effects non-centred (most replicate groups have samples of one study)
sampler <- greta::hmc(Lmin = 15, Lmax = 30)

# the collection classes fitted
classes <- c("larvae", "house resting or exit trap",
             "resting outdoors", "human landing or light trap")

# ---- data ----

samples <- read.csv("data/kdr_samples.csv") |>
  dplyr::filter(collection %in% classes) |>
  replicate_groups()

# ---- fit ----

model <- build_rho_collection_model(samples, classes)
draws <- greta::mcmc(model, sampler = sampler, verbose = FALSE)

# ---- estimates ----

rho_draws <- as.matrix(draws)[, sprintf("rho[%i,1]", seq_along(classes))]
colnames(rho_draws) <- classes
counts <- dplyr::summarise(samples, groups = dplyr::n_distinct(group),
                           samples = dplyr::n(), .by = collection)
estimates <- summarise_matrix(rho_draws) |>
  dplyr::rename(collection = parameter) |>
  dplyr::left_join(counts, by = dplyr::join_by(collection)) |>
  dplyr::relocate(groups, samples, .after = collection) |>
  dplyr::mutate(p_above_host_seeking = dplyr::if_else(
    collection == "human landing or light trap", NA,
    colMeans(rho_draws > rho_draws[, "human landing or light trap"])
  ))
print(estimates, digits = 3)

write.csv(estimates, "output/kdr_rho_collection_estimates.csv",
          row.names = FALSE)
saveRDS(rho_draws, "output/kdr_rho_collection_fit.RDS")
ggplot2::ggsave("figures/kdr_rho_collection.png",
                plot_rho_collection(estimates),
                width = 6, height = 3, dpi = 200, bg = "white")

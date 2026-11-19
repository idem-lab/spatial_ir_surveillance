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

# ==== kdr predictor (R/fit_kdr_predictor.R) ====

# bioassay intra-class correlation of mortality by insecticide: the posterior
# means in ir_cube's outputs/bioassay_rho_hierarchical.csv (untracked; from its
# final run of 30 Oct 2026 on the deduplicated bioassays our data come from),
# made by R/fig_illustrate_bioassay_variability.R, which fits rho
# hierarchically to all replicate bioassays in a pixel, year and insecticide.
# ir_cube's two-stage model uses the same values
rho <- c("Deltamethrin" = 0.1492, "Permethrin" = 0.1401)

# beta-binomial distribution in greta with mean p and intra-class correlation
# rho (ir_cube's, from its R/functions.R)
betabinomial_p_rho <- function(N, p, rho) {
  a <- p * (1 / rho - 1)
  b <- a * (1 - p) / p
  greta::beta_binomial(size = N, alpha = a, beta = b)
}

# pyrethroids as spelled in the IR Mapper survivor/dead file, and as in the
# bioassay data
pair_pyrethroids <- c(deltamethrin = "Deltamethrin",
                      permethrin = "Permethrin")

# ---- data ----

# kdr samples in the same pixel as a bioassay pixel-year (pixels of mask),
# collected in that year, one row per (pixel-year, sample) pair. A sample
# collected over several years that match bioassays in the same pixel is
# assigned to the one nearest the middle of its collection years (the earlier
# on a tie), so each sample is used once
find_same_pixel_samples <- function(bioassays, kdr_samples, mask) {
  pixel_years <- dplyr::distinct(bioassays, pixel_year, pixel_id, year)
  kdr_samples |>
    dplyr::mutate(pixel_id = terra::cellFromXY(mask,
                                               cbind(longitude, latitude))) |>
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
                                  right = FALSE) |>
                    as.character())
  dplyr::inner_join(bioassays, dplyr::select(sampled, pixel_year, selection),
                    by = dplyr::join_by(pixel_year))
}

# genotype counts of survivors and dead in IR Mapper's survivor/dead pairs
# tested with deltamethrin or permethrin, other than An. arabiensis: one row
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

# ---- model ----

# greta model, from the selected bioassays, their kdr samples (from
# find_same_pixel_samples()) and the survivor/dead pairs (from
# read_pair_genotypes()). For bioassay i of insecticide j in pixel-year
# (l, t), with kdr frequency f:
#   logit(1 - pi_i) = alpha_j + u_lt + 2 f_lt (beta_j + b_lt)
# Returns the model; the greta array of each bioassay's expected survival,
# 1 - pi_i (for posterior predictive checks); the insecticides, in the order
# of alpha and beta; and the greta arrays of each pixel-year's latent
# variables, u, b and f, with the pixel-years in their order
build_model <- function(bioassays, samples, pairs) {
  # ---- indices ----
  pixel_years <- unique(bioassays$pixel_year)
  lt_i <- match(bioassays$pixel_year, pixel_years)
  insecticides <- sort(unique(bioassays$insecticide_type))
  j_i <- match(bioassays$insecticide_type, insecticides)
  # the pairs inform the kdr effects of the bioassays' insecticides
  stopifnot(all(pairs$insecticide %in% insecticides))
  rho_i <- unname(rho[bioassays$insecticide_type])
  stopifnot(!anyNA(rho_i))
  n_lt <- length(pixel_years)

  # ---- priors ----
  alpha <- greta::normal(0, 3, dim = length(insecticides))
  # kdr effect for each insecticide
  beta <- greta::normal(0, 1, dim = length(insecticides))
  sigma_u <- greta::normal(0, 1, truncation = c(0, Inf))
  sigma_b <- greta::normal(0, 1, truncation = c(0, Inf))

  # ---- kdr frequency ----
  f <- greta::beta(1, 1, dim = n_lt)
  # kdr alleles A_k in each sample k of M_k mosquitoes: beta-binomial around
  # its pixel-year's frequency, with intra-class correlation rho_kdr (uniform
  # prior)
  rho_kdr <- greta::beta(1, 1)
  lt_k <- match(samples$pixel_year, pixel_years)
  stopifnot(!anyNA(lt_k))
  f_k <- f[lt_k]
  alleles <- 2 * samples$n_tested
  kdr_alleles <- greta::as_data(round(alleles * samples$kdr_frequency))
  greta::distribution(kdr_alleles) <- betabinomial_p_rho(alleles, f_k,
                                                         rho_kdr)

  # ---- bioassays ----
  # logit survival is linear in the mean number of kdr alleles, 2f
  u <- sigma_u * greta::normal(0, 1, dim = n_lt)
  b <- sigma_b * greta::normal(0, 1, dim = n_lt)
  survival <- greta::ilogit(
    alpha[j_i] + u[lt_i] + 2 * f[lt_i] * (beta[j_i] + b[lt_i])
  )
  died <- greta::as_data(bioassays$died)
  greta::distribution(died) <- betabinomial_p_rho(bioassays$mosquito_number,
                                                  1 - survival, rho_i)

  # ---- survivor/dead pairs ----
  # genotyped survivors L of L + D with a kdr alleles in pair h: binomial with
  # logit(theta_ha) = gamma_h + (beta_j + b'_h) a, where gamma_h is a nuisance
  # intercept and b'_h the pair's deviation in the kdr effect
  pair_ids <- unique(pairs$pair)
  h <- match(pairs$pair, pair_ids)
  j_h <- match(pairs$insecticide, insecticides)
  gamma <- greta::normal(0, 3, dim = length(pair_ids))
  b_pair <- sigma_b * greta::normal(0, 1, dim = length(pair_ids))
  # probability that a genotyped mosquito is a survivor
  theta <- greta::ilogit(gamma[h] + (beta[j_h] + b_pair[h]) * pairs$a)
  survivors <- greta::as_data(pairs$survivors)
  greta::distribution(survivors) <- greta::binomial(
    pairs$survivors + pairs$dead, theta
  )

  list(model = greta::model(alpha, beta, sigma_u, sigma_b, rho_kdr),
       survival = survival, insecticides = insecticides,
       latents = list(u = u, b = b, f = f),
       pixel_years = pixel_years)
}

# ---- summaries ----

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
# factor (R-hat) and effective sample size, naming alpha and beta's elements
# (greta's alpha[1,1]...) by insecticide. Followed by the same parameters as
# bioassay mortality for each insecticide, for a population with average
# deviations (u = b = 0): without kdr, 1 - ilogit(alpha), and with kdr fixed,
# 1 - ilogit(alpha + 2 beta)
summarise_draws <- function(draws, insecticides) {
  matrix_draws <- as.matrix(draws)
  pattern <- "^(alpha|beta)\\[(\\d+),1\\]$"
  vector_element <- grepl(pattern, colnames(matrix_draws))
  colnames(matrix_draws)[vector_element] <- paste0(
    sub(pattern, "\\1", colnames(matrix_draws)[vector_element]), "_",
    insecticides[as.integer(sub(pattern, "\\2",
                                colnames(matrix_draws)[vector_element]))]
  )
  estimates <- summarise_matrix(matrix_draws) |>
    dplyr::mutate(
      r_hat = coda::gelman.diag(draws, multivariate = FALSE,
                                autoburnin = FALSE)$psrf[, 1],
      effective_samples = coda::effectiveSize(draws)
    )

  alpha <- matrix_draws[, paste0("alpha_", insecticides), drop = FALSE]
  beta <- matrix_draws[, paste0("beta_", insecticides), drop = FALSE]
  mortality <- cbind(1 - plogis(alpha), 1 - plogis(alpha + 2 * beta))
  colnames(mortality) <- c(paste0("mortality_kdr_absent_", insecticides),
                           paste0("mortality_kdr_fixed_", insecticides))
  dplyr::bind_rows(estimates, summarise_matrix(mortality))
}

# summaries of bioassay mortality, for posterior predictive checks
mortality_summaries <- function(mortality, pixel_year) {
  pixel_year_mortality <- tapply(mortality, pixel_year, mean)
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
# bioassays, with the central 95% of the simulated summaries and a two-sided
# posterior predictive p-value. A check fails if its p-value is below
# min_p_threshold, the Westfall-Young min-p threshold for a 5% family-wise
# error rate across all the checks: the 5th percentile of the smallest
# p-value across the checks of each simulated dataset, which allows for the
# checks being correlated
posterior_predictive_check <- function(survival_draws, bioassays, group) {
  # beta shapes a and b sum to 1 / rho - 1, as in betabinomial_p_rho()
  phi_i <- 1 / unname(rho[bioassays$insecticide_type]) - 1
  n_i <- bioassays$mosquito_number
  simulated <- apply(survival_draws, 1, function(survival) {
    p_dead <- rbeta(length(survival), (1 - survival) * phi_i,
                    survival * phi_i)
    rbinom(length(survival), n_i, p_dead) / n_i
  })
  observed <- bioassays$died / n_i
  groups <- split(seq_len(nrow(bioassays)), group)
  # summaries of the simulated data: checks by draws
  simulated_summaries <- lapply(groups, function(rows) {
    apply(simulated[rows, , drop = FALSE], 2, mortality_summaries,
          pixel_year = bioassays$pixel_year[rows])
  })
  simulated_summaries <- do.call(rbind, simulated_summaries)
  observed_summaries <- unlist(lapply(groups, function(rows) {
    mortality_summaries(observed[rows], bioassays$pixel_year[rows])
  }))
  stopifnot(length(observed_summaries) == nrow(simulated_summaries))

  p_upper <- rowMeans(simulated_summaries >= observed_summaries)
  p_lower <- rowMeans(simulated_summaries <= observed_summaries)
  p_value <- pmin(1, 2 * pmin(p_upper, p_lower))
  # smallest p-value across the checks of each simulated dataset
  min_p <- apply(apply(simulated_summaries, 1, two_sided_p_within), 1, min)
  min_p_threshold <- unname(quantile(min_p, 0.05))

  data.frame(
    group = rep(names(groups), each = nrow(simulated_summaries) /
                  length(groups)),
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

# plot posterior predictive checks (output of posterior_predictive_check(),
# with groups named "band: insecticide"): for each summary of bioassay
# mortality (columns) and band of kdr frequency (rows), the observed summary
# for each insecticide (rows within a panel; point) against the central 95% of
# summaries of simulated data (whisker), with the two-sided posterior
# predictive p-value. Checks that fail (p-value below the min-p threshold for
# a 5% family-wise error rate) are in black
plot_checks <- function(checks) {
  # summaries, in the order and with the labels to plot them
  summaries <- c(
    mean = "mean mortality",
    share_at_100 = "share at 100%",
    share_below_98 = "share below 98%",
    share_below_90 = "share below 90%",
    percentile_10 = "10th percentile",
    pixel_year_sd = "SD between pixel-years"
  )
  bands <- c("kdr below 0.3", "kdr 0.3 to 0.7", "kdr 0.7 or above")

  checks <- checks |>
    dplyr::mutate(
      band = factor(sub(": .*", "", group), levels = bands),
      # first insecticide (alphabetically) at the top
      insecticide = factor(sub(".*: ", "", group),
                           levels = rev(sort(unique(sub(".*: ", "", group))))),
      summary = factor(summaries[summary], levels = summaries),
      p_label = dplyr::case_when(
        p_value < 0.001 ~ "p < 0.001",
        p_value < 0.01 ~ sprintf("p = %.3f", p_value),
        .default = sprintf("p = %.2f", p_value)
      )
    )
  stopifnot(!anyNA(checks$band), !anyNA(checks$summary))

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
    # kdr bands labelled on the left, summaries on top, each summary with its
    # own horizontal scale
    ggplot2::facet_grid(band ~ summary, scales = "free_x", switch = "y") +
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

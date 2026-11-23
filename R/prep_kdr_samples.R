# kdr allele frequencies in genotyped samples of the An. gambiae complex, from
# the Vector Atlas and the IR Mapper analysis-ready data.
#
# kdr here is either of the two Vgsc mutations at codon 995 (1014 in the older
# numbering): 995F (L1014F) or 995S (L1014S). For each sample we compute the
# frequency of kdr alleles, from genotype counts where they are given and from
# reported allele percentages otherwise.
#
# We use only sources that give genotype counts or allele percentages in fixed
# fields. IR Mapper's other exports give frequencies as free text and are left
# out. The Vector Atlas and IR Mapper overlap, and so do IR Mapper's species
# and complex files, so some samples appear more than once: these copies are
# removed (see "copies of the same sample", below).
#
# Each sample also has its collection class (collection_class(), below: how
# its mosquitoes were collected), its pixel, species group, collection year
# and study-year (see "pixel, species group and year", below). Samples of
# An. arabiensis are left out.
#
# Inputs (R/get_data.R): data/vector_atlas/va_full_data.csv,
#   data/irmapper/vgsc_species.csv, data/irmapper/vgsc_complex.csv,
#   data/ir_cube/raster_mask.tif
# Output: data/kdr_samples.csv

source("R/packages.R")
source("R/functions.R")

# the number of mosquitoes closest to typical, within range, whose 2n
# alleles give a whole number of kdr alleles at frequency f (typical if none
# does), as ir_cube infers bioassays' missing sample sizes
infer_sample_size <- function(f, typical, range) {
  candidates <- seq(range[1], range[2])
  candidates <- candidates[order(abs(candidates - typical))]
  vapply(f, function(f) {
    alleles <- 2 * candidates * f
    whole <- abs(alleles - round(alleles)) < 1e-9
    if (isTRUE(any(whole))) candidates[which(whole)[1]] else typical
  }, numeric(1))
}

# species labels in the An. gambiae complex, as they appear in the two sources
gambiae_complex <- "gambiae|coluzzii|arabiensis|melas|merus|quadriannulatus"

# number genotyped and kdr allele frequency from a source's genotype count
# columns (in the order of genotypes), both NA if no counts are given (or
# they are all 0, as when only one genotype's count is filled in, as 0)
from_genotypes <- function(counts) {
  counts <- as.matrix(counts)
  counts[is.na(counts)] <- 0
  n <- rowSums(counts)
  n[n == 0] <- NA
  data.frame(n_genotyped = n,
             from_counts = drop(counts %*% genotypes$kdr_alleles) / (2 * n))
}

# the collection class of each collection method, as the two sources record
# them (lower case): "unknown" for "?", NA if not in these classes
method_class <- function(method) {
  dplyr::case_when(
    method %in% "?" ~ "unknown",
    grepl("larva", method) ~ "larvae",
    grepl("^hri$|house resting|exit trap|winexit|veranda", method) ~
      "house resting or exit trap",
    grepl("^ro\\b", method) ~ "resting outdoors",
    # host-seeking: human landing catches, human-baited nets, light traps and
    # tent traps (IR Mapper's, the human-baited Ifakara tent trap)
    grepl("^hlc|^hbn|man biting|human baited|light trap|cdc-lt|tent trap",
          method) ~ "human landing or light trap"
  )
}

# how the mosquitoes of each sample were collected, from a source's
# collection method columns (one row per sample): the class of its methods,
# "pooled" if they fall in more than one class (mosquitoes collected in
# different ways, or perhaps so where one method is "?"), or "unknown" if no
# method is recorded or it is "?". A method in none of the classes counts as
# a class of its own, NA if it is the sample's only method
collection_class <- function(methods) {
  methods <- tolower(trimws(as.matrix(methods)))
  methods[methods %in% ""] <- NA
  classes <- matrix(method_class(methods), nrow = nrow(methods))
  vapply(seq_len(nrow(methods)), function(i) {
    found <- unique(classes[i, !is.na(methods[i, ])])
    if (length(found) > 1) "pooled"
    else if (length(found) == 1) found
    else "unknown"
  }, character(1))
}

# ---- Vector Atlas ----

va <- read.csv("data/vector_atlas/va_full_data.csv")
va_samples <- va |>
  dplyr::mutate(record_id = paste0("va_", dplyr::row_number())) |>
  dplyr::filter(
    # records flagged confidential are not used
    !confidentiality_status %in% "yes",
    grepl(gambiae_complex, species),
    # genotypes of only the survivors or only the dead of a bioassay don't
    # estimate the population's kdr frequency
    !tolower(bioassay_subsample_used_in_genotypic_test) %in%
      c("alive", "alive only", "dead", "dead only")
  ) |>
  dplyr::transmute(
    source = "vector_atlas",
    record_id,
    citation = citation_doi,
    month = month_start,
    latitude = latitude_1,
    longitude = longitude_1,
    year_start,
    year_end,
    species,
    from_genotypes(dplyr::pick(dplyr::all_of(genotypes$vector_atlas))),
    # allele percentage: 995F plus 995S where either is given, otherwise the
    # combined kdr percentage
    percent = dplyr::if_else(
      !is.na(vgsc995f_percent) | !is.na(vgsc995s_percent),
      dplyr::coalesce(vgsc995f_percent, 0) +
        dplyr::coalesce(vgsc995s_percent, 0),
      kdr_percent),
    n_tested = dplyr::coalesce(n_genotyped, vgsc_number_of_mosquitoes_tested),
    collection = collection_class(
      dplyr::pick(dplyr::all_of(paste0("sampling_occurrence_", 1:4)))
    )
  )

# ---- IR Mapper analysis-ready data ----

read_ir_mapper <- function(file) {
  # missing values are blank, NA, NR or NF
  records <- read.csv(file, check.names = FALSE, fileEncoding = "latin1",
                      na.strings = c("", "NA", "NR", "NF"))
  # some column names end in (non-breaking) spaces
  names(records) <- trimws(gsub("\u00a0", " ", names(records)))
  records
}
# records with a species column added
ir_mapper_samples <- function(records) {
  records |>
    dplyr::transmute(
      source = "ir_mapper",
      # one test can be split by species, so the ID includes the species
      record_id = paste("irm", `kdr test ID`, species, sep = "_"),
      citation = `Source 1 citation`,
      month = `Start month`,
      latitude = Latitude,
      longitude = Longitude,
      year_start = `Start year`,
      year_end = `End year`,
      species,
      from_genotypes(dplyr::pick(dplyr::all_of(genotypes$ir_mapper))),
      percent = dplyr::if_else(
        !is.na(`L1014F (%)`) | !is.na(`L1014S (%)`),
        dplyr::coalesce(`L1014F (%)`, 0) + dplyr::coalesce(`L1014S (%)`, 0),
        NA_real_),
      n_tested = dplyr::coalesce(n_genotyped, `No. mosquitoes tested`),
      collection = collection_class(
        dplyr::pick(dplyr::all_of(paste("Capture method", 1:4)))
      )
    )
}

by_species <- read_ir_mapper("data/irmapper/vgsc_species.csv")
by_complex <- read_ir_mapper("data/irmapper/vgsc_complex.csv")
stopifnot(!anyDuplicated(by_complex$`kdr test ID`))
irm_samples <- dplyr::bind_rows(
  by_species |>
    dplyr::mutate(species = tolower(Species)) |>
    dplyr::filter(grepl(gambiae_complex, species)) |>
    ir_mapper_samples(),
  # records for the whole complex whose test isn't in the species file (many
  # are still sums of species records filed under other test IDs, which step
  # 1 below removes)
  by_complex |>
    dplyr::filter(`Complex/Subgroup` == "Gambiae Complex",
                  !`kdr test ID` %in% by_species$`kdr test ID`) |>
    dplyr::mutate(species = "gambiae complex") |>
    ir_mapper_samples()
)

# ---- papers ----

# A key for each sample's paper, the same in both sources where possible. For
# a Vector Atlas sample it is the DOI. For an IR Mapper sample it is the DOI of
# the Vector Atlas paper whose title best matches its reference, or otherwise
# the reference itself. A title matches if at least min_title_match of its
# words (of three or more letters) appear in the reference, and so do at least
# min_title_match of its distinctive words: those in no more than
# max_word_share of the titles (not "anopheles", "resistance" and the like).
# Titles with fewer than min_title_words words (like "PLoS One", given as a
# title) aren't matched
min_title_match <- 0.8
max_word_share <- 0.05
min_title_words <- 5
normalise_doi <- function(x) {
  # without a resolver prefix, a PLoS table or figure suffix (".t002"), or
  # anything but letters, digits, "." and "/" (one DOI is given with spaces
  # where another copy has hyphens)
  x <- sub("^(https?://(dx\\.)?doi\\.org/|doi:)", "", tolower(trimws(x)))
  x <- gsub("[^a-z0-9./]", "", x)
  x <- sub("\\.[tg][0-9]+$", "", x)
  dplyr::if_else(is.na(x) | x == "", NA_character_, x)
}
# the distinct words of at least three letters in each string
words <- function(x) {
  x <- gsub("[^a-z0-9]+", " ",
            tolower(iconv(x, "", "ASCII//TRANSLIT", sub = " ")))
  lapply(strsplit(trimws(x), " +"), function(w) unique(w[nchar(w) >= 3]))
}
# the titles of the Vector Atlas papers with kdr samples
va_titles <- va |>
  dplyr::mutate(doi = normalise_doi(citation_doi)) |>
  dplyr::filter(doi %in% normalise_doi(va_samples$citation),
                !is.na(article_title), article_title != "") |>
  dplyr::distinct(doi, .keep_all = TRUE)
# without notes in square brackets (e.g. "[version 1; referees: 2 approved]")
title_words <- words(gsub("\\[[^]]*\\]", " ", va_titles$article_title))
title_words[lengths(title_words) < min_title_words] <- list(NULL)
word_share <- table(unlist(title_words)) / length(title_words)
common_words <- names(word_share)[word_share > max_word_share]
# the DOI of the Vector Atlas paper matching each reference, or NA
matching_doi <- function(reference) {
  reference_words <- words(reference)
  vapply(reference_words, function(r) {
    share <- vapply(title_words, function(t) {
      distinctive <- setdiff(t, common_words)
      if (length(distinctive) == 0) return(0)
      min(mean(t %in% r), mean(distinctive %in% r))
    }, numeric(1))
    if (max(share) >= min_title_match) va_titles$doi[which.max(share)]
    else NA_character_
  }, character(1))
}
paper_key <- function(source, citation) {
  doi <- rep(NA_character_, length(citation))
  is_va <- source == "vector_atlas"
  doi[is_va] <- normalise_doi(citation[is_va])
  references <- unique(citation[!is_va])
  doi[!is_va] <- matching_doi(references)[match(citation[!is_va], references)]
  dplyr::coalesce(doi, paste("citation:", citation))
}

# ---- combine ----

kdr_samples <- dplyr::bind_rows(va_samples, irm_samples) |>
  dplyr::mutate(paper = paper_key(source, citation)) |>
  # a paper gives allele fractions in the percentage field (e.g. 0.86 for
  # 86%) if its percentages are all at most 1, with some between 0 and 1. In
  # these data that is four papers (Niang 2016 in IR Mapper; Badolo 2012,
  # Ibrahim 2019 and Koffi 2012 in the Vector Atlas), each checked against the
  # paper or the other source's copy
  dplyr::group_by(source, paper) |>
  dplyr::mutate(percent_as_fraction =
                  any(percent > 0 & percent < 1, na.rm = TRUE) &
                  all(percent <= 1, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    percent = dplyr::if_else(percent_as_fraction, percent * 100, percent),
    frequency_from = dplyr::if_else(!is.na(from_counts), "counts", "percent"),
    kdr_frequency = dplyr::coalesce(from_counts, percent / 100)
  ) |>
  dplyr::filter(
    !is.na(kdr_frequency),
    # a few reported percentages exceed 100, which are data-entry errors
    kdr_frequency <= 1,
    !is.na(latitude), !is.na(longitude), !is.na(year_start),
    # a few Vector Atlas records have impossible coordinates (e.g. longitude
    # -42833, presumably a mistyped -4.2833)
    abs(latitude) <= 90, abs(longitude) <= 180
  ) |>
  dplyr::mutate(
    year_end = dplyr::coalesce(year_end, year_start),
    place_and_time = paste(round(latitude, 3), round(longitude, 3), month,
                           year_start, year_end)
  )

# ---- copies of the same sample ----

# drop the records for which drop is TRUE, reporting how many and why
drop_records <- function(data, drop, rule) {
  cat(sum(drop), "records dropped:", rule, "\n")
  data[!drop, ]
}

# 1. totals for the whole complex that are the sum of the species records of
#    the same test (same source, place, start month, years and citation),
#    mostly IR Mapper complex records that sum its species records:
#    the same number tested, or within 5% (or 1) of it with the same pooled
#    frequency (within 0.01), as when a few mosquitoes weren't identified
kdr_samples$test <- with(kdr_samples, paste(source, place_and_time, citation))
species_sums <- kdr_samples |>
  dplyr::filter(species != "gambiae complex") |>
  dplyr::group_by(test) |>
  dplyr::summarise(f_species = sum(n_tested * kdr_frequency) / sum(n_tested),
                   n_species = sum(n_tested))
complex_totals <- kdr_samples |>
  dplyr::filter(species == "gambiae complex") |>
  dplyr::inner_join(species_sums, by = dplyr::join_by(test)) |>
  dplyr::filter(n_tested == n_species |
                  abs(n_tested - n_species) <= pmax(1, 0.05 * n_species) &
                  abs(kdr_frequency - f_species) <= 0.01)
kdr_samples <- drop_records(
  kdr_samples, kdr_samples$record_id %in% complex_totals$record_id,
  "complex total that sums its species records")

# 2. repeats within a source: the same paper, place, start month, years,
#    species, number tested and frequency (keeping one with genotype counts,
#    if any)
kdr_samples <- kdr_samples |>
  dplyr::arrange(frequency_from != "counts") |>
  dplyr::mutate(copy = paste(source, paper, place_and_time, species, n_tested,
                             round(kdr_frequency, 4)))
kdr_samples <- drop_records(kdr_samples, duplicated(kdr_samples$copy),
                            "repeat within a source")

# 3. a paper entered twice in the Vector Atlas: within a paper, species,
#    start year and start month, a record within 0.05 degrees of an earlier
#    one with the same number tested and frequency (within 0.01). The number
#    tested must be equal where both come from genotype counts (records of
#    the two kdr mutations typed in the same mosquitoes can otherwise match),
#    and within 1 otherwise. Frequencies of 0 or 1 are left alone, since
#    separate small samples often share them
is_entered_twice <- function(latitude, longitude, n, frequency, counted) {
  copy <- rep(FALSE, length(n))
  for (i in seq_along(n)[-1]) {
    earlier <- seq_len(i - 1)
    copy[i] <- frequency[i] > 0 && frequency[i] < 1 &&
      any(!copy[earlier] &
            abs(latitude[earlier] - latitude[i]) <= 0.05 &
            abs(longitude[earlier] - longitude[i]) <= 0.05 &
            abs(n[earlier] - n[i]) <= 1 - (counted[earlier] & counted[i]) &
            abs(frequency[earlier] - frequency[i]) <= 0.01, na.rm = TRUE)
  }
  copy
}
kdr_samples <- kdr_samples |>
  dplyr::group_by(source, paper, species, year_start, month) |>
  dplyr::mutate(entered_twice = if (dplyr::first(source) == "vector_atlas") {
    is_entered_twice(latitude, longitude, n_tested, kdr_frequency,
                     frequency_from == "counts")
  } else {
    FALSE
  }) |>
  dplyr::ungroup()
kdr_samples <- drop_records(kdr_samples, kdr_samples$entered_twice,
                            "a paper entered twice in the Vector Atlas")

# 4. papers in both sources: the version with more distinct samples the model
#    can use (not An. arabiensis; distinct in species, year, number tested and
#    frequency, so that a sample entered twice counts once) is kept, IR
#    Mapper's on a tie. Each source can lack samples the other has (e.g. IR
#    Mapper having only a paper's An. arabiensis samples)
usable <- kdr_samples |>
  dplyr::filter(!grepl("arabiensis", species)) |>
  dplyr::distinct(paper, source, species, year_start, n_tested,
                  frequency = round(kdr_frequency, 2)) |>
  dplyr::count(paper, source, name = "samples")
ir_mapper <- dplyr::filter(usable, source == "ir_mapper")
vector_atlas <- dplyr::filter(usable, source == "vector_atlas")
in_both <- intersect(unique(kdr_samples$paper[kdr_samples$source ==
                                                "ir_mapper"]),
                     unique(kdr_samples$paper[kdr_samples$source ==
                                                "vector_atlas"]))
versions <- data.frame(paper = in_both) |>
  dplyr::mutate(
    ir_mapper = dplyr::coalesce(
      ir_mapper$samples[match(paper, ir_mapper$paper)], 0),
    vector_atlas = dplyr::coalesce(
      vector_atlas$samples[match(paper, vector_atlas$paper)], 0),
    kept = dplyr::if_else(vector_atlas > ir_mapper, "vector_atlas",
                          "ir_mapper")
  )
cat("papers in both sources, by version kept:\n")
print(table(versions$kept))
cat("papers whose Vector Atlas version is kept:\n")
versions |>
  dplyr::filter(kept == "vector_atlas") |>
  print()
kept_version <- versions$kept[match(kdr_samples$paper, versions$paper)]
kdr_samples <- drop_records(
  kdr_samples, !is.na(kept_version) & kdr_samples$source != kept_version,
  "the other source's version of a paper in both")

kdr_samples <- drop_records(kdr_samples, kdr_samples$n_tested %in% 0,
                            "0 mosquitoes recorded")

# the number of mosquitoes: a recorded number rounded (a few are not whole
# numbers), or a missing one inferred around the median recorded number,
# within 40% either side of it (kdr samples are often small subsamples, so
# not ir_cube's 100 for bioassays)
recorded <- kdr_samples$n_tested[!is.na(kdr_samples$n_tested)]
typical_n <- round(median(recorded))
kdr_samples <- kdr_samples |>
  dplyr::mutate(n_tested = dplyr::if_else(
    is.na(n_tested),
    infer_sample_size(kdr_frequency, typical_n, round(typical_n * c(0.6, 1.4))),
    round(n_tested)
  )) |>
  dplyr::select(source, record_id, citation, paper, latitude, longitude,
                year_start, year_end, species, n_tested, kdr_frequency,
                collection)

# samples of An. arabiensis are left out, because its kdr effect differs from
# that of An. gambiae and An. coluzzii. They are kept until here because
# finding copies (step 1) and inferring missing numbers of mosquitoes use them
kdr_samples <- drop_records(kdr_samples,
                            grepl("arabiensis", kdr_samples$species),
                            "An. arabiensis")

# ---- pixel, species group and year ----

# species group of each sample, from the species labels of the two sources:
# An. coluzzii is the M form and An. gambiae the S form
species_group_of <- function(species) {
  dplyr::case_when(
    grepl("coluzzii/gambiae|s form_m form", species) ~ "coluzzii and gambiae",
    grepl("coluzzii|m form", species) ~ "coluzzii",
    grepl("gambiae_s form|anopheles gambiae", species) ~ "gambiae",
    .default = sub("anopheles ", "", species)
  )
}

# each sample's pixel (pixel_id_at(), as for bioassays). The collection year
# is the middle of the sample's collection years, rounded, and a study-year
# is a paper's samples of one species group in one collection year
kdr_samples <- kdr_samples |>
  dplyr::mutate(
    pixel_id = pixel_id_at(longitude, latitude),
    species_group = species_group_of(species),
    collection_year = round((year_start + year_end) / 2),
    study_year = paste(species_group, paper, collection_year)
  )

stopifnot(
  !anyDuplicated(kdr_samples$record_id),
  all(kdr_samples$kdr_frequency >= 0),
  !anyNA(kdr_samples$collection),
  !anyNA(kdr_samples$pixel_id)
)

write.csv(kdr_samples, "data/kdr_samples.csv", row.names = FALSE)

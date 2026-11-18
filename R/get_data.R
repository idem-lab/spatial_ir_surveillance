# Gather the input data. The data files are not in the repo (data/ is
# gitignored). Each file listed below that is missing from data/ is copied
# from a local clone of ir_cube (https://github.com/idem-lab/ir_cube) at
# ir_cube_dir on Nick's machine; without a clone, ask Nick for the files and
# put them where listed. The script then stops if any file's md5 sum differs
# from the one recorded, so that everyone runs the analysis on the same data
# and a changed file can't slip in unnoticed. To update a file deliberately,
# update its source and checksum here.
#
# Inputs:
#   - from ir_cube's data prep at commit
#     2e7e4a62a35bc8905b9213f5e7384849b251cb9f: the cleaned bioassays (made by
#     ir_cube's R/prep_bioassays.R) and the 0.0417 degree (about 5 km) grid
#     of pixels it uses (made by ir_cube's R/prep_rasters.R)
#   - raw downloads kept in the ir_cube clone but not in that repo: the Vector
#     Atlas full data download of 21 July 2025, with Vgsc (kdr) genotype
#     counts and allele frequencies, and IR Mapper's analysis-ready Vgsc data
#     (Moyes et al. 2019, Scientific Data): allele frequencies by species and
#     for the whole complex, and genotype counts in bioassay survivors and
#     dead, in matched pairs

ir_cube_dir <- "~/Dropbox/github/ir_cube"
ir_mapper_dir <- "data/raw/IRMapper Analysis ready data sets"

# each file's place in data/, its source in the ir_cube clone, and md5 sum
inputs <- data.frame(
  file = c(
    "data/ir_cube/all_gambiae_complex_data.RDS",
    "data/ir_cube/raster_mask.tif",
    "data/vector_atlas/va_full_data.csv",
    "data/irmapper/vgsc_species.csv",
    "data/irmapper/vgsc_complex.csv",
    "data/irmapper/vgsc_survivors_dead.csv"
  ),
  source = c(
    "data/clean/all_gambiae_complex_data.RDS",
    "data/clean/raster_mask.tif",
    "data/raw/VA_FULL_DATA_20250721.csv",
    file.path(ir_mapper_dir, "7_Vgsc-allele-freq_species.csv"),
    file.path(ir_mapper_dir, "6_Vgsc-allele-freq_complex-subgroup.csv"),
    file.path(ir_mapper_dir, "8_Vgsc-allele-freq_survivors-dead.csv")
  ),
  md5 = c(
    "44772bd1da61751c30d7e6e1f020eaf7",
    "54e89d34918b62c0296cfec92c890585",
    "595a850085e6bf9c5f0843718fa48abb",
    "0f82466c6ed88c6b2bcfe979e5f3a857",
    "702378a6670a91a198fd71bffb0932fc",
    "77d857dd182eee79b3ca131bfafe7b70"
  )
)

for (i in seq_len(nrow(inputs))) {
  file <- inputs$file[i]
  if (!file.exists(file)) {
    source_file <- file.path(ir_cube_dir, inputs$source[i])
    if (!file.exists(source_file)) {
      stop(file, " is missing and there is no ir_cube clone at ",
           ir_cube_dir, ". Ask Nick for the file")
    }
    dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
    file.copy(source_file, file)
  }
  if (tools::md5sum(file) != inputs$md5[i]) {
    stop(file, " does not match the version this analysis was run with")
  }
}

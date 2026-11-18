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

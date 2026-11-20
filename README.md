# spatial_ir_surveillance

Simulation study comparing sentinel-site bioassays with spatially dispersed
genetic sampling for estimating district-level insecticide resistance in African
malaria vectors.

Work is tracked in the repository's milestones and issues.

## Running the code

Data files are not in the repository (`data/` is gitignored). Before running
anything, ask Nick for these files and put them where shown. The first two
come from [ir_cube](https://github.com/idem-lab/ir_cube); the rest from the
Vector Atlas and IR Mapper. `R/get_data.R` lists their sources and checksums.

- `data/ir_cube/all_gambiae_complex_data.RDS`
- `data/ir_cube/raster_mask.tif`
- `data/vector_atlas/va_full_data.csv`
- `data/irmapper/vgsc_species.csv`
- `data/irmapper/vgsc_complex.csv`
- `data/irmapper/vgsc_survivors_dead.csv`

Open the RStudio project and run `renv::restore()` to install the packages at
the versions recorded in `renv.lock`. Run it again whenever `renv.lock`
changes. Then run the scripts from the project root, in this order:

1. `R/get_data.R`: checks the data files match the versions the analysis was
   run with, copying them from a local clone of ir_cube if they are missing.
2. `R/prep_bioassays.R`: subsets the bioassays and adds district and pixel IDs.
   Downloads admin boundaries from GADM on the first run.
3. `R/prep_kdr_samples.R`: kdr allele frequencies in genotyped samples.
4. `R/fit_kdr_predictor.R`: a Bayesian model (fitted with greta) of
   pyrethroid bioassay survival as a function of kdr allele frequency, with
   posterior predictive checks. greta is installed from its GitHub repository
   (the version in `renv.lock`) and needs Python and TensorFlow: run
   `greta::greta_set_python("uv")` once before the first run.
5. `R/fit_bioassay_variance.R`: a Bayesian model (fitted with greta) of how
   much bioassay mortality varies between pixels within a district in one
   year, after the variation between districts and over time, with
   posterior predictive checks and districts simulated from its priors and
   fit. Some bioassays are excluded after checking them against their
   sources: `metadata/bioassay_exclusions.csv` gives the decision and reason
   for each.

Functions used by these scripts are in `R/functions.R`.

**How the code is written.** Much of the code here is written with AI coding
agents and checked locally by a person before it is pushed. Agent involvement
shows as a co-author on each commit. You're welcome to contribute however suits
you: comments, issues and pull requests, with or without an agent.

# spatial_ir_surveillance

Simulation study comparing sentinel-site bioassays with spatially dispersed
genetic sampling for estimating district-level insecticide resistance in African
malaria vectors.

Work is tracked in the repository's milestones and issues.

## Running the code

Data files are not in the repository (`data/` is gitignored). Before running
anything, ask Nick for two files made by
[ir_cube](https://github.com/idem-lab/ir_cube) and put them in `data/ir_cube/`:

- `all_gambiae_complex_data.RDS`
- `raster_mask.tif`

Open the RStudio project and run `renv::restore()` to install the packages at
the versions recorded in `renv.lock`. Run it again whenever `renv.lock`
changes. Then run the scripts from the project root, in this order:

1. `R/prep_bioassays.R`: subsets the bioassays and adds district and pixel IDs.
   Downloads admin boundaries from GADM on the first run.

**How the code is written.** Much of the code here is written with AI coding
agents and checked locally by a person before it is pushed. Agent involvement
shows as a co-author on each commit. You're welcome to contribute however suits
you: comments, issues and pull requests, with or without an agent.

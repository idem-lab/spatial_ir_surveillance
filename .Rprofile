# Keep the renv package library outside the project, so that Dropbox doesn't
# sync it. This must come before renv is activated.
Sys.setenv(
  RENV_PATHS_LIBRARY_ROOT = file.path(tools::R_user_dir("renv", "data"),
                                      "library")
)
source("renv/activate.R")

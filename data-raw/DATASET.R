## code to prepare `DATASET` dataset goes here

# Mass is corrected by 10^omega to convert units from hours to seconds and 
# mass from kg to g.
omega <- 5.945

usethis::use_data(omega, overwrite = TRUE, internal = TRUE)

################################################
# Install packges (run once)
install.packages(c("devtools",
                   "tidyverse",
                   "activity",
                   "Distance"))
devtools::install_github("inbo/camtraptor")

################################################
# Load packages (run before each session)
library(tidyverse)
library(camtraptor)
library(activity)
library(Distance)
devtools::source_url("https://raw.githubusercontent.com/MarcusRowcliffe/camtraptor_density_example/main/rem_functions.R")

################################################
# Load data
package <- read_camtrap_dp2("./example_data/datapackage.json")

################################################
# One step analysis
result <- rem_estimate(package)
# Check outputs
result$estimates
result$species
result$data
plot(result$activity_model)
plot(result$radius_model, pdf=TRUE)
plot(result$angle_model)

################################################
# Building the analysis yourself
species <- "Vulpes vulpes"
# Fit auxilaiary parameter models
spd <- fit_speedmodel(package, species=species)
act <- fit_actmodel(package, species=species)
rad <- fit_detmodel(radius~1, package, species=species,
                    order=0, transect="point", truncation=10)
ang <- fit_detmodel(angle~1, package, species=species, order=0)
# Examine models
plot(rad, pdf=TRUE)
plot(ang)
plot(act)
# Generate trap rate and parameter data tables
trdata <- get_rem_data(package, species)
params <- get_parameter_table(rad, ang, spd, act)
# Estimate density
rem(trdata, params)

################################################
# Install packges (run once)
install.packages(c("devtools",
                   "tidyverse",
                   "jpeg",
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
package <- read_camtrap_dp("./data/datapackage.json")

################################################
# One step analysis
result <- rem_estimate(package)
# Inspect outputs
result$estimates
result$species
result$data
plot(result$activity_model)
plot(result$radius_model, pdf=TRUE)
plot(result$angle_model)

################################################
# Building the analysis yourself
species <- "Sus scrofa"
# Fit auxiliary parameter models
spd <- fit_speedmodel(package, species=species)
act <- fit_actmodel(package, species=species)
rad <- fit_detmodel(radius~1, package, species=species, order=0, truncation=10)
ang <- fit_detmodel(angle~1, package, species=species, order=0, unit="radian")
# Examine models
plot(rad, pdf=TRUE)
plot(ang)
plot(act)
# Generate trap rate and parameter data tables
data <- get_rem_data(package, species)
param <- get_parameter_table(rad, ang, spd, act)
# Estimate density
rem(data, param)


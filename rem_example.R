library(tidyverse)
library(camtraptor)
library(activity)
library(Distance)

dir <- "C:/Users/rowcliffe.m/Downloads/camera-model-calibration-20221114172535"
package <- read_camtrap_dp(file.path(dir, "datapackage.json"))
package$data$observations
View(package$data$observations)
package$data$observations <- rename(package$data$observations, 
                                    speed=X22, 
                                    radius=X23, 
                                    angle=X24)

est <- rem_estimate(pkg, check_deployments=FALSE)
est$estimates
est$species
est$data

species="Vulpes vulpes"
spd <- fit_speedmodel(pkg, species=sp)
act <- fit_actmodel(pkg, species=sp)
rad <- fit_detmodel(radius~1, pkg, species=sp,
                    order=0, transect="point", truncation=10)
ang <- fit_detmodel(angle~1, pkg, species=sp, order=0)

plot(rad, pdf=TRUE)
plot(est$radius_model, pdf=TRUE)
plot(est$angle_model)
plot(est$activity_model)

param <- head(est$estimates, -2)[,1:3]
param[2, 2:3] <- param[2, 2:3]*pi/180
rem(est$data, param)
write.csv(est$data, "rem_data.csv", row.names = FALSE)
write.csv(param, "rem_param.csv", row.names = FALSE)

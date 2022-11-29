#' Read a camera trap datapackage
#'
#' A temporary replacement for camtraptor::read_camtrap_dp to read
#' data packages that contain animal position data.
#' 
#' @param file Path or URL to a datapackage.json file.
#' @return As for camtraptor::read_camtrap_dp.
read_camtrap_dp2 <- function(file){
  package <- read_camtrap_dp(file)
  package$data$observations <- rename(package$data$observations, 
                                      speed=X22, 
                                      radius=X23, 
                                      angle=X24)
  package
}


#' Select a species name
#'
#' Presents a table of species names with observation count for each
#' and allows the user to interactively select one.
#'
#' @param package Camera trap data package object, as returned by
#'   `read_camtrap_dp()`.
#' @return A character string, scientific species name.
#' @family density estimation functions
#' @export
select_species <- function(package){
  tab <-  get_species(package)
  tab <- tab[, grepl("Name", names(tab))]
  n <- table(package$data$observations$scientificName)
  tab$n_observations <- n
  
  print.data.frame(tab)
  i <- NA
  while(is.na(i) || i<1 || i>nrow(tab))
    i <- readline("Enter row number of species to analyse: ") %>%
    as.numeric() %>%
    suppressWarnings()
  as.character(tab$scientificName[i])
}

#' Check deployment calibration models
#' 
#' Displays deployment calibration model diagnostic plots and allows 
#' users to record interactively whether each deployment is reliable.
#' 
#' @param package Camera trap data package object, as returned by
#'   `read_camtrap_dp()`.
#' @return The original package with logical column `useDeployment`
#'  added to deployments and observations data.
#' @family density estimation functions
#' @export
check_deployment_models <- function(package){
  plot_folder <- file.path(package$directory, "positioning_plots")
  if(!dir.exists(plot_folder)) stop("The specified folder does not exist.")
  plot_dirs <- list.dirs(plot_folder, recursive=FALSE)
  plots <- file.path(plot_dirs, "ratio.jpeg") %>%
    lapply(jpeg::readJPEG)
  names(plots) <- basename(plot_dirs)
  
  depdat <- package$data$deployments %>%
    dplyr::select(deploymentID, locationName)
  depdat$useDeployment <- FALSE
  
  for(i in 1:nrow(depdat)){
    dep <- depdat$deploymentID[i]
    if(dep %in% names(plots)){
      img <- plots[[dep]]
      imdim <- dim(img)
      p <- ggplot() + annotation_raster(img, 1, imdim[2], 1, imdim[1]) + 
        xlim(1, imdim[2]) + ylim(1, imdim[1]) +
        theme_void() + ggtitle(depdat$locationName[i])
      print(p)
      answer <- NA
      while(is.na(answer) || !answer %in% c("y", "n"))
        answer <- readline("Use deployment? (y/n): ") %>% tolower()
      if(answer=="y") depdat$useDeployment[i] <- TRUE
    }
  }
  
  if("useDeployment" %in% names(package$data$observations))
    package$data$observations <- select(package$data$observations, -useDeployment)
  package$data$observations <- package$data$observations %>%
    left_join(depdat, by="deploymentID")
  
  package$data$deployments$useDeployment <- depdat$useDeployment
  package
}

#' Estimate average animal speed
#'
#' Calculates harmonic mean and standard error of animal speed while active
#' from a data package.
#'
#' @param package Camera trap data package object, as returned by
#'   `read_camtrap_dp()`.
#' @param species A character string indicating species subset to analyse; user
#'   input required if NULL.
#' @return List with elements:
#'  - speed: a one row dataframe containing mean (estimate) and standard error 
#'    (se) speed while active
#'  - data: a numeric vector of the data from which the estimate is derived
#' @family density estimation functions
#' @export
fit_speedmodel <- function(package, 
                           species=NULL, 
                           distUnit=c("m", "km", "cm"),
                           timeUnit=c("second", "minute", "hour", "day")){
  distUnit <- match.arg(distUnit)
  timeUnit <- match.arg(timeUnit)
  
  if(is.null(species)) species <- select_species(package)
  obs <- package$data$observations %>%
    subset(scientificName==species & speed > 0.01 & speed < 10)
  if("useDeployment" %in% names(obs)) obs <- subset(obs, useDeployment)
  mn <- 1/mean(1/obs$speed, na.rm=FALSE)
  se <- mn^2 * sqrt(var(1/obs$speed, na.rm=FALSE)/nrow(obs))
  list(speed=data.frame(estimate=mn, se=se), 
       data=obs$speed,
       distUnit=distUnit,
       timeUnit=timeUnit)
}

#' Fit an activity model
#'
#' Fits an activity model to data package data and estimates activity 
#' level (proportion of time spent active).
#'
#' @param package Camera trap data package object, as returned by
#'   `read_camtrap_dp()`.
#' @param species A character string indicating species subset to analyse; user
#'   input required if NULL.
#' @param obsdef Observation definition, either individual or sequence.
#' @param reps Number of bootstrap replicates to run.
#' @param ... Arguments passed to fitact.
#' @return An `actmod` list.
#' @seealso \code{\link{activity::fitact}}
#' @family density estimation functions
#' @export
fit_actmodel <- function(package, 
                         species=NULL, 
                         reps=999,
                         obsdef=c("individual", "sequence"),
                         ...){
  obsdef <- match.arg(obsdef)
  if(is.null(species)) species <- select_species(package)
  deps <- package$data$deployments
  obs <- package$data$observations %>%
    dplyr::filter(scientificName==species) %>%
    dplyr::select(deploymentID, sequenceID, timestamp, count)
  i <- switch(obsdef,
              individual = rep(1:nrow(obs), obs$count),
              sequence = !duplicated(obs$sequenceID))
  obs <- obs[i, ]

  if(nrow(obs)>1){ 
    obs <- deps %>%
      dplyr::select(deploymentID, latitude, longitude) %>%
      dplyr::right_join(obs, by="deploymentID") %>%
      dplyr::select(-count)
    suntimes <- insol::daylength(obs$latitude, obs$longitude, 
                                 insol::JD(obs$timestamp), 0)
    timeshift <- pi - mean(suntimes[, 1] + suntimes[,3]/2) * pi/12
    obs$solartime <- obs %>%
      with(activity::solartime(timestamp, latitude, longitude, 0)) %>%
      .$solar %>%
      + timeshift %>%
      activity::wrap()
      activity::fitact(obs$solartime,
                       adj = 1.5, sample = "data", reps = reps, ...)
  } else
    NULL
}

#' Fit a detection function model
#' 
#' Fits a detection function to a data package and estimates effective 
#' detection distance (EDD).
#' 
#' @param formula A two sided formula relating radius or angle data 
#'   to covariates.
#' @param package Camera trap data package object, as returned by
#'   `read_camtrap_dp()`.
#' @param species A character string indicating species subset to analyse; user
#'   input required if NULL.
#' @param ... Arguments passed to ds.
#' @return A `ddf` detection function model list, with additional element
#'   `edd`, a vector with estimated and standard error effective detection 
#'   distance, or the `newdata` dataframe with EDD estimate and se added.
#' @seealso \code{\link{Distance::ds}}
#' @family density estimation functions
#' @export
fit_detmodel <- function(formula, 
                         package, 
                         species=NULL, 
                         newdata=NULL,
                         unit=c("m", "km", "cm", "degree", "radian"),
                         ...){
  unit <- match.arg(unit)
  
  # get and check model variables
  allvars <- all.vars(formula)
  depvar <- allvars[1]
  covars <- tail(allvars, -1)
  data <- package$data$observations
  if(!all(allvars %in% names(data))) stop("Can't find all model variables in data")
  if("distance" %in% covars) stop("Cannot use \"distance\" as a covariate name - rename and try again")
  
  # set up data
  if(is.null(species)) species <- select_species(package)
  data <- data %>%
    subset(scientificName==species) %>%
    dplyr::select(all_of(allvars)) %>%
    tidyr::drop_na() %>%
    as.data.frame()
  if("useDeployment" %in% names(data)) data <- subset(data, useDeployment)
  
  classes <- dplyr::summarise_all(data, class)
  if(classes[depvar]=="numeric"){
    data <- data %>%
      dplyr::rename(distance=all_of(depvar)) %>%
      dplyr::mutate(distance=abs(distance))
  } else{
    cats <- strsplit(as.character(dplyr::pull(data, depvar)), "-")
    data$distbegin <- unlist(lapply(cats, function(x) as.numeric(x[1])))
    data$distend <- unlist(lapply(cats, function(x) as.numeric(x[2])))
    data$distance <- (data$distbegin + data$distend) / 2
  }
  
  # model fitting
  type <- if(unit %in% c("m", "km", "cm")) "point" else "line"
  args <- c(data=list(data), formula=formula[-2], transect=type, list(...))
  mod <- suppressWarnings(suppressMessages(do.call(ds, args)$ddf))
  
  # esw prediction
  if(length(covars)==0) 
    newdata <- data.frame(x=0) else{
      if(is.null(newdata)){
        newdata <- data %>% dplyr::select(all_of(covars)) %>%
          lapply(function(x) 
            if(is.numeric(x)) mean(x, na.rm=T) else sort(unique(x)))  %>%
          expand.grid()
      } else{
        if(!all(covars %in% names(newdata))) stop("Can't find all model covariates in newdata")
      }}
  prdn <- predict(mod, newdata, esw=TRUE, se.fit=TRUE)
  if(mod$meta.data$point){
    prdn$se.fit <- 0.5 * prdn$se.fit / (pi * prdn$fitted)^0.5
    prdn$fitted <- sqrt(prdn$fitted/pi)
  }
  ed <- cbind(estimate=prdn$fitted, se=prdn$se.fit)
  if(length(covars)>=1) ed <- cbind(newdata, ed)
  mod$edd <- ed
  mod$unit <- unit
  mod
}

#' Get REM data from a camtrap-dp datapackage
#'
#' Extracts a data table of observation counts and effort for each
#' camera location in a camtrap-dp data package.
#' 
#' @param package Camera trap data package object, as returned by
#'   `read_camtrap_dp()`.
#' @param species A character string indicating species subset to extract
#'   data for; user input required if NULL.
#' @param unit The time unit in which to return camera effort.
#' @return A tibble with columns:
#'   - locationName: name of the camera location
#'   - effort: the camera time for the location
#'   - unit: the effort time unit
#'   - scientificName: the scientific name of the species data extracted
#'   - n: the observation counts
#' @family density estimation functions
#' @export
get_rem_data <- function(package, species=NULL, 
                         unit=c("second", "minute", "hour", "day")){
  unit <- match.arg(unit)
  if(is.null(species)) species <- select_species(package)
  dep <- package$data$deployments %>%
    dplyr::select(deploymentID, locationName)
  eff <- package %>%
    get_effort(unit=unit) %>%
    dplyr::select(deploymentID, effort)
  res <- package %>%
    get_n_individuals(species=species) %>%
    suppressMessages() %>%
    dplyr::left_join(dep, by="deploymentID") %>%
    dplyr::left_join(eff, by="deploymentID") %>%
    dplyr::group_by(locationName) %>%
    dplyr::summarise(n = sum(n), effort=sum(effort))
  res$effort_unit <- unit
  res$species <- species
  res
}

#' Get a unit multiplier
#' 
#' Gives the value by which to multiply a value in one unit to another in
#' one of three types: distance, time and angle.
#'
#' @param unitIN A text value giving the unit of the input; must be one of
#'   "cm", "m", "km" for distances, "second", "minute", "hour", "day" for 
#'   times, or "radian" "degree for angles.
#' @param unitOUT The same for output; must be of the same type as unitIN.
#' @return A number giving the amount by which to multipy input values
#'   to arrive a unit-converted output.
get_multiplier <- function(unitIN, unitOUT){
  dunits <- c("cm", "m", "km")
  dmult <- c(1, 1e2, 1e5)
  tunits <- c("second", "minute", "hour", "day")
  tmult <- c(1, 60, 60^2, 24*60^2)
  aunits <- c("radian", "degree")
  amult <- c(1, pi/180)
  if(unitIN %in% dunits & unitOUT %in% dunits){
    u <- dunits
    m <- dmult
    n <- length(dunits)
  } else
    if(unitIN %in% tunits & unitOUT %in% tunits){
      u <- tunits
      m <- tmult
      n <- length(tunits)
    } else
      if(unitIN %in% aunits & unitOUT %in% aunits){
        u <- aunits
        m <- amult
        n <- length(aunits)
      } else
        stop("Units not of the same type or not recognised")
  
  tab <- data.frame(from = rep(u, each=n),
                    to = rep(u, n),
                    mult = rep(m, each=n) / rep(m, n))
  tab$mult[tab$from==unitIN & tab$to==unitOUT]
}

#' Change the units of an REM parameter table
#' 
#' Changes the units of parameters from their current setting to new 
#' user-definedunits.
#' 
#' @param param An REM parameter dataframe with columns: estimate, se and unit,
#'   and rows named with parameter names in "radius", "angle", "speed".
#' @param units A named vector of character unit names giving the units to which
#'   each parameter is to be converted. 
#' @return A new parameter dataframe with estimate, se and unit values modified
#'   to change units from current to units
re_unit <- function(param,
                    units = c(radius="m",
                              angle="degree",
                              speedDist="km",
                              speedTime="hour")){

  spd_units <- unlist(strsplit(param["speed", "unit"], "/"))
  rm <- get_multiplier(param["radius", "unit"], units["radius"])
  am <- get_multiplier(param["angle", "unit"], units["angle"])
  sdm <- get_multiplier(spd_units[1], units["speedDist"])
  stm <- get_multiplier(spd_units[2], units["speedTime"])

  # Modify table
  j <- c("estimate", "se")
  param["radius", j] <- param["radius", j] * rm
  param["angle", j] <- param["angle", j] * am
  param["speed", j] <- param["speed", j] * sdm / stm
  param["radius", "unit"] <- units["radius"]
  param["angle", "unit"] <- units["angle"]
  param["speed", "unit"] <- paste(units["speedDist"], units["speedTime"], sep="/")
  param
}

#' Create a parameter table from a set of models
#' 
#' Creates a table of REM parameters taken from models for detection radius, 
#' angle, speed and activity level. 
#' 
#' @param radius_model A detection radius model fitted using fit_detmodel
#' @param angle_model A detection angle model fitted using fit_detmodel
#' @param speed_model A speed model fitted using fit_speedmodel
#' @param speed_model An activity  model fitted using fit_actmodel
#' @param units A named vector of character unit names giving the units to which
#'   each parameter is to be converted. 
#' @return A dataframe of parameter estimates with columns: estimate, se and
#'   unit.
#' @family density estimation functions
#' @export
get_parameter_table <- function(radius_model, angle_model, 
                                speed_model, activity_model=NULL,
                                units = c(radius="m",
                                          angle="degree",
                                          speedDist="km",
                                          speedTime="hour")){
  expectedNames <- c("radius", "angle", "speedDist", "speedTime")
  expectedDist <- c("cm", "m", "km")
  expectedTime <- c("second", "minute", "hour", "day")
  expectedAngle <- c("radian", "degree")
  if(!all(expectedNames %in% names(units)))
    stop(paste("Names of units must include all of:",
     paste(expectedNames, collapse = ", ")))
  if(!all(units[c("radius", "speedDist")] %in% expectedDist))
    stop(paste("Distance units must be one of: ",
               paste(expectedDist, collapse = ", ")))
  if(!all(units["speedTime"] %in% expectedTime))
    stop(paste("Time units must be one of: ",
               paste(expectedTime, collapse = ", ")))
  if(!all(units["angle"] %in% expectedAngle))
    stop(paste("Angle units must be one of: ",
               paste(expectedAngle, collapse = ", ")))

  act_val <- if(is.null(activity_model)) c(1,0) else activity_model@act[1:2]
  res <- data.frame(rbind(radius_model$edd, 
                          angle_model$edd * 2,
                          speed_model$speed, 
                          act_val),
                    unit = c(radius_model$unit, 
                             angle_model$unit, 
                             paste(speed_model$distUnit,
                                   speed_model$timeUnit, 
                                   sep="/"),
                             "none"))
  rownames(res) <- c("radius", "angle", "speed", "activity")
  re_unit(res, units)
}

#' Harmonise parameter units
#' 
#' Changes REM parameter values and units to ensure that time and distance 
#' are expressed in consistent units across parameters.
#' 
#' @param param An REM parameter dataframe with columns: estimate, se and unit,
#'   and rows named with parameter names in "radius", "angle", "speed", 
#'   obtained using get_parameter_table.
#' @param data An rem data object obtained using get_rem_data.
#' @param distUnit A character defining the ditance unit to which to harmonise
#' @param timeUnit A character defining the time unit to which to harmonise
#' @family density estimation functions
#' @export
harmonise_units <- function(param, 
                            data,
                            distUnit=c("km","m", "cm"), 
                            timeUnit=c("day", "hour", "minute", "second")){
  distUnit <- match.arg(distUnit)
  timeUnit <- match.arg(timeUnit)

  spd_units <- unlist(strsplit(param["speed", "unit"], "/"))
  units <- c(radius=distUnit, angle="radian", 
                speedDist=distUnit, speedTime=timeUnit)
  paramOUT <- re_unit(param, units)

  dm <- get_multiplier(data$effort_unit[1], timeUnit)
  data$effort <- data$effort * dm
  data$effort_unit <- timeUnit
  list(data=data, param=paramOUT)
}


#' Fit a random encounter model
#'
#' Estimates REM density given dataframes of trap rate and auxiliary 
#' parameter data
#'
#' @param data A dataframe containing a row per sampling location and columns:
#'  - observations: the number of animal contact events
#'  - effort: the amount of camera time
#'  If `stratum_areas` is provided, additional column required:
#'  - stratumID: key identifying which stratum each location sits in
#' @param param A dataframe containing REM parameter estimates with columns
#'  parameter (parameter name), estimate(parameter standard error) and se 
#'  (parameter standard error); use one row per parameter, with the following
#'  names:
#'  Mandatory
#'  - radius: effective detection radius
#'  - angle: effective detection angle
#'  - speed: average animal speed
#'  Optionally
#'  - activity: activity level (proportion of time spent active)
#'  If activity is provided, speed is assumed to be average speed while active,
#'  otherwise it is taken to be day range (distance traveled per day)
#' @param stratum_areas A dataframe with one row per stratum and columns:
#' - stratumID: stratum ID key, matched with the same key in data
#' - area: stratum areas (or proportional coverage of the study area)
#' @param reps Number of bootstrap replicates for error estimation. 
#' @return A dataframe with the original parameters plus trap rate and density
#'  estimates and standard errors.
#' @param ... Arguments passed to harmonise_units.
#' @details The function makes no assumptions about units. It is up to the user to ensure 
#'  that these are harmonised across data and parameters.
#' @family density estimation functions
#' @export
rem <- function(data, param, stratum_areas=NULL, reps=999, ...){
  
  traprate <- function(data){
    if(is.null(stratum_areas)){
      sum(data$n) / sum(data$effort)
    } else{
      local_density <- sapply(stratum_areas$stratumID, function(stratum){
        i <- data$stratumID==stratum
        sum(data$n[i]) / sum(data$effort[i])
      })
      sum(local_density * stratum_areas$area) / sum(stratum_areas$area)
    }
  }
  
  sampled_traprate <- function(){
    i <- if(is.null(stratum_areas)) 
      sample(1:nrow(data), replace=TRUE) else
        as.vector(sapply(stratum_areas$stratumID, function(stratum){
          sample(which(data$stratumID==stratum), replace=TRUE)
        }))
    traprate(data[i, ])
  }
  
  if(!all(c("effort", "n") %in% names(data)))
    stop("data must contain (at least) columns effort and observations")
  if(!all(c("speed", "radius", "angle") %in% rownames(param)))
    stop("param must contain (at least) parameters speed, radius and angle")
  if(!is.null(stratum_areas)){
    if(!"stratumID" %in% names(data))
      stop("data must contain column stratumID for stratified analysis")
    if(!all(c("stratumID", "area") %in% names(stratum_areas)))
      stop("stratum_areas must contain columns stratumID and area")
    if(!all(data$stratumID %in% stratum_areas$stratumID)) 
      stop("Not all strata in data are present in stratum_areas")
  }  
  
  if(!"activity" %in% rownames(param)) 
    param <- rbind(param, activity=data.frame(estimate=1, se=0))
  pkg <- harmonise_units(param, data, ...)
  param <- pkg$param
  data <- pkg$data
  add <- ifelse(rownames(param) == "angle", 2, 0)
  multiplier <- pi / prod(param$estimate + add)
  
  tr_sample <- replicate(reps, sampled_traprate())
  tr <- data.frame(parameter="traprate", estimate=traprate(data), se=sd(tr_sample))
  density <- multiplier * tr$estimate
  Es <- c(tr$estimate, param$estimate + add)
  SEs <- c(tr$se, param$se)
  SE <- density * sqrt(sum((SEs/Es)^2))
  res <- data.frame(estimate = c(traprate(data), density),
                    se = c(sd(tr_sample), SE),
                    unit=c(paste0("n/", pkg$data$effort_unit[1]), 
                           paste0("n/", param["radius", "unit"], "2")))
  rownames(res) <- c("traprate", "density")
  res
}


#' Integrated random encounter model density estimate
#' 
#' Estimates animal density for a given species given a camtrap DP datapackage.
#' Models for detection radius and angle, speed and activity level can be 
#' fitted externally and provided as arguments, or are fitted internally if not 
#' provided. Input units are assumed to be distance in m and time in seconds.
#' 
#' @param package Camera trap data package object, as returned by
#'   `read_camtrap_dp()`.
#' @param check_deployments Logical indicating whether to check deployment
#' calibration model diagnostic plots. If `TRUE` (default) runs 
#' `check_deployment_models`; radius, angle and speed data from any excluded 
#' deployments are then dropped from analysis. If `FALSE` all data are used.
#' @param activity_model An activity model fitted using `activity::fitact` or
#' `fit_actmodel`; fitted internally if `NULL`.
#' @param radius_model A detection function model for radii fitted using 
#' `fitdf` or `fit_detmodel` with argument `transect="point"`; fitted 
#' internally if `NULL`.
#' @param angle_model A detection function model for angles fitted using 
#' `fitdf` or `fit_detmodel`; fitted internally if `NULL`.
#' @param speed_model A named vector with elements `estimate` and `se`
#' (giving mean and standard error of speed), as derived from `fit_speedmodel`; 
#' fitted internally if `NULL`.
#' @param species A character string indicating species subset to analyse; user
#'   input required if NULL.
#' @param reps Number of bootstrap replicates for error estimation. 
#' @return A dataframe with .
#' @seealso \code{\link{Distance::ds}}
#' @family density estimation functions
#' @export
rem_estimate <- function(package,
                         check_deployments=TRUE,
                         radius_model=NULL,
                         angle_model=NULL,
                         speed_model=NULL,
                         activity_model=NULL,
                         species=NULL,
                         reps=999){
  
  if(check_deployments) package <- check_deployment_models(package)
  if(is.null(species)) species <- select_species(package)
  
  if(is.null(activity_model)) 
    activity_model <- fit_actmodel(package, species, reps)
  
  if(is.null(radius_model)) 
    radius_model <- fit_detmodel(radius~1, package, species,
                                 order=0, truncation=10)
  
  if(is.null(angle_model))
    angle_model <- fit_detmodel(angle~1, package, species, 
                                order=0, unit="radian")
  
  if(is.null(speed_model))
    speed_model <- fit_speedmodel(package, species)
  
  data <- get_rem_data(package, species, unit="day")
  param <- get_parameter_table(radius_model, angle_model, 
                               speed_model, activity_model)
  res <- rbind(param, rem(data, param))
  res$cv <- 100 * res$se / res$estimate
  res$n <- c(nrow(radius_model$data),
             nrow(angle_model$data),
             length(speed_model$data),
             length(activity_model@data),
             nrow(data),
             NA)
  list(species=species, data=data, estimates=res,
       speed_model=speed_model, activity_model=activity_model, 
       radius_model=radius_model, angle_model=angle_model)
}

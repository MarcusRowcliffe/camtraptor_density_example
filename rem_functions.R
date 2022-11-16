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
fit_speedmodel <- function(package, species=NULL){
  if(is.null(species)) species <- select_species(package)
  obs <- package$data$observations %>%
    subset(scientificName==species & speed > 0.01 & speed < 10)
  if("useDeployment" %in% names(obs)) obs <- subset(obs, useDeployment)
  mn <- 1/mean(1/obs$speed, na.rm=FALSE)
  se <- mn^2 * sqrt(var(1/obs$speed, na.rm=FALSE)/nrow(obs))
  list(speed=data.frame(estimate=mn, se=se), data=obs$speed)
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
#' @param reps Number of bootstrap replicates to run.
#' @return An `actmod` list.
#' @seealso \code{\link{activity::fitact}}
#' @family density estimation functions
#' @export
fit_actmodel <- function(package, 
                         species=NULL, 
                         reps=999){
  if(is.null(species)) species <- select_species(package)
  deps <- package$data$deployments
  obs <- package$data$observations
  if(sum(obs$scientificName==species, na.rm=TRUE)>1){ 
    obs <- left_join(obs, dplyr::select(deps, deploymentID, latitude, longitude),
                     by="deploymentID")
    suntimes <- insol::daylength(obs$latitude, obs$longitude, 
                                 insol::JD(obs$timestamp), 0)
    timeshift <- pi - mean(suntimes[, 1] + suntimes[,3]/2) * pi/12
    obs$solartime <- obs %>%
      with(activity::solartime(timestamp, latitude, longitude, 0)) %>%
      .$solar %>%
      + timeshift %>%
      activity::wrap()
    obs %>%
      subset(scientificName==species) %>%
      .$solartime %>%
      activity::fitact(adj = 1.5, sample = "data", reps = reps)
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
#' @return A `ddf` detection function model list, with additional element
#'   `edd`, a vector with estimated and standard error effective detection 
#'   distance, or the `newdata` dataframe with EDD estimate and se added.
#' @seealso \code{\link{Distance::ds}}
#' @family density estimation functions
#' @export
fit_detmodel <- function(formula, 
                         package, 
                         species=NULL, 
                         newdata=NULL, ...){
  
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
  args <- c(data=list(data), formula=formula[-2], list(...))
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
  mod
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
#' @return A `ddf` detection function model list, with additional element
#'   `edd`, a vector with estimated and standard error effective detection 
#'   distance, or the `newdata` dataframe with EDD estimate and se added.
#' @seealso \code{\link{Distance::ds}}
#' @family density estimation functions
#' @export
fit_detmodel <- function(formula, 
                         package, 
                         species=NULL, 
                         newdata=NULL, ...){
  
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
  args <- c(data=list(data), formula=formula[-2], list(...))
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
  mod
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
                         activity_model=NULL,
                         radius_model=NULL,
                         angle_model=NULL,
                         speed_model=NULL,
                         species=NULL,
                         reps=999){
  
  if(check_deployments) package <- check_deployment_models(package)
  if(is.null(species)) species <- select_species(package)
  
  if(is.null(activity_model)) 
    activity_model <- fit_actmodel(package, species, reps)
  
  if(is.null(radius_model)) 
    radius_model <- fit_detmodel(radius~1, package, species,
                                 transect="point", order=0, truncation=10)
  
  if(is.null(angle_model))
    angle_model <- fit_detmodel(angle~1, package, species, order=0)
  
  if(is.null(speed_model))
    speed_model <- fit_speedmodel(package, species)
  
  data <- data.frame(
    observations = get_n_individuals(package, species=species)$n,
    effort = get_effort(package, unit="second")$effort) %>%
    suppressMessages()
  param <- data.frame(parameter=c("radius", "angle", "speed", "activity"),
                      rbind(radius_model$edd, 
                            angle_model$edd * 2,
                            speed_model$speed, 
                            activity_model@act[1:2]))
  rownames(param) <- NULL
  res <- rem(data, param) %>%
    dplyr::mutate(estimate = estimate * c(1, 180/pi, 1, 1, 86400, 1e6),
                  se = se * c(1, 180/pi, 1, 1, 86400, 1e6))
  res$CV <- 100 * res$se / res$estimate
  res$n <- c(nrow(radius_model$data),
             nrow(angle_model$data),
             length(speed_model$data),
             length(activity_model@data),
             nrow(data),
             NA)
  res$unit = c("m", "deg", "m/s", "none", "n/d", "n/km2")
  list(species=species, data=data, estimates=res,
       speed_model=speed_model, activity_model=activity_model, 
       radius_model=radius_model, angle_model=angle_model)
}

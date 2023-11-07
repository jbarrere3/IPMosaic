#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien Barrere
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to get the path of a file, and create directories if they don't exist
#' @param file.in character: path of the file, filename included (ex: "plot/plot.png")
create_dir_if_needed <- function(file.in){

  path.in <- strsplit(file.in, "/")[[1]]
  if(length(path.in) > 1){
    for(i in 1:(length(path.in)-1)){
      if(i == 1) path.in_i <- path.in[i]
      else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
      if(!dir.exists(path.in_i)) dir.create(path.in_i)
    }
  }
}


#' Function to load demographic parameters for several species
#' @param species.names character vector of species ("Genus_species" format)
load_param_demo <- function(species.in){

  # Initialize list of fits (one element per species)
  fit.list <- list()

  # Loop on all species to load and store demographic parameters
  for(i in 1:length(species.in)){
    eval(parse(text=paste0("fit.list$", species.in[i], " <- fit_", species.in[i])))
  }

  # Return the list
  return(fit.list)

}



#' Function to create a list of climates. So far, just take optimum of sp in input
#' @param species.in Species for which to get the optimum
make_climate_list = function(species.in){

  list(clim1 = subset(
    matreex::climate_species, N == 2 & sp == species.in, select = -sp))

}

#' Function to create a list of IPM to run
#' @param climate.list list of climates
#' @param species.in vector of species to include in all climates
make_species_list = function(climate.list, species.in){

  expand.grid(species = species.in, climate = names(climate.list)) %>%
    mutate(file = paste0("rds/", climate, "/species/", species, ".rds")) %>%
    mutate(ID.species = c(1:dim(.)[1]),
           climate = as.character(climate)) %>%
    dplyr::select(ID.species, climate, species, file)


}

#' Function to make a species object, save it as rds and return filename
#' @param species_list table containing all species to create
#' @param climate.list list of climates object
#' @param param.demo demographic parameter for all species
#' @param ID.species.in ID of the species to make in species.list
make_species_rds = function(
  species.list, climate.list, param.demo, ID.species.in){

  # Make IPM
  IPM.in = make_IPM(
    species = as.character(species.list$species[ID.species.in]),
    climate = climate.list[[species.list$climate[ID.species.in]]],
    fit =  param.demo[[species.list$species[ID.species.in]]],
    clim_lab = species.list$climate[ID.species.in],
    mesh = c(m = 700, L = 100, U = as.numeric(
      param.demo[[species.list$species[ID.species.in]]]$info[["max_dbh"]]) * 1.1),
    BA = 0:200, verbose = TRUE, correction = "none"
  )

  dist_coef <- matreex::disturb_coef %>%
      dplyr::filter(species ==  as.character(species.list$species[ID.species.in]))

  # Create species object
  species.out = species(
    IPM.in, init_pop = def_initBA(20), harvest_fun = def_harv,
    disturb_fun = def_disturb, disturb_coef = dist_coef)

  # Save species object in a rdata
  create_dir_if_needed(species.list$file[ID.species.in])
  saveRDS(species.out, species.list$file[ID.species.in])

  # Return the file created
  return(species.list$file[ID.species.in])

}


#' Disturbance function
#'
#' @param x population state distribution at time t
#' @param species The species class object of interest to get mesh and RDIcoef
#' values from. RDIcoef is a one line dataframe with RDI coefficient for one
#' species.
#' @param disturb Disturbance parameters. Highly depend on the disturbance
#' impact parameters given to the species.
#' @param ... Not used in this case.
#' \describe{
#' \item{qmd}{Forest Quadratic Mean Diameter}
#' }
#' @author Maxime Jeaunatre
#'
disturb_fun <- function(x, species, disturb = NULL, ...){

  dots <- list(...)
  qmd <- dots$qmd
  size <- species$IPM$mesh
  coef <- species$disturb_coef
  if(any(disturb$type %in% coef$disturbance)){
    coef <- subset(coef, disturbance == disturb$type)
  } else {
    stop(sprintf("The species %s miss this disturbance type (%s) parameters",
                 sp_name(species), disturb$type))
  }

  # edits for delay
  size[size == 0] <- min(size[size !=0])

  logratio <-  log(size / qmd)
  dbh.scaled = coef$dbh.intercept + size * coef$dbh.slope
  logratio.scaled = coef$logratio.intercept + logratio * coef$logratio.slope
  Pkill <- plogis(coef$a0 + coef$a1 * logratio.scaled +
                    coef$b * disturb$intensity ^(coef$c * dbh.scaled))

  return(x* Pkill) # always return the mortality distribution
}



#' generate a dataframe listing all forests to simulate
#' @param species.list dataframe listing all species object to make
make_forest_list = function(species.list){

  # Loop on all climates
  for(j in 1:length(unique(species.list$climate))){
    # Name of the climate
    clim.j = unique(species.list$climate)[j]
    # Identify the list of species for this climate
    species.in.j = unique(subset(species.list, climate == clim.j)$species)
    # Initialize the vector of species combinations for climate j
    combi.j = c()
    # Loop on all level of richness
    for(r in 1:length(species.in.j)){
      # Loop on all combinations in the level of richness r
      for(i in 1:dim(t(combn(species.in.j, r)))[1]){
        combi.j = c(combi.j, paste(t(combn(species.in.j, r))[i, ], collapse = "."))
      }
    }
    # Make output data.frame for climate j
    out.j = data.frame(climate = clim.j, combi = combi.j)
    # Add to final output
    if(j == 1) out = out.j
    else out = rbind(out, out.j)
  }

  # Final formatting
  out = out %>%
    mutate(ID.forest = c(1:dim(.)[1]),
           climate = as.character(climate)) %>%
    dplyr::select(ID.forest, climate, combi)

  # return output
  return(out)
}



#' Function to make a list of simulations till equilibrium
#' @param species_list df with information on all species object
#' @param forest_list df with information on all forest to simulate
#' @param species vector containing all species rds files created
#' @param ID.forest.in ID of the forest to simulate in forest_list
make_simulations_equilibrium = function(
  species.list, forest.list, species, ID.forest.in){

  # Harvesting rules
  harv_rules.ref = c(Pmax = 0.25, dBAmin = 3, freq = 5, alpha = 1)

  # Identify ID of the climate
  climate.in = forest.list$climate[ID.forest.in]

  # Identify species combination i
  combination.in = forest.list$combi[ID.forest.in]

  # vector of species in forest i
  species.in = unlist(strsplit(combination.in, "\\."))

  # Initialize list of species to create
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in

  # Loop on all species
  for(i in 1:length(species.in)){

    # Identify the file in species containing species i
    species.file.i = species[(species.list %>%
                                filter(species == species.in[i]) %>%
                                filter(climate == climate.in))$ID.species]

    # Store the file in the list
    list.species[[i]] = readRDS(species.file.i)

  }


  # Make forest
  forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)

  # Run simulation till equilibrium
  sim.in = sim_deter_forest(
    forest.in, tlim = 5000, equil_time = 5000, equil_dist = 1000,
    equil_diff = 0.5, harvest = "default", SurfEch = 0.03, verbose = TRUE)

  # Save simulation in a rdata
  file.out = paste0("rds/", climate.in, "/sim_equil/", combination.in, ".rds")
  create_dir_if_needed(file.out)
  saveRDS(sim.in, file.out)

  # Return output list
  return(file.out)
}


#' Function to make a list of simulations till equilibrium
#' @param species_list df with information on all species object
#' @param forest_list df with information on all forest to simulate
#' @param species vector containing all species rds files created
#' @param dist boolean to indicate whether we should include disturbance regimes or not
#' @param ID.forest.in ID of the forest to simulate in forest_list
make_simulations_reg = function(
  species.list, forest.list, species, sim_equilibrium, dist, ID.forest.in){

  # Harvesting rules
  harv_rules.ref = c(Pmax = 0.25, dBAmin = 3, freq = 5, alpha = 1)

  # Identify ID of the climate
  climate.in = forest.list$climate[ID.forest.in]

  # Identify species combination i
  combination.in = forest.list$combi[ID.forest.in]

  # vector of species with non null abundance in forest i
  species.in = unlist(strsplit(combination.in, "\\."))

  # vector of all species in forest i (including regional pool)
  species.in.all = as.character(unique(subset(species.list, climate == climate.in)$species))

  # Vector of migration rate
  migration_rate.in = rep(0.1, length(species.in.all))
  names(migration_rate.in) = species.in.all

  # Get relative abundance at equilibrium
  # -- Identify forest with all species together
  ID.forest.equil = (forest.list %>%
                       filter(climate == climate.in) %>%
                       filter(combi == paste(
                         species.in.all, collapse = ".")))$ID.forest
  # -- read corresponding simulation at equil
  sim.equil = readRDS(sim_equilibrium[ID.forest.equil])
  # -- Extract basal area at equilibrium per species
  equil.BA = dplyr::filter(sim.equil, equil, var == "BAsp")$value
  names(equil.BA) = dplyr::filter(sim.equil, equil, var == "BAsp")$species
  # -- Extract size distribution at equilibrium
  equil_dist <- dplyr::filter(sim.equil, equil, var == "n") %>%
    dplyr::group_by(species) %>%
    dplyr::group_split() %>% map(pull, value)

  # Initialize list of species to create
  list.species <- vector("list", length(species.in.all))
  names(list.species) = species.in.all

  # Loop on all species
  for(i in 1:length(species.in.all)){

    # Identify the file in species containing species i
    species.file.i = species[(species.list %>%
                                filter(species == species.in.all[i]) %>%
                                filter(climate == climate.in))$ID.species]

    # Read species
    species.i = readRDS(species.file.i)

    # If species is not in the combination, set its abundance to 0
    if(!(species.in.all[i] %in% species.in)){
      species.i$init_pop = def_init_k(equil_dist[[i]]*0)}

    # Update disturbance function if disturbances are included
    if(dist) species.i$disturb_fun = disturb_fun

    # Store the file in the list
    list.species[[i]] = species.i

  }


  # Make forest
  forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref,
                         regional_abundance = equil.BA,
                         migration_rate = migration_rate.in)

  # Run simulation
  # -- Time of simulation
  tsim = 5000
  # -- In the absence of disturbances
  if(!dist){
    sim.in = sim_deter_forest(
      forest.in, tlim = tsim, equil_time = tsim, equil_dist = 1000,
      equil_diff = 0.5, harvest = "default", SurfEch = 0.03, verbose = TRUE)}
  # -- In the absence of disturbances
  if(dist){
    sim.in = sim_deter_forest(
      forest.in, tlim = tsim, equil_time = tsim, equil_dist = 1000,
      equil_diff = 0.5, harvest = "default", SurfEch = 0.03, verbose = TRUE,
      disturbance = make_disturbance.df(t = c(1:tsim), freq = 0.01))}


  # Save simulation in a rdata
  if(dist) file.out = paste0("rds/", climate.in, "/sim_reg_dist/", combination.in, ".rds")
  else file.out = paste0("rds/", climate.in, "/sim_reg/", combination.in, ".rds")
  create_dir_if_needed(file.out)
  saveRDS(sim.in, file.out)

  # Return output list
  return(file.out)
}



#' Function to create disturbance dataframe
#' @param type type of disturbance ("storm", "fire", or "biotic")
#' @param freq frequency of the disturbance
#' @param I.param parameter of beta distribution for intensity
#' @param t time vector
make_disturbance.df = function(
  type = "storm", freq = 0.01, I.param = c(1.97, 5.82), t = c(1:5000)){

  data.frame(
    type = type,
    intensity = rbeta(length(t), shape1 = I.param[1], shape2 = I.param[2]),
    IsSurv = FALSE, t = t)[which(rbinom(length(t), 1, freq) == 1), ]

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## matreex functions for regional forest class -------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' species are initiated with their own functions and only the regional abundance
#' is set in the forest object as well as the migration rate

#' Constructor for forest class
#'
#' Only used in the matreex package
#'
#' @param species List of species created with matreex package.
#' @param harv_rules Vector for harvest rules at the scale of the forest.
#' \describe{
#'   \item{Pmax}{maximum proportion of BAcut / BA}
#'   \item{dBAmin}{the minimum BA to perform cut}
#'   \item{freq}{Frequence at which the harvest will be executed.}
#' }
#' @param regional_abundance list of vector with size distribution for each species.
#' This is not a direct species relative abundance, but I don't know how to implement this...help ,
#' @param migration_rate numeric vector with a migration rate in percentage between 1 et and 0.
#'
#' @importFrom purrr map_chr
#'
#' @keywords internal
#' @export
new_forest <- function(species = list(),
                       harv_rules = c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1),
                       regional_abundance = NULL,
                       migration_rate = NULL
){

  sp <- map_chr(species, sp_name)
  names(species) <- sp
  if(!is_null(regional_abundance)){
    names(regional_abundance) <- sp
    names(migration_rate) <- sp
  }
  forest <- list(
    species = species, harv_rules = harv_rules,
    info = list(species = sp,
                clim_lab = map_chr(species, climatic)),
    regional_abundance = regional_abundance,
    migration_rate = migration_rate
  )

  if(!is_null(regional_abundance)){
    class(forest) <- c("forest", "reg_forest")
  } else {
    class(forest) <- "forest"
  }

  return(forest)
}

#' validator for forest class.
#'
#' @param x forest class object
#'
#' @import checkmate
#'
#' @noRd
validate_forest <- function(x){

  regional <- inherits(x, "reg_forest")
  values <- unclass(x)
  names <- attr(x, "names")

  #map(values$species, validate_species)
  # TODO check forest harv rules

  clim_lab <- values$info$clim_lab
  if(length(unique(clim_lab)) > 1){
    clim_ipm <- clim_lab[clim_lab != "mu_gr"]
    if(length(clim_ipm) > 1){ # D & F
      stop(paste0("Some ipm species are not defined with the same climatic name.",
                  "Check it with : map_chr(species, climatic)"))
    }
  }

  # check the regional pool settings
  if(regional){

    assertNumeric(values$migration_rate, len = length(values$species),
                  lower = 0, upper = 1)
    if(all(values$migration_rate == 0)){
      warning("All migration rate are 0, the regional pool of this forest is deleted")
      x$regional_abundance <- NULL
      x$migration_rate <- NULL
      class(x) <- "forest"

      return(invisible(x))
    }

    assertSubset(names(values$migration_rate), names(values$species))
    # length_meshs <- map_dbl(values$species, ~ length(.x$IPM$mesh))

    # assertList(values$regional_abundance, types = "numeric",
    # len = length(values$species))
    # if(any(lengths(values$regional_abundance) != length_meshs)){
    # stop("regional abundance numeric vector should be the length of the species mesh.")
    # }

    assertSubset(names(values$regional_abundance), names(values$species))
  }

  invisible(x)
}

#' Create a new forest for simulation
#'
#' A forest is a group of one of multiple species to silumate along time using
#' the IPM defined for each species and harvest rules.
#'
#' @inheritParams new_forest
#'
#' @export
forest <- function(species = list(),
                   harv_rules = c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1),
                   regional_abundance = NULL,
                   migration_rate = NULL
){

  res <- validate_forest(new_forest(
    species = species,
    harv_rules = harv_rules,
    regional_abundance = regional_abundance,
    migration_rate = migration_rate
  ))

  return(res)
}



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
    mutate(ID.species = c(1:dim(.)[1])) %>%
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
    clim_lab = paste0("IDspecies_", ID.species.in),
    mesh = c(m = 700, L = 100, U = as.numeric(
      param.demo[[species.list$species[ID.species.in]]]$info[["max_dbh"]]) * 1.1),
    BA = 0:200, verbose = TRUE, correction = "none"
  )
  
  # Create species object 
  species.out = species(
    IPM.in, init_pop = def_initBA(20), harvest_fun = def_harv, disturb_fun = def_disturb)
  
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
    forest.in, tlim = 2000, equil_time = 2000, equil_dist = 1000, 
    equil_diff = 0.5, harvest = "default", SurfEch = 0.03, verbose = TRUE)
  
  # Save simulation in a rdata
  file.out = paste0("rds/", climate.in, "/sim_equil/", combination.in, ".rds")
  create_dir_if_needed(file.out)
  saveRDS(sim.in, file.out)
  
  # Return output list
  return(file.out)
}




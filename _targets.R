#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name _targets.R  
#' @description R script to launch the target pipeline
#' @author Julien BARRERE
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options and packages ----------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load targets
library(targets)
# Load functions
lapply(grep("R$", list.files("R"), value = TRUE), 
       function(x) source(file.path("R", x)))
# install if needed and load packages
packages.in <- c("dplyr", "ggplot2", "matreex", "tidyr", "data.table", 
                 "rnaturalearth", "cowplot", "future")
for(i in 1:length(packages.in)){
  if(!(packages.in[i] %in% rownames(installed.packages()))){
    install.packages(packages.in[i])
  } 
} 
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient")
future::plan(future::multisession, workers = 6)
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  
  # Reference climate for the first simulations (silver fir optimum)
  tar_target(climate.list, make_climate_list("Abies_alba")),

  # Vector of all species to include in the simulations
  tar_target(species.in, c("Abies_alba", "Fagus_sylvatica", "Picea_abies")),

  # Get demographic parameters of the species included in the simulations
  tar_target(param.demo, load_param_demo(species.in)),

  # Make a list of species object to create
  tar_target(species.list, make_species_list(climate.list, species.in)), 
  
  # Make species objects
  # -- First: Make a vector of ID for each species to make
  tar_target(ID.species, species.list$ID.species), 
  # -- Make the species via branching over ID.species
  tar_target(species, make_species_rds(
    species.list, climate.list, param.demo, ID.species.in = ID.species), 
    pattern = map(ID.species), iteration = "vector", format = "file"), 
  
  # Prepare forests
  # -- Make the forest_list dataframe from species.list
  tar_target(forest.list, make_forest_list(species.list)), 
  # -- Make a vector of ID for each forest to simulate
  tar_target(ID.forest, forest.list$ID.forest),
  
  # Make simulations till equilibrium
  tar_target(sim_equilibrium, make_simulations_equilibrium(
    species.list, forest.list, species, ID.forest), 
    pattern = map(ID.forest), iteration = "vector", format = "file"), 
  
  # Plot the simulations
  tar_target(fig_simequil, plot_sim(
    sim_equilibrium, "output/fig/fig_sim_equil.jpg"), format = "file")
  
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_plot.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien Barrere
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Plot change over time in basal area for a list of simulations saved as rds
#' @param sim.in vector of character containing path to simulations saved as rds
#' @param file.out Name including path of the file to save
plot_sim <- function(sim.in, file.out){
  
  # Create directory if needed
  create_dir_if_needed(file.out)
  
  # Loop on all species combinations
  for(i in 1:length(sim.in)){
    
    # Printer
    print(paste0(i, "/", length(sim.in)))
    
    # Read simulation
    sim.i = readRDS(sim.in[i])
    
    # Vector of species for simulation i
    sp.vec.i <- unique(sim.i$species)
    
    # Format data for plotting simulation i
    data.i <- sim.i %>%
      # Remove unused lines and columns
      filter(var == "BAsp") %>%
      dplyr::select(species, time, value) %>%
      distinct() %>%
      # Calculate the sum of basal area
      spread(key = "species", value = "value") %>%
      mutate(all = rowSums(.[, c(2:dim(.)[2])])) %>%
      gather(key = "species", value = "value", c(sp.vec.i, "all")) %>%
      # Add simulation ID
      mutate(ID = i)
    
    # Add to final dataset
    if(i == 1) data = data.i
    else data <- rbind(data, data.i)
  }
  
  # Make the plot 
  plot.out <- data %>%
    mutate(species = factor(
      species, levels = c("all", unique((filter(data, species != "all"))$species)))) %>%
    ggplot(aes(x = time, y = value, group = species, color = species)) + 
    geom_line() + 
    facet_wrap(~ ID) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank(),
          legend.title = element_blank(), 
          legend.key = element_blank()) + 
    xlab("Time (years)") + ylab("Basal area")
  
  # - Save the plot
  ggsave(file.out, plot.out, width = 21, height = 14, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.out)
}
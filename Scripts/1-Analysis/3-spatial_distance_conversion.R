# This script should map scaled climatic differences to spatial distance

setwd(project_directory)

for(i in 1:length(results_list)){
  spc_name <- names(results_list)[i]
  spec_i <- results_list[[spc_name]]
  raw.climate_B <- spec_i[["pca_B"]][["pca_source"]]
  temp.dist_B <- mean(raw.climate_B[raw.climate_B$ID == "Novel", "MWMT"]) - 
    mean(raw.climate_B[raw.climate_B$ID == "Both", "MWMT"])
  
  raw.climate_N <- spec_i[["pca_N"]][["pca_source"]]
  temp.dist_N <- mean(raw.climate_N[raw.climate_N$ID == "Novel", "MWMT"]) - 
    mean(raw.climate_N[raw.climate_N$ID == "Both", "MWMT"])
  
  lapse.rate <- 6.5 # ºC/km
  # ISO International Standard 2533–1975. 
  # “Standard Atmosphere First Edition,” 
  # Corrigendum 1, 1978. ISO, Geneva, Switzerland.
  
  alt.dist_B <- -1 * temp.dist_B / lapse.rate * 1000 # -1 because cooler temperatures are uphill
  alt.dist_N <- -1 * temp.dist_N / lapse.rate * 1000
  
  results_list[[spc_name]][["alt.dist_B"]] <- alt.dist_B
  results_list[[spc_name]][["alt.dist_N"]] <- alt.dist_N
}



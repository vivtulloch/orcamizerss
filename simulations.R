# Salish Sea Ecosystem Model using mizer
# Marine ecosystem modeling with focus on orcas, salmon, and marine mammals
# Author: Viv Tulloch
# Date: June 12 2025

# Load required libraries
remotes::install_github("sizespectrum/mizerExperimental")
library(mizerExperimental)
remotes::install_github("sizespectrum/mizerMR")
library(mizerMR)
library(mizer)

library(ggplot2)
library(plotly)
library(dplyr)
library(reshape2)
library(abind)
library(tidyverse)


# =============================================================================
# SECTION 1: MODEL SETUP AND PARAMETER TUNING
# =============================================================================

# 1.1 Baseline Model Assessment
# Explore unexploited community state (zero fishing effort)
assess_unexploited_state <- function(model, t_max = 100) {
  sim_unexploited <- project(model, 
                             effort = 0,  
                             t_max = t_max,
                             return_sim = TRUE)
  
  plotlyBiomassRelative(sim_unexploited)
  return(sim_unexploited)
}

# 1.2 Save Model Parameters
save_model_parameters <- function(model, filename_base, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save RDS file
  rds_path <- file.path(output_dir, paste0(filename_base, ".rds"))
  saveParams(model, rds_path)
  
  # Save species parameters as CSV
  species_csv_path <- file.path(output_dir, paste0(filename_base, "_species.csv"))
  write.csv(species_params(model), species_csv_path, row.names = FALSE)
  
  # Save gear parameters as CSV
  gear_csv_path <- file.path(output_dir, paste0(filename_base, "_gears.csv"))
  write.csv(gear_params(model), gear_csv_path, row.names = FALSE)
  
  cat("Model parameters saved to:", output_dir, "\n")
}

# 1.3 Tune Marine Mammal Parameters
tune_marine_mammals <- function(base_model) {
  # Create copy of base model
  tuned_model <- base_model
  
  # Adjust marine mammal reproduction efficiency and maximum intake
  # Indices 19-21: Marine mammal species (adjust based on your species list)
  species_params(tuned_model)$erepro[19:21] <- 0.008
  species_params(tuned_model)$h[19:21] <- species_params(tuned_model)$h[19:21] * 0.5
  
  # Adjust predation parameters for specific species
  # Index 22: Large predator (adjust based on your species list)
  species_params(tuned_model)$gamma[22] <- species_params(tuned_model)$gamma[22] * 3
  species_params(tuned_model)$h[22] <- species_params(tuned_model)$h[22] * 1.5
  
  # Index 4: Another key species
  species_params(tuned_model)$erepro[4] <- 0.7
  species_params(tuned_model)$gamma[4] <- species_params(tuned_model)$gamma[4] * 10
  species_params(tuned_model)$h[4] <- species_params(tuned_model)$h[4] * 1.2
  
  # Re-calibrate model
  tuned_model <- setParams(tuned_model)
  tuned_model <- setFishing(tuned_model, effort = 1)
  tuned_model <- tuned_model |> steady() |> matchBiomasses() |> steady()
  
  return(tuned_model)
}

# 1.4 Fine-tune Ecosystem Parameters
fine_tune_ecosystem <- function(base_model) {
  tuned_model <- base_model
  
  # Adjust feeding and metabolism parameters for key species
  # Marine mammals
  species_params(tuned_model)$h[19] <- species_params(tuned_model)$h[19] * 1
  species_params(tuned_model)$gamma[19] <- species_params(tuned_model)$gamma[19] * 0.3
  
  species_params(tuned_model)$h[20:21] <- species_params(tuned_model)$h[20:21] * 1.4
  species_params(tuned_model)$gamma[20:21] <- species_params(tuned_model)$gamma[20:21] * 1.3
  
  # Large predators
  species_params(tuned_model)$gamma[22] <- species_params(tuned_model)$gamma[22] * 0.5
  species_params(tuned_model)$gamma[4] <- species_params(tuned_model)$gamma[4] / 4
  species_params(tuned_model)$h[4] <- species_params(tuned_model)$h[4] * 1.2
  
  # Mid-level predators
  species_params(tuned_model)$gamma[14] <- species_params(tuned_model)$gamma[14] * 1.4
  species_params(tuned_model)$h[14] <- species_params(tuned_model)$h[14] * 1.5
  
  # Adjust resource capacity
  current_capacity <- resource_capacity(tuned_model)
  tuned_model <- setResource(tuned_model, resource_capacity = current_capacity * 1.2)
  
  # Re-calibrate
  tuned_model <- setParams(tuned_model)
  tuned_model <- setFishing(tuned_model, effort = 1)
  tuned_model <- tuned_model |> steady() |> matchBiomasses() |> steady()
  
  return(tuned_model)
}

# =============================================================================
# SECTION 2: HISTORICAL SIMULATIONS
# =============================================================================

# 2.1 Run Historical Simulation
run_historical_simulation <- function(model, effort_matrix, dt = 0.1) {
  # Project to steady state without fishing first
  model_no_fishing <- projectToSteady(model, t_max = 30, effort = 0)
  
  # Enhance resource capacity and adjust marine mammal parameters
  current_capacity <- resource_capacity(model_no_fishing)
  model_enhanced <- setResource(model_no_fishing, 
                                resource_capacity = current_capacity * 2)
  
  # Adjust marine mammal mortality and reproduction
  species_params(model_enhanced)$z0[19:21] <- species_params(model_enhanced)$z0[19:21] * 0.9
  species_params(model_enhanced)$z0[20:21] <- species_params(model_enhanced)$z0[20:21] * 0.2
  species_params(model_enhanced)$erepro[19:21] <- species_params(model_enhanced)$erepro[19:21] * 1.5
  
  # Adjust pinniped parameters
  species_params(model_enhanced)$erepro[24:25] <- species_params(model_enhanced)$erepro[24:25] * 2
  species_params(model_enhanced)$z0[24:25] <- species_params(model_enhanced)$z0[24:25] * 0.2
  
  model_enhanced <- setParams(model_enhanced)
  
  # Run historical simulation
  sim_historical <- project(model_enhanced,
                            effort = effort_matrix, 
                            dt = dt, 
                            t_save = 1)
  
  return(list(model = model_enhanced, simulation = sim_historical))
}

# 2.2 Compare with Survey Data
compare_with_surveys <- function(simulation, survey_file_path) {
  # Load survey data
  survey_data <- read.csv(survey_file_path)
  
  # Define marine mammal species
  marine_mammal_species <- c("Harbor_seals", "Steller_sealions", "California_sealions", 
                             "Resident_Orca_S", "Resident_Orca_N", "Transient_Orca")
  
  # Extract biomass predictions
  biomass_matrix <- getN(simulation)
  biomass_selected <- biomass_matrix[, marine_mammal_species]
  
  # Calculate relative biomass (relative to 2011 baseline)
  initial_biomass <- biomass_selected["2011", ]
  relative_biomass <- sweep(biomass_selected, 2, initial_biomass, "/")
  
  # Convert to long format for plotting
  relative_biomass_long <- melt(relative_biomass, 
                                varnames = c("time", "sp"), 
                                value.name = "relative_biomass")
  
  # Rename species for better plotting
  species_name_mapping <- c(
    "Harbor_seals" = "Harbor Seals",
    "Steller_sealions" = "Steller Sea Lions",
    "California_sealions" = "California Sea Lions",
    "Resident_Orca_S" = "SRKW",
    "Resident_Orca_N" = "NRKW",
    "Transient_Orca" = "TKW"
  )
  
  relative_biomass_long <- relative_biomass_long %>%
    mutate(sp = recode(sp, !!!species_name_mapping))
  
  # Prepare survey data
  survey_long <- melt(survey_data, id.vars = "time", 
                      variable.name = "sp", value.name = "value")
  survey_long <- survey_long %>%
    mutate(sp = recode(sp, !!!species_name_mapping))
  
  # Create comparison plot
  comparison_plot <- ggplot(relative_biomass_long) +
    geom_line(aes(x = time, y = relative_biomass), size = 1) +
    geom_point(data = survey_long, aes(x = time, y = value), 
               color = "blue", size = 1.5) +
    geom_smooth(data = survey_long, aes(x = time, y = value),
                method = "lm", color = "blue", linetype = "dashed", 
                size = 0.5, se = FALSE) +
    facet_wrap(~sp, scales = "free", nrow = 2) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Year", y = "Relative Biomass",
         title = "Model Predictions vs Survey Data",
         subtitle = "Blue points and dashed line show survey data trends")
  
  return(list(plot = comparison_plot, 
              predictions = relative_biomass_long,
              survey = survey_long))
}

# =============================================================================
# SECTION 3: FUTURE PROJECTIONS
# =============================================================================

# 3.1 Create Future Effort Scenarios
create_effort_scenarios <- function(historical_effort_matrix) {
  # Extract years and sectors
  historical_years <- rownames(historical_effort_matrix)
  sectors <- colnames(historical_effort_matrix)
  
  # Get 2011 reference values
  reference_2011 <- historical_effort_matrix["2011", ]
  
  # Scenario 1: Status quo (continue 2011 levels)
  scenario1 <- create_status_quo_scenario(historical_effort_matrix, reference_2011)
  
  # Scenario 2: Reduced salmon fishing
  scenario2 <- create_reduced_salmon_scenario(historical_effort_matrix)
  
  return(list(status_quo = scenario1, reduced_salmon = scenario2))
}

create_status_quo_scenario <- function(historical_effort, reference_levels) {
  # Subset historical data to 1986-2011
  effort_base <- historical_effort[as.character(1986:2011), , drop = FALSE]
  
  # Create future projections (2012-2050)
  n_years_future <- length(2012:2050)
  effort_future <- array(rep(reference_levels, each = n_years_future),
                         dim = c(n_years_future, dim(historical_effort)[2]),
                         dimnames = list(as.character(2012:2050), 
                                         dimnames(historical_effort)[[2]]))
  
  # Combine historical and future
  scenario <- abind::abind(effort_base, effort_future, along = 1)
  return(scenario)
}

create_reduced_salmon_scenario <- function(historical_effort) {
  # Base scenario (status quo)
  scenario <- create_status_quo_scenario(historical_effort, historical_effort["2011", ])
  
  # Define fishing effort reductions for different sectors
  effort_reductions <- list(
    # Salmon-targeting sectors (50% reduction)
    salmon_sectors = c("Chinook_Hatch_Adult_PS", "Chinook_Wild_Adult_PS", 
                       "Coho_Hatch_Adult_PS", "Coho_Wild_Adult_PS",
                       "Chum_Adult_PS", "Sockeye_Adult_PS",
                       "SOG_salmon", "SOG_recreational_salmon"),
    salmon_multiplier = 0.5,
    
    # Other commercial sectors (slight increase)
    other_sectors = c("Invertebrates_PS", "Pacific_Herring_PS", "Hake_PS", 
                      "Lingcod_PS", "Rockfish_PS", "SOG_bottom_trawl",
                      "SOG_hook_and_line", "SOG_midwater_trawl"),
    other_multiplier = 1.1,
    
    # Bycatch and incidental mortality (eliminate)
    bycatch_sectors = c("Seal_mort", "Pinniped_bycatch", "Cetacean_bycatch", "Incidental"),
    bycatch_multiplier = 0
  )
  
  # Apply reductions for future years (2012-2050)
  future_years <- as.character(2012:2050)
  
  # Salmon sectors
  for (sector in effort_reductions$salmon_sectors) {
    if (sector %in% colnames(scenario)) {
      scenario[future_years, sector] <- effort_reductions$salmon_multiplier
    }
  }
  
  # Other commercial sectors
  for (sector in effort_reductions$other_sectors) {
    if (sector %in% colnames(scenario)) {
      scenario[future_years, sector] <- effort_reductions$other_multiplier
    }
  }
  
  # Bycatch sectors
  for (sector in effort_reductions$bycatch_sectors) {
    if (sector %in% colnames(scenario)) {
      scenario[future_years, sector] <- effort_reductions$bycatch_multiplier
    }
  }
  
  return(scenario)
}

# 3.2 Run Future Projections
run_future_projections <- function(model, scenarios, dt = 0.25) {
  projections <- list()
  
  for (scenario_name in names(scenarios)) {
    cat("Running scenario:", scenario_name, "\n")
    projections[[scenario_name]] <- project(model, 
                                            effort = scenarios[[scenario_name]], 
                                            dt = dt)
  }
  
  return(projections)
}

# =============================================================================
# SECTION 4: ANALYSIS AND VISUALIZATION
# =============================================================================

# 4.1 Calculate Biodiversity Reference Points
calculate_reference_points <- function(unexploited_simulation) {
  final_time <- dim(getSSB(unexploited_simulation))[1]
  ssb_unexploited <- getSSB(unexploited_simulation)[final_time, ]
  
  reference_points <- list(
    unexploited = ssb_unexploited,
    reference_level = ssb_unexploited * 0.5  # 50% of unexploited as reference
  )
  
  return(reference_points)
}

# 4.2 Compare Scenarios
compare_scenarios <- function(projections, reference_points, species_list) {
  # Prepare data for comparison
  scenario_data <- list()
  
  for (scenario_name in names(projections)) {
    ssb_data <- melt(getSSB(projections[[scenario_name]]))
    scenario_data[[scenario_name]] <- transform(ssb_data, scenario = scenario_name)
  }
  
  # Combine all scenarios
  combined_data <- do.call(rbind, scenario_data)
  
  # Add reference lines
  years <- unique(combined_data$time)
  reference_data <- expand.grid(sp = species_list, time = years)
  reference_data$value <- rep(reference_points$reference_level, each = length(years))
  reference_data$scenario <- "Reference Level"
  
  unexploited_data <- expand.grid(sp = species_list, time = years)
  unexploited_data$value <- rep(reference_points$unexploited, each = length(years))
  unexploited_data$scenario <- "Unexploited"
  
  # Create final dataset
  all_data <- rbind(combined_data, reference_data, unexploited_data)
  
  # Create comparison plot
  colors <- c("status_quo" = "red", "reduced_salmon" = "green", 
              "Reference Level" = "purple", "Unexploited" = "blue")
  
  comparison_plot <- ggplot(all_data) +
    geom_line(aes(x = time, y = value, color = scenario), size = 1) +
    facet_wrap(~sp, scales = "free", nrow = 5) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    labs(x = "Year", y = "Spawning Stock Biomass",
         title = "Future Projections: Species Recovery Under Different Scenarios",
         color = "Scenario") +
    theme(legend.position = "bottom")
  
  return(list(plot = comparison_plot, data = all_data))
}

# 4.3 Export Model Outputs
export_model_outputs <- function(simulation, output_dir, filename_prefix) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Export different metrics
  metrics <- list(
    biomass = getBiomass(simulation),
    abundance = getN(simulation),
    prop_large_fish = getProportionOfLargeFish(simulation),
    ssb = getSSB(simulation),
    predation_mortality = getM2(simulation),
    yield = getYield(simulation)
  )
  
  for (metric_name in names(metrics)) {
    filename <- file.path(output_dir, paste0(filename_prefix, "_", metric_name, ".csv"))
    write.csv(metrics[[metric_name]], filename, row.names = TRUE)
  }
  
  cat("Model outputs exported to:", output_dir, "\n")
}

# =============================================================================
# SECTION 5: MAIN ANALYSIS WORKFLOW
# =============================================================================

# Main function to run complete analysis
run_complete_analysis <- function(base_model, historical_effort, survey_data_path, output_dir) {
  cat("Starting Salish Sea ecosystem analysis...\n")
  
  # Step 1: Model tuning
  cat("Step 1: Tuning model parameters...\n")
  tuned_model1 <- tune_marine_mammals(base_model)
  final_model <- fine_tune_ecosystem(tuned_model1)
  
  # Step 2: Assess unexploited state
  cat("Step 2: Assessing unexploited state...\n")
  unexploited_sim <- assess_unexploited_state(final_model)
  
  # Step 3: Historical simulation
  cat("Step 3: Running historical simulation...\n")
  historical_results <- run_historical_simulation(final_model, historical_effort)
  
  # Step 4: Compare with survey data
  if (!is.null(survey_data_path) && file.exists(survey_data_path)) {
    cat("Step 4: Comparing with survey data...\n")
    survey_comparison <- compare_with_surveys(historical_results$simulation, survey_data_path)
    print(survey_comparison$plot)
  }
  
  # Step 5: Future projections
  cat("Step 5: Creating future scenarios...\n")
  scenarios <- create_effort_scenarios(historical_effort)
  projections <- run_future_projections(historical_results$model, scenarios)
  
  # Step 6: Analysis and comparison
  cat("Step 6: Analyzing results...\n")
  reference_points <- calculate_reference_points(unexploited_sim)
  species_list <- dimnames(getSSB(projections[[1]]))[[2]]
  scenario_comparison <- compare_scenarios(projections, reference_points, species_list)
  print(scenario_comparison$plot)
  
  # Step 7: Export results
  cat("Step 7: Exporting results...\n")
  save_model_parameters(final_model, "final_model", output_dir)
  export_model_outputs(historical_results$simulation, output_dir, "historical")
  
  for (scenario_name in names(projections)) {
    export_model_outputs(projections[[scenario_name]], output_dir, 
                         paste0("projection_", scenario_name))
  }
  
  cat("Analysis complete! Results saved to:", output_dir, "\n")
  
  return(list(
    model = final_model,
    historical_simulation = historical_results$simulation,
    projections = projections,
    reference_points = reference_points,
    scenario_comparison = scenario_comparison
  ))
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# 
# # Load your base fully paramaterized model and effort/fishing mortality data
# base_model <- readRDS("path/to/your/base_model.rds")
# historical_effort <- read.csv("path/to/effort_data.csv")  # or load from RDS


#All_fisheries_effort_array <- as(read.csv("~/Dropbox/UBC/Projects/Orca_salmon_modelling/07_sizespectra/Viv_sizespectra_model/Salish_Sea/fisheries/time series fisheries/All_fisheries_effort_array_250521.csv", row.names = 1), "matrix")
#relative_effort <- sweep(All_fisheries_effort_array, 2, All_fisheries_effort_array["2011", ], "/")
#new_rel_effort<-head(relative_effort,-7)
#new_rel_effort[as.character(1986:2011), ]

# 
# # Run complete analysis
 results <- run_complete_analysis(
   base_model = base_model,
   historical_effort = new_rel_effort,
   survey_data_path = "path/to/survey_data.csv",
   output_dir = "results/"
 )
# 
# # Access results
final_model <- results$model
historical_sim <- results$historical_simulation
future_projections <- results$projections
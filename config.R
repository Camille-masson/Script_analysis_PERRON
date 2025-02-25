### PROJECT CONFIGURATION ###

# Define paths
root_dir <- getwd()  # Automatically retrieves the project path
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "outputs")
raster_dir <- file.path(root_dir, "raster")
functions_dir <- file.path(root_dir, "Functions")

# Load libraries
library(tidyverse)  # Includes ggplot2, dplyr, etc.
library(lubridate)
library(sf)
library(sp)
library(terra)
library(viridisLite)

# Automatically load functions
# Automatically load function files
function_files <- list.files(functions_dir, pattern = "\\.R$", full.names = TRUE)

# Check if files exist before sourcing them
if (length(function_files) > 0) {
  for (file in function_files) {
    tryCatch(
      source(file),
      error = function(e) {
        message(paste("Error loading:", file))
        message(e)
      }
    )
  }
} else {
  warning("No function files found in Functions/")
}
sapply(function_files, source)

# Global parameters
ncores <- parallel::detectCores() / 3  # Optimal CPU core utilization

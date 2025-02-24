
#### 0. LIBRARIES AND CONSTANTS ####
#----------------------------------#
gc()

# Loading the configuration
source("config.R")

# Defining the analysis year
YEAR <- 9999
ALPAGES_TOTAL <- list(
  "9999" = c("Alpage_demo"),
  "2022" = c("Ane-et-Buyant", "Cayolle", "Combe-Madame", "Grande-Fesse", "Jas-des-Lievres", "Lanchatra", "Pelvas", "Sanguiniere", "Viso"),
  "2023" = c("Cayolle", "Crouzet", "Grande-Cabane", "Lanchatra", "Rouanette", "Sanguiniere", "Vacherie-de-Roubion", "Viso"),
  "2024" = c("Viso", "Cayolle", "Sanguiniere")
)
ALPAGES <- ALPAGES_TOTAL[[as.character(YEAR)]]

#### 1. DATA SIMPLIFICATION IN GPKG ####
#--------------------------------------#

if (FALSE) {  # Set to TRUE to execute
  library(terra)
  source(file.path(functions_dir, "Functions_filtering.R"))
  
  ## INPUT ##
  # A folder containing the raw trajectories, in CSV format from the catlog collars, stored in subfolders named after their pastures
  raw_data_dir <- file.path(data_dir, paste0("Collars_", YEAR, "_raw"))
  
  # The pastures to be processed
  alpages <- c("Alpage_demo")
  
  ## OUTPUT ##
  
  # Creation of the output subfolder: GPS_simple_GPKG
  gps_output_dir <- file.path(output_dir, "GPS_simple_GPKG")
  if (!dir.exists(gps_output_dir)) {
    dir.create(gps_output_dir, recursive = TRUE)
  }
  # Creation of the output GPKG named: Raw_data_9999_Alpage_demo_simplified.gpkg
  output_file <- file.path(gps_output_dir, paste0("Raw_data_", YEAR, "_", alpages, "_simplified.gpkg"))
  
  
  ## CODE ##
  
  lapply(alpages, function(alpage) {
    collar_dir <- file.path(raw_data_dir, alpage) # GPS folder path containing the files of the pasture
    collar_files <- list.files(collar_dir, full.names = TRUE) # List of CSV files in this folder (the different collars)
    lapply(collar_files, function(collar_f) {
      collar_ID <- strsplit(basename(collar_f), split = "_")[[1]][1] # Extracting the collar ID (format = C00_000000 -> C00)
      load_catlog_data(collar_f) %>% # Loading raw GPS data
        slice(which(row_number() %% 30 == 10)) %>% # Keep one GPS point every 30 measurements
        mutate(ID = collar_ID, date = lubridate::format_ISO8601(date)) %>% # Convert date and time to the correct format using lubridate
        vect(geom = c("lon", "lat"), crs = CRS_WSG84) # Transform the dataframe into a spatial object (crs = WSG84)
    }) %>% do.call(rbind, .) # Merge data from different collars
  }) %>% do.call(rbind, .) %>%
    writeVector(filename = output_file, overwrite = TRUE) # Data in GPKG format
}


##GOOD##

#### 2. BJONERAAS FILTER CALIBRATION ####
#--------------------------------------#
if (F) {
  # Loading necessary functions
  source(file.path(functions_dir, "Functions_filtering.R"))
  source(file.path(functions_dir, "Functions_map_plot.R"))
  
  ## INPUTS ##
  # A folder containing the raw trajectories, in CSV format from Catlog collars,
  # stored in subfolders named after their pastures
  raw_data_dir = file.path(data_dir, paste0("Collars_", YEAR, "_raw"))
  
  # A data.frame containing the collar deployment and removal dates
  # Must contain the columns "alpage", "deployment_date" and "removal_date"
  AIF <- file.path(raw_data_dir, paste0(YEAR, "_pasture_info.csv"))
  
  # Checking the file format (often misformatted, be careful with UTF-8 CSV)
  AIF_data <- read.csv(AIF, sep = ",", header = TRUE, row.names = NULL, check.names = FALSE, encoding = "UTF-8")
  
  # The pasture to be processed
  alpage = "Alpage_demo"
  
  ## OUTPUT ##
  # Creating a subfolder to store the results of the Bjorneraas filter
  filter_output_dir <- file.path(output_dir, "Bjorneraas_filter")
  if (!dir.exists(filter_output_dir)) {
    dir.create(filter_output_dir, recursive = TRUE)
  }
  # Output PDF file for filtering visualization
  pdf(file.path(filter_output_dir, paste0("Filtering_calibration_", YEAR, "_", alpage, ".pdf")), width = 9, height = 9)
  
  ## CODE ##
  # List of raw data files for the pasture
  files <- list.files(file.path(raw_data_dir, alpage), full.names = TRUE)
  files <- files[1:3]  # Selecting the first three files (avoids memory overload)
  
  # Loading and concatenating data from the selected files
  data <- do.call(rbind, lapply(files, function(file) { 
    data <- load_catlog_data(file) # Loading CSV files with the dedicated function
    data$ID <- file  # Adding the file identifier to trace its origin
    return(data) 
  }))
  
  # Retrieving collar deployment and removal dates
  beg_date = as.POSIXct(get_alpage_info(alpage, AIF, "deployment_date"), tz="GMT", format="%d/%m/%Y %H:%M")
  end_date = as.POSIXct(get_alpage_info(alpage, AIF, "removal_date"), tz="GMT", format="%d/%m/%Y %H:%M")
  data = date_filter(data, beg_date, end_date) # Filtering data based on validity dates
  
  # Projecting data into Lambert93 (EPSG:2154) from WGS84 (EPSG:4326)
  data_xy <- data %>%
    terra::vect(crs="EPSG:4326") %>%
    terra::project("EPSG:2154") %>%
    as.data.frame(geom = "XY")
  
  # Histogram of time intervals between GPS points
  temps <- diff(data_xy$date)
  temps <- as.numeric(temps, units = "mins")
  hist(temps, nclass = 30)
  
  # Histogram of distances traveled between two successive positions
  dist <- sqrt(diff(data_xy$x)^2+diff(data_xy$y)^2)
  h <- hist(dist, nclass = 30, xlab='Distance (m)', xaxt="n")
  
  ### Testing different filters ###
  # Defining test values for the Bjorneraas filter
  medcrits = c(750, 500, 750) # Median threshold for abnormal distances
  meancrits = c(500, 500, 350) # Mean threshold for abnormal distances
  spikesps = c(1500, 1500, 1500) # Maximum speed threshold (m/h)
  spikecoss = c(-0.95, -0.95, -0.95) # Direction change threshold (cosine of angle)
  
  for (i in 1:length(medcrits)) {
    # Applying the Bjorneraas filter with each parameter combination
    trajectories <- position_filter(data, medcrit=medcrits[i], meancrit=meancrits[i], spikesp=spikesps[i], spikecos=spikecoss[i])
    
    # Determining spatial limits for the map
    minmax_xy = get_minmax_L93(trajectories[!(trajectories$R1error | trajectories$R2error ),], buffer = 100)
    
    # Assigning error codes for visualization
    trajectories$errors = 1
    trajectories$errors[trajectories$R1error] = 2
    trajectories$errors[trajectories$R2error] = 3
    
    # Defining the color palette for the map
    pal <- c("#56B4E9", "red", "black")
    
    # Displaying GPS trajectories with detected errors
    print(ggplot(trajectories, aes(x, y, col = errors)) +
            geom_path(size = 0.2) +
            geom_point(size = 0.3) +
            coord_equal() +
            xlim(minmax_xy$x_min, minmax_xy$x_max) + ylim(minmax_xy$y_min, minmax_xy$y_max) +
            ggtitle(paste0("medcrit = ", medcrits[i], ", meancrit = ", meancrits[i], ", spikesp = ", spikesps[i], ", spikecos = ", spikecoss[i])) +
            scale_colour_gradientn(colors=pal, guide="legend", breaks = c(1, 2, 3), labels = c("OK", "R1error", "R2error")))
    
    # Zooming in on the first five days after deployment
    trajectories <- trajectories %>%
      filter(date < beg_date + 3600*24*5)
    print(ggplot(trajectories, aes(x, y, col = errors)) +
            geom_path(size = 0.2) +
            geom_point(size = 0.3) +
            coord_equal() +
            ggtitle(paste0("5 days only, medcrit = ", medcrits[i], ", meancrit = ", meancrits[i], ", spikesp = ", spikesps[i], ", spikecos = ", spikecoss[i])) +
            scale_colour_gradientn(colors=pal, guide="legend", breaks = c(1, 2, 3), labels = c("OK", "R1error", "R2error")))
  }
  
  # Closing the PDF file containing the visualizations
  dev.off()
}

##GOOD##

#### 3. FILTERING CATLOG DATA ####
#--------------------------------#

if (F) {
  # Loading necessary functions
  source(file.path(functions_dir, "Functions_filtering.R"))
  
  ## INPUTS ##
  # A folder containing the raw trajectories, in CSV format from Catlog collars, stored in subfolders named after their pastures. Coordinates are in WSG84. The file names, in the format "ID_something.csv", will be used to determine the collar ID, which must have 3 characters.
  raw_data_dir = file.path(data_dir, paste0("Collars_", YEAR, "_raw"))
  
  # A file containing information on each equipped individual, including the deployment and removal dates of the collars, as well as the proportion of time the collars are programmed to be active (18h per day = 0.75). Must contain the columns "Collier", "Alpage", "Espece", "Race", "date_pose", "date_retrait", and "proportion_jour_allume".
  IIF = file.path(raw_data_dir, paste0(YEAR, "_collars_deployed.csv"))
  
  # Load and verify collar deployment data (format)
  IIF_data <- read.csv(IIF, stringsAsFactors = FALSE, encoding = "UTF-8")
  
  # Pastures to be processed
  alpages = c("Alpage_demo")
  
  ## OUTPUTS ##
  filter_output_dir <- file.path(output_dir, "Bjorneraas_filter")
  if (!dir.exists(filter_output_dir)) {
    dir.create(filter_output_dir, recursive = TRUE)
  }
  
  # An .RDS file containing the filtered trajectories (new trajectories are appended to previously processed trajectories). Coordinates are in Lambert93.
  output_rds_file = file.path(filter_output_dir, paste0("Catlog_", YEAR, "_filtered_", alpages, ".rds"))
  
  # A .csv file containing collar performance metrics (percentages of points removed at each step, defective collars, etc.).
  indicator_file = file.path(filter_output_dir, paste0(YEAR, "_filtering_", alpages, ".csv"))
  
  ## CODE ##
  
  for (alpage in alpages) {
    print(paste("WORKING ON PASTURE", alpage))
    collar_dir <- file.path(raw_data_dir, alpage)
    collar_file <- list.files(collar_dir, pattern = ".csv", full.names = TRUE)
    
    medcrit = get_alpage_info(alpage, AIF, "medcrit")
    meancrit = get_alpage_info(alpage, AIF, "meancrit")
    spikesp = get_alpage_info(alpage, AIF, "spikesp")
    spikecos = as.numeric(gsub(",", ".", get_alpage_info(alpage, AIF, "spikecos")))
    print(paste0("Bjorneraas filter parameters: medcrit=", medcrit, ", meancrit=", meancrit, ", spikesp=", spikesp, ", spikecos=", spikecos))
    print(collar_files)
    
    # Filtering trajectories and calculating indicators
    indicators <- lapply(collar_files, function(collar) {
      filter_one_collar(
        load_catlog_data(collar),  
        basename(collar),  # Only pass the file name
        output_rds_file, alpage, beg_date, end_date, IIF,
        bjoneraas.medcrit = medcrit,
        bjoneraas.meancrit = meancrit,
        bjoneraas.spikesp = spikesp,
        bjoneraas.spikecos = spikecos
      )
    }) %>%
      do.call(rbind, .)
    
    indicators_tot = indicators %>%
      filter(worked_until_end == 1) %>% # To compute performance indicators at the pasture level, we remove defective collars
      add_row(name = paste("TOTAL", alpage), worked_until_end = sum(.$worked_until_end), nloc = NA,
              R1error = NA, R2error = NA,
              error_perc = sum(.$nloc * .$error_perc) / sum(.$nloc), localisation_rate = mean(.$localisation_rate))
    indicators = rbind(indicators, indicators_tot[nrow(indicators_tot),])
    
    write.table(indicators, file=indicator_file, append = T, sep=',', row.names=F, col.names=F)
  }
}

##GOOD##


#### 4. HMM FITTING #### 
#----------------------#
if (F) {
  library(snow)
  library(stats)
  # Movement modeling packages
  library(momentuHMM)
  library(adehabitatLT)
  library(adehabitatHR)
  # RMarkdown libraries
  library(knitr)
  library(rmarkdown)
  source(file.path(functions_dir, "Functions_HMM_fitting.R"))
  
  ## INPUTS ##
  # An .RDS file containing the filtered trajectories
  input_rds_file <- file.path(output_dir, "Bjorneraas_filter", paste0("Catlog_", YEAR, "_filtered_", alpage, ".rds"))
  
  # A data.frame containing the correspondence between collars and pastures. Must contain the columns "ID", "Alpage", and "Sampling Period".
  individual_info_file <- file.path(data_dir, paste0("Collars_", YEAR, "_raw"), paste0(YEAR, "_collars_deployed.csv"))
  individual_info_file_data <- read.csv(individual_info_file, stringsAsFactors = FALSE, encoding = "UTF-8")
  
  # Pastures to be processed
  alpages = ALPAGES
  
  ## OUTPUTS ##
  # Creating the subfolder to store the HMM behavior results
  filter_output_dir <- file.path(output_dir, "HMM_behavior")
  if (!dir.exists(filter_output_dir)) {
    dir.create(filter_output_dir, recursive = TRUE)
  }
  
  # An .RDS file containing the trajectories categorized by behavior (new trajectories are appended to previously processed ones)
  output_rds_file = file.path(output_dir, "HMM_behavior", paste0("Catlog_", YEAR, "_", alpage, "_viterbi.rds"))
  
  ### LOADING DATA FOR ANALYSIS ###
  data = readRDS(input_rds_file)
  data = data[data$species == "brebis",]
  data = data[data$alpage %in% alpages,]
  
  ### HMM FIT ###
  run_parameters = list(
    # Model
    model = "HMM",
    
    # Resampling
    resampling_ratio = 5,
    resampling_first_index = 0,
    rollavg = FALSE,
    rollavg_convolution = c(0.15, 0.7, 0.15),
    knownRestingStates = FALSE,
    
    # Observation distributions (step lengths and turning angles)
    dist = list(step = "gamma", angle = "vm"),
    # Design matrices to be used for the probability distribution parameters of each data stream
    DM = list(angle=list(mean = ~1, concentration = ~1)),
    # Covariate formula
    covariants = ~cos(hour * 3.141593 / 12), # ~1 if no covariates are used
    
    # 3-state HMM
    Par0 = list(step = c(10, 25, 50, 10, 15, 40), angle = c(tan(pi/2), tan(0/2), tan(0/2), log(0.5), log(0.5), log(3))),
    fixPar = list(angle = c(tan(pi/2), tan(0/2), tan(0/2), NA, NA, NA))
  )
  run_parameters = scale_step_parameters_to_resampling_ratio(run_parameters)
  
  startTime = Sys.time()
  results = par_HMM_fit_test(data, run_parameters, ncores = ncores, individual_info_file, sampling_period = 120, output_dir)
  endTime = Sys.time()
  
  # Check internet connection
  
  ## SUMMARIZE MODEL FITTING BY PASTURE in rmarkdown PDFs (To be reviewed)
  parameters_df <- parameters_to_data.frame(run_parameters)
  
  # To be reviewed!
  individual_IDs <- sapply(results, function(hmm) hmm$data$ID[1])
  individual_alpages <- get_individual_alpage(individual_IDs, individual_info_file)
  for (alpage in unique(individual_alpages)) {
    res_index  <- individual_alpages == alpage
    data_alpage <- do.call("rbind", lapply(results[res_index], function(result) result$data))
    results_df <- do.call("rbind", lapply(results[res_index], hmm_result_to_data.frame))
    
    runs_to_pdf(alpage, parameters_df, results_df, data_alpage , paste(round(difftime(endTime, startTime, units='mins'),2), "min"), output_dir, paste0(output_dir,"1 HMM Fit/HMM_fitting_",alpage,".pdf"), show_performance_indicators = FALSE)
  }
  # To be reviewed! : function Functions/HMM_fit_template.Rmd does not exist
  
  ### SAVE RESULTING TRAJECTORIES ###
  data_hmm <- do.call("rbind", lapply(results, function(result) result$data))
  viterbi_trajectory_to_rds(data_hmm, output_rds_file, individual_info_file)
}


##GOOD##

#### 5. FLOCK STOCKING RATE (Charge) BY DAY AND BY STATE ####
#-------------------------------------------------------------#
library(adehabitatHR)
library(data.table)
library(snow)
source(file.path(functions_dir, "Functions_map_plot.R"))
source(file.path(functions_dir, "Functions_flock_density.R"))

## INPUTS ##
# An .RDS file containing the trajectories categorized by behavior
input_rds_file <- file.path(output_dir, "HMM_behavior", paste0("Catlog_", YEAR, "_", alpage, "_viterbi.rds"))

# A data.frame containing flock sizes and their changes over time based on date
flock_size_file <- file.path(raw_data_dir, paste0(YEAR, "_herd_sizes.csv"))
flock_size_file_data <- read.csv(flock_size_file, stringsAsFactors = FALSE, encoding = "UTF-8")

# Pastures to be processed
alpages <- ALPAGES

## OUTPUTS ##
# Output folder
save_dir <- file.path(output_dir, "Flock_stocking_rate")

# An .RDS per pasture containing daily stocking rates per behavior state
state_daily_rds_prefix <- paste0("by_day_and_state_", YEAR, "_")
# An .RDS per pasture containing daily stocking rates
daily_rds_prefix <- paste0("by_day_", YEAR, "_")
# An .RDS per pasture containing stocking rates by behavior state
state_rds_prefix <- paste0("by_state_", YEAR, "_")
# An .RDS per pasture containing total stocking rate over the entire season
total_rds_prefix <- paste0("total_", YEAR, "_")

h <- 25 # Characteristic distance for calculating stocking rate

for (alpage in alpages) {
  flock_sizes <- get_flock_size_through_time(alpage, flock_size_file)
  prop_time_collar_on <- get_alpage_info(alpage, AIF, "proportion_active_day")
  
  # Load filtered data for the pasture
  data <- readRDS(input_rds_file)
  data <- data[data$alpage == alpage,]
  
  # Load the phenology raster with the correct path
  raster_file <- file.path(raster_dir, paste0("ndvis_", YEAR, "_", alpage, "_pheno_metrics.tif"))
  pheno_t0 <- get_raster_cropped_L93(raster_file, get_minmax_L93(data, 100), reproject = TRUE, band = 2, as = "SpatialPixelDataFrame")
  
  # Define the storage folder specific to the pasture
  alpage_save_dir <- file.path(save_dir, paste0(alpage, "_", YEAR))
  if (!dir.exists(alpage_save_dir)) dir.create(alpage_save_dir, recursive = TRUE)
  
  ## BY DAY AND BY STATE ##
  flock_load_by_day_and_state_to_rds_kernelbb_test_11(
    data, 
    pheno_t0, 
    alpage_save_dir,  
    state_daily_rds_prefix, 
    flock_sizes, 
    prop_time_collar_on
  )
  
  gc()
  
  # Merge individual files
  merged_file <- flock_merge_rds_files_0(alpage_save_dir, state_daily_rds_prefix)
  
  rm(data)
  
  charge <- readRDS(file.path(alpage_save_dir, paste0(state_daily_rds_prefix, alpage, ".rds")))
  unique(charge$state)
  
  ## BY STATE ##
  charge_state <- charge %>%
    group_by(x, y, state) %>%
    summarise(Charge = sum(Charge, na.rm = TRUE), .groups = 'drop') %>%
    as.data.frame()
  saveRDS(charge_state, file.path(alpage_save_dir, paste0(state_rds_prefix, alpage, ".rds")))
  rm(charge_state)
  
  ## BY DAY ##
  charge_day <- lapply(unique(charge$day), function(d) {
    charge %>%
      filter(day == d) %>%
      group_by(x, y, day) %>%
      summarise(Charge = sum(Charge, na.rm = TRUE), .groups = 'drop')
  })
  charge_day <- as.data.frame(rbindlist(charge_day, use.names = TRUE))
  saveRDS(charge_day, file.path(alpage_save_dir, paste0(daily_rds_prefix, alpage, ".rds")))
  rm(charge_day)
  
  ## TOTAL STOCKING RATE ##
  charge_tot <- charge %>%
    group_by(x, y) %>%
    summarise(Charge = sum(Charge, na.rm = TRUE), .groups = 'drop') %>%
    as.data.frame()
  saveRDS(charge_tot, file.path(alpage_save_dir, paste0(total_rds_prefix, alpage, ".rds")))
  rm(charge_tot)
  
  rm(charge)
}


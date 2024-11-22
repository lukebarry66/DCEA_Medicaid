#------------------------------------------------------------------------------#
#---------------------Purpose: MSM DCEA Medicaid-CVD---------------------------#
#---------------------------------GLOBAL---------------------------------------#
#------------------------------------------------------------------------------#


################# Run main analysis across computer clusters ###################

#### Model iterations by cluster 
pacman::p_load(cleaR, here)
clear()

# number of bootstrapped model iterations; default 1000
bootstart    <- 1
bootsize     <- 10
# Adjust this size based on your system capacity to do bootstrapping in chunks
chunk_size   <- ((bootsize-bootstart)+1)/2
# Number of imputed datasets; default m = 10 (use current_m to analyze subsets of imputed datasets)
current_m    <- 1
m            <- 1

source(here("scripts", "5-Analysis.R"))

############################## Combine results #################################

# Define the path to the data directory
data_dir <- here("data", "output_data")

# Define the dataset numbers you want to process (1, 2, 3, etc.)
dataset_numbers <- current_m:m  # Adjust this based on the range you want

# Initialize an empty list to store the results
all_datasets <- list()

# Loop through each dataset number and process the corresponding RDS files
for (dataset_number in dataset_numbers) {
  # Generate the file pattern for each datase
  file_pattern <- paste0("^output_Dataset", dataset_number, "_boot")
  
  # Get the list of RDS files matching the current pattern
  files <- list.files(data_dir, pattern = file_pattern, full.names = TRUE)
  
  # Initialize an empty list to store all items for this dataset
  all_items <- list()
  
  # Process each RDS file for the current dataset
  for (file in files) {
    # Read the RDS file
    data <- readRDS(file)
    
    # Construct the dataset key 
    dataset_key <- paste0("Dataset_", dataset_number)
    
    # Check if the dataset key exists in the data
    if (dataset_key %in% names(data)) {
      # Combine the chunks into one list and remove the chunk names
      combined_chunks <- unlist(data[[dataset_key]], recursive = FALSE)
      names(combined_chunks) <- sub("^[^\\.]+\\.", "", names(combined_chunks))  # Remove chunk names from the names of items
      
      # Add the replication and dataset_index variables to each item
      combined_chunks <- lapply(names(combined_chunks), function(item_name) {
        # Extract the replication value from the item name
        replication_value <- as.numeric(sub(".*\\.(\\d+)$", "\\1", item_name))  # Get the number after the dot
        
        # Add 'replication' and 'dataset_index' fields to each item
        item <- combined_chunks[[item_name]]
        item$replication <- replication_value  # Assign the actual item number as 'replication'
        item$dataset_index <- dataset_number  # Assign the dataset index (dataset_number)
        return(item)
      })
      
      # Add the combined chunks to the all_items list
      all_items <- c(all_items, combined_chunks)
    } else {
      message(dataset_key, " not found in file: ", file)
    }
  }
  
  # Assign the combined list to the appropriate key in all_datasets
  all_datasets[[paste0("Dataset", dataset_number)]] <- all_items
}

# Now let's combine all the items from all datasets into a single data frame
# Initialize an empty list to store the combined data
combined_df_list <- list()

# Loop through all datasets in all_datasets
for (dataset_name in names(all_datasets)) {
  # Extract the current dataset's items
  dataset_items <- all_datasets[[dataset_name]]
  
  # Convert the list of items to a data frame
  dataset_df <- bind_rows(dataset_items)
  
  # Add the dataset's data frame to the combined_df_list
  combined_df_list[[dataset_name]] <- dataset_df
}

# Combine all dataset data frames into a single data frame
combined_output <- bind_rows(combined_df_list) %>%
  filter(!is.na(variable)) %>%
  mutate(  diff_mi        = mi_trt - mi_ntrt
           , diff_stroke    = stroke_trt - stroke_ntrt
           , diff_CVDdeath  = CVD_death_trt - CVD_death_ntrt
           , diff_nCVDdeath = nCVD_death_trt - nCVD_death_ntrt
           , inc_cost       = tc_trt - tc_ntrt
           , inc_effect     = te_trt - te_ntrt
           , nhb_ntrt       = (te_ntrt) - ((tc_ntrt) / wtp)
           , nhb_trt        = (te_trt) - ((tc_trt) / wtp)
           , icer           = (tc_trt - tc_ntrt) / (te_trt - te_ntrt)
           , inhb           = (te_trt - te_ntrt) - ((tc_trt - tc_ntrt) / wtp)
           , inmb           = ((te_trt - te_ntrt) * wtp) - (tc_trt - tc_ntrt)
  )

## Save the final combined data frame to a new RDS file
saveRDS(combined_output, here("data", "output_data", glue("combined_output", "_{today()}.rds")))


############################# Summarize results ################################

#source the combined_output
combined_output <- readRDS(here("data", "output_data", glue("combined_output", "_{today()}.rds")))

# run results synthesis
source(here("scripts", "6-Results.R"))



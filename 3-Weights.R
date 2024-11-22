#------------------------------------------------------------------------------#
#---------------------Purpose: MSM DCEA Medicaid-CVD---------------------------#
#-----------------------------Weighted counts----------------------------------#
#------------------------------------------------------------------------------#

set.seed(seed)

#### Weighted counts ----

.weight <- "WTSAF2YR"
.psu    <- "SDMVPSU"
.strata <- "SDMVSTRA"
options(survey.lonely.psu = "adjust") 

#--------------------Creating functions--------------------------
#--Create function to display point estimate and confidence interval
ci_fct <- function(pe, low, high, multiplier=1, digits=1, scale="%", count=FALSE){
  if(count==F){
    res <- paste(multiplier*formattable::digits(pe,digits),scale, " " ,"(", 
                 multiplier*formattable::digits(low,digits),",", 
                 multiplier*formattable::digits(high,digits), ")", sep="")
  } else {
    res <- paste(formattable::comma(pe, digits=digits), " " ,"(", 
                 formattable::comma(low, digits=digits),", ", 
                 formattable::comma(high, digits=digits), ")", sep="")
  }
  return(res)
}

#--Test the function
ci_fct(pe=2, low=1, high=3, multiplier=100, scale="%")
ci_fct(pe=200000, low=200000, high=200000, multiplier=NULL, digits=0, scale=NULL, count=T)

#--Create function to calculate the count

count_function <- function(.data, variable){
  .data %>% 
    count(!!sym(variable)) %>% 
    rename(value = variable,
           unweighted_count=n) %>% 
    mutate(variables=variable) %>% 
    relocate(variables, value, unweighted_count)
}

#to create a binary variable
d <- function(expr){
  as.numeric(expr)
}

#--------------------Obtaining weighted and unweighted counts-------------------

##data prep
nhanes_wt <- nhanes_full %>%
  filter(!is.na(.data[[.weight]])
         , age <= .stop_age
         , age >= .start_age
         , wave_year == "H_13") %>% # Restrict to wave relelvant to Medicaid expansion
  # create a family income variable to correspond to MEPS cost and utility calculations and according to the poverty line 
  #see, https://meps.ahrq.gov/survey_comp/hc_technical_notes.shtml
  mutate(medicaid_exp = ifelse(insurance  == "uninsured" & fpl < .fpl & age < 65, 1, 0)) %>%
  #na.omit() %>%
  filter(medicaid_exp == 1 | .ATT != 1) %>%
  select(  .data[[.weight]]
         , .data[[.psu]]
         , .data[[.strata]]) %>%
  mutate(  total = 1) 

##Estimating weighted weighted_count----
weighted_count_data_ci <- nhanes_wt %>% 
  as_survey_design(weights = .data[[.weight]]
                   , strata = .data[[.strata]]
                   , ids = .data[[.psu]]
                   , nest = TRUE
                   ) %>% 
  summarize_all(survey_total, vartype = "ci", na.rm=T) %>% 
  t() %>%
  data.frame(weighted_count=.) %>%
  rownames_to_column(var = "variables")

weighted_count_data_low <- weighted_count_data_ci %>% 
  filter(str_detect(variables, '_low')) %>% 
  mutate(variables=str_replace_all(variables,"_low", "")) %>%
  rename(weighted_count_low=weighted_count)

weighted_count_data_high <- weighted_count_data_ci %>% 
  filter(str_detect(variables, '_upp')) %>% 
  mutate(variables=str_replace_all(variables,"_upp", "")) %>%
  rename(weighted_count_high=weighted_count)

weighted_count_data_pe <- weighted_count_data_ci %>% 
  filter(variables %in% c(
    "total")) %>% 
  rename(weighted_count_pe=weighted_count)

weighted_count_data <- weighted_count_data_pe %>% 
  full_join(weighted_count_data_low, by=c("variables")) %>% 
  full_join(weighted_count_data_high, by=c("variables")) 

table1_us_weighted_count <- weighted_count_data %>% 
  transmute(
    `Pop. characteristic` = case_when(variables== "total" ~ "Total"
                                      , TRUE	~	NA_character_),
    `NHANES weighted_count`=ci_fct(pe=weighted_count_pe, low=weighted_count_low, 
                                  high=weighted_count_high, multiplier=NULL, digits=0, scale=NULL, count=T),
    `NHANES weighted_count pe`= paste(formattable::comma(weighted_count_pe, digits = 0), sep=""))

##Outputing table 1 - weighted ----
table1_us_weighted_count %>% 
  as_hux()

write_xlsx(table1_us_weighted_count, 
           here("tables", "table1_us_overall_weighted_count.xlsx"))

#### Population counts ----

# create weighted population counts to merge
pop_counts <- weighted_count_data_pe %>%
  as_tibble() %>%
  rename(category = variables,
         pop_counts = weighted_count_pe) %>%
  mutate(category = case_when(  category == "total" ~ "Total"),
         category = as_factor(category)) 

#### Replication weights for bootrapped analysis

if (.rw == 1) {
  # Create survey object using HCUP complex survey weights
  nhanes_rw1 <- svydesign(id=~SDMVPSU
                          , weights=~WTSAF
                          , strata=~SDMVSTRA
                          , nest=T
                          , survey.lonely.psu = "adjust"
                          , data=nhanes_full %>% filter(!is.na(WTSAF)))
  svymean(~age, nhanes_rw1)
  
  # Create the bootstrap replicate weights
  nhanes_rw2 <- as.svrepdesign(nhanes_rw1
                               , type = "bootstrap"
                               , replicates = n.rw
                               , compress = TRUE)
  rm(nhanes_rw1)
  
  # Extract the replicate weights and convert to a data frame
  replicates <- weights(nhanes_rw2, type = "analysis")
  replicates <- as.data.frame(replicates)
  replicates <- as_tibble(replicates) %>%
    mutate(wave_year = nhanes_rw2$variable$wave_year
           , SEQN = nhanes_rw2$variable$SEQN)
  rm(nhanes_rw2)
  
  write_rds(replicates, here("data", "output_data", glue("replicate_weights_n{n.rw}.rds")))
}

## test replicate weights
#execution_time <- system.time({
#  # test and compare with previous estimates
#  boot_results <- list()
#  for (i in 1:1000) {
#    # Randomly select a column from V1 to V100
#    random_col <- sample(paste0("V", 1:n.rw), 1)
#    temp <- nhanes %>%
#      right_join(replicates %>% select(SEQN, wave_year, !!sym(random_col))
#                 , by = c("SEQN", "wave_year")) %>% 
#      slice_sample(., prop = 1, weight_by = !!sym(random_col), replace = T)
#    
#    boot_results[[i]] <- mean(temp$age, na.rm = T)
#  }
#  # Combine the bootstrap results into a single vector
#  boot_bind <- unlist(boot_results)
#  print(mean(boot_bind))
#  print(sd(boot_bind))
#})
## Print the result and execution time
#print(execution_time)#
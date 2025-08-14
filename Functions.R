### VKT Elasticities - Functions

# Load libraries
library(haven)
library(lmtest)
library(sandwich)
library(dplyr)
library(plm)
library(fixest)
library(tigris)
library(tidyr)
library(stargazer)
library(readr)
library(mapview)
library(sf)
library(ivreg)
library(modelsummary)

# Set working directory to folder location
setwd("/Users/Username/Downloads/VMT_Elasticities/Data/")

# Suppress dplyr summaries
options(dplyr.summarise.inform = FALSE)

# Function to clean HPMS data
clean_hpms <- function(filepath_universe,
                       state_code_name_universe, 
                       county_code_name_universe, 
                       section_id_name_universe, 
                       rural_urban_name_universe, 
                       f_class_name_universe, 
                       length_name_universe, 
                       aadt_name_universe, 
                       through_lanes_name_universe,
                       filepath_sample,
                       state_code_name_sample, 
                       county_code_name_sample, 
                       section_id_name_sample, 
                       rural_urban_name_sample, 
                       f_class_name_sample, 
                       length_name_sample, 
                       aadt_name_sample, 
                       through_lanes_name_sample,
                       expansion_factor_name_sample,
                       year) {
  
  # Read universe data
  hpms_universe <- read.csv(filepath_universe) %>%
    # Rename variable names to be consistent (they vary by year)
    rename(State_Code = state_code_name_universe,
           County_Code = county_code_name_universe,
           Section_ID = section_id_name_universe,
           Rural_Urban = rural_urban_name_universe,
           F_Class = f_class_name_universe,
           Section_Length = length_name_universe,
           AADT = aadt_name_universe,
           Through_Lanes = through_lanes_name_universe) %>%
    # Calculate VKT and Lane Kilometers
    mutate(VKT = AADT * Section_Length * 1.609344,
           Lane_Km = Through_Lanes * Section_Length * 1.609344,
           # Categorize segment by urban/rural classification and functional classification
           Segment_Category = case_when(
             (Rural_Urban == 3 | Rural_Urban == 4) & (F_Class == 1 | F_Class == 11) ~ "Urban_Class1",
             (Rural_Urban == 3 | Rural_Urban == 4) & (F_Class == 2 | F_Class == 12) ~ "Urban_Class2",
             (Rural_Urban == 3 | Rural_Urban == 4) & (F_Class == 6 | F_Class == 14) ~ "Urban_Class3",
             (Rural_Urban == 3 | Rural_Urban == 4) & (F_Class == 7 | F_Class == 16) ~ "Urban_Class4",
             (Rural_Urban == 3 | Rural_Urban == 4) & (F_Class == 8 |F_Class == 9 | F_Class == 17 | F_Class == 19) ~ "Urban_Local",
             (Rural_Urban == 1 | Rural_Urban == 2) & (F_Class == 1 | F_Class == 11) ~ "NonUrban_Class1",
             (Rural_Urban == 1 | Rural_Urban == 2) & (F_Class == 2 | F_Class == 12) ~ "NonUrban_Class2",
             (Rural_Urban == 1 | Rural_Urban == 2) & (F_Class == 6 | F_Class == 14) ~ "NonUrban_Class3",
             (Rural_Urban == 1 | Rural_Urban == 2) & (F_Class == 7 | F_Class == 16) ~ "NonUrban_Class4",
             (Rural_Urban == 1 | Rural_Urban == 2) & (F_Class == 8 | F_Class == 9 | F_Class == 17 | F_Class == 19) ~ "NonUrban_Local",
             TRUE ~ "Other"
           )) %>%
    select(State_Code, County_Code, Section_ID, F_Class, Section_Length, AADT, Through_Lanes, VKT, Lane_Km, Segment_Category) %>%
    # Create a unique segment id tat incorporates the existing segment ID and state/county codes
    mutate(unique_id = paste0(State_Code, "_", County_Code, "_", Section_ID)) %>%
    # Group by the unique_id and segment category and sum VKT/Lane Kilometers
    group_by(unique_id, Segment_Category) %>%
    summarise(VKT = sum(VKT, na.rm = T),
              Lane_Km = sum(Lane_Km, na.rm = T)) %>%
    ungroup()
  
  print(paste0("Finished cleaning ", year, " hpms universe data."))
  
  # Read sample data
  hpms_sample <- read.csv(filepath_sample) %>%
    # Rename variable names to be consistent (they vary by year)
    rename(State_Code = state_code_name_sample,
           County_Code = county_code_name_sample,
           Section_ID = section_id_name_sample,
           Rural_Urban = rural_urban_name_sample,
           F_Class = f_class_name_sample,
           Section_Length = length_name_sample,
           AADT = aadt_name_sample,
           Through_Lanes = through_lanes_name_sample,
           Expansion_Factor = expansion_factor_name_sample) %>%
    # Calculate VKT and Lane Kilometers, incorporating the sample's expansion factor
    mutate(VKT_Sample = AADT * (Section_Length * Expansion_Factor) * 1.609344,
           Lane_Km_Sample = Through_Lanes * (Section_Length * Expansion_Factor) * 1.609344,
           # Categorize segment by urban/rural classification and functional classification
           Segment_Category_Sample = case_when(
             (Rural_Urban == 3 | Rural_Urban == 4) & (F_Class == 1 | F_Class == 11) ~ "Urban_Class1",
             (Rural_Urban == 3 | Rural_Urban == 4) & (F_Class == 2 | F_Class == 12) ~ "Urban_Class2",
             (Rural_Urban == 3 | Rural_Urban == 4) & (F_Class == 6 | F_Class == 14) ~ "Urban_Class3",
             (Rural_Urban == 3 | Rural_Urban == 4) & (F_Class == 7 | F_Class == 16) ~ "Urban_Class4",
             (Rural_Urban == 3 | Rural_Urban == 4) & (F_Class == 8 | F_Class == 9 | F_Class == 17 | F_Class == 19) ~ "Urban_Local",
             (Rural_Urban == 1 | Rural_Urban == 2) & (F_Class == 1 | F_Class == 11) ~ "NonUrban_Class1",
             (Rural_Urban == 1 | Rural_Urban == 2) & (F_Class == 2 | F_Class == 12) ~ "NonUrban_Class2",
             (Rural_Urban == 1 | Rural_Urban == 2) & (F_Class == 6 | F_Class == 14) ~ "NonUrban_Class3",
             (Rural_Urban == 1 | Rural_Urban == 2) & (F_Class == 7 | F_Class == 16) ~ "NonUrban_Class4",
             (Rural_Urban == 1 | Rural_Urban == 2) & (F_Class == 8 | F_Class == 9 | F_Class == 17 | F_Class == 19) ~ "NonUrban_Local",
             TRUE ~ "Other"
           )) %>%
    select(State_Code, County_Code, Section_ID, F_Class, Section_Length, AADT, Through_Lanes, VKT_Sample, Lane_Km_Sample, Segment_Category_Sample) %>%
    # Create a unique segment id tat incorporates the existing segment ID and state/county codes
    mutate(unique_id = paste0(State_Code, "_", County_Code, "_", Section_ID)) %>%
    # Group by the unique_id and segment category and sum VKT/Lane Kilometers
    group_by(unique_id, Segment_Category_Sample) %>%
    summarise(VKT_Sample = sum(VKT_Sample, na.rm = T),
              Lane_Km_Sample = sum(Lane_Km_Sample, na.rm = T)) %>%
    ungroup()
  
  print(paste0("Finished cleaning ", year, " hpms sample data."))
  
  # Merge hpms universe and sample data
  hpms <- merge(hpms_universe,
                hpms_sample,
                by = "unique_id",
                all = T) %>%
    # Calculate a variable for the max value of VKT and Lane Kilometers across universe and sample data
    mutate(VKT_Max = pmax(VKT, VKT_Sample, na.rm = T),
           Lane_Km_Max = pmax(Lane_Km, Lane_Km_Sample, na.rm = T),
           # If using sample data, use the sample segment category value
           Segment_Category = ifelse(is.na(Segment_Category), Segment_Category_Sample, Segment_Category)) %>%
    # Group by the unique_id and segment category and find the maximum values for VKT and Lane Kilometers
    # There should only be one value per unique_id/Segment_Category pairing, so this is mostly unecessary
    group_by(unique_id, Segment_Category) %>%
    summarise(VKT = max(VKT_Max, na.rm = T),
              Lane_Km = max(Lane_Km_Max, na.rm = T)) %>%
    ungroup() %>%
    # Add state/county FIPS codes from unique_id
    mutate(State = sub("_.*", "", unique_id),
           County = sub("^[^_]+_([^_]+).*", "\\1", unique_id),
           County_ID = paste0(State, "_", County)) %>%
    # Group by MSA and Segment Category and calculate total VKT/Lane Kilometers
    group_by(County_ID, Segment_Category) %>%
    summarise(VKT = sum(VKT, na.rm = T),
              Lane_Km = sum(Lane_Km, na.rm = T)) %>%
    ungroup() %>%
    # Reshape the dataframe so that each row is equal to a single MSA
    pivot_wider(names_from = Segment_Category, values_from = c(VKT, Lane_Km)) %>%
    select(County_ID,
           VKT_Urban_Class1, VKT_Urban_Class2, VKT_Urban_Class3, VKT_Urban_Class4, VKT_Urban_Local, VKT_NonUrban_Class1, VKT_NonUrban_Class2, VKT_NonUrban_Class3, VKT_NonUrban_Class4, VKT_NonUrban_Local,
           Lane_Km_Urban_Class1, Lane_Km_Urban_Class2, Lane_Km_Urban_Class3, Lane_Km_Urban_Class4, Lane_Km_Urban_Local, Lane_Km_NonUrban_Class1, Lane_Km_NonUrban_Class2, Lane_Km_NonUrban_Class3, Lane_Km_NonUrban_Class4, Lane_Km_NonUrban_Local)
  
  hpms[is.na(hpms)] <- 0
  
  hpms_out <- hpms %>%
    mutate(VKT_Class1_Total = VKT_Urban_Class1 + VKT_NonUrban_Class1,
           VKT_Class2_Total = VKT_Urban_Class2 + VKT_NonUrban_Class2,
           VKT_Class3_Total = VKT_Urban_Class3 + VKT_NonUrban_Class3,
           VKT_Class4_Total = VKT_Urban_Class4 + VKT_NonUrban_Class4,
           VKT_Local_Total = VKT_Urban_Local + VKT_NonUrban_Local,
           
           VKT_Urban_Total = VKT_Urban_Class1 + VKT_Urban_Class2 + VKT_Urban_Class3 + VKT_Urban_Class4 + VKT_Urban_Local,
           VKT_NonUrban_Total = VKT_NonUrban_Class1 + VKT_NonUrban_Class2 + VKT_NonUrban_Class3 + VKT_NonUrban_Class4 + VKT_NonUrban_Local,
           
           VKT_Total = VKT_Class1_Total + VKT_Class2_Total + VKT_Class3_Total + VKT_Class4_Total + VKT_Local_Total,
           
           Lane_Km_Class1_Total = Lane_Km_Urban_Class1 + Lane_Km_NonUrban_Class1,
           Lane_Km_Class2_Total = Lane_Km_Urban_Class2 + Lane_Km_NonUrban_Class2,
           Lane_Km_Class3_Total = Lane_Km_Urban_Class3 + Lane_Km_NonUrban_Class3,
           Lane_Km_Class4_Total = Lane_Km_Urban_Class4 + Lane_Km_NonUrban_Class4,
           Lane_Km_Local_Total = Lane_Km_Urban_Local + Lane_Km_NonUrban_Local,
           
           Lane_Km_Urban_Total = Lane_Km_Urban_Class1 + Lane_Km_Urban_Class2 + Lane_Km_Urban_Class3 + Lane_Km_Urban_Class4 + Lane_Km_Urban_Local,
           Lane_Km_NonUrban_Total = Lane_Km_NonUrban_Class1 + Lane_Km_NonUrban_Class2 + Lane_Km_NonUrban_Class3 + Lane_Km_NonUrban_Class4 + Lane_Km_NonUrban_Local,
           
           Lane_Km_Total = Lane_Km_Class1_Total + Lane_Km_Class2_Total + Lane_Km_Class3_Total + Lane_Km_Class4_Total + Lane_Km_Local_Total,
           
           Year = year)
  
  print(paste0("Finished combining ", year, " hpms data."))
  
  return(hpms_out)
  
}




# Function to return hpms df with full year coverage

full_coverage <- function(input_df, study_area, number_of_periods, ...) {
  
  cols <- enquos(...)
  
  coverage_df <- input_df %>%
    mutate(VKT_Urban_Class1 = ifelse(VKT_Urban_Class1 > 0, 1, 0),
           VKT_Urban_Class2 = ifelse(VKT_Urban_Class2 > 0, 1, 0),
           VKT_Urban_Class3 = ifelse(VKT_Urban_Class3 > 0, 1, 0),
           VKT_Urban_Class4 = ifelse(VKT_Urban_Class4 > 0, 1, 0),
           VKT_Urban_Local = ifelse(VKT_Urban_Local > 0, 1, 0),
           VKT_NonUrban_Class1 = ifelse(VKT_NonUrban_Class1 > 0, 1, 0),
           VKT_NonUrban_Class2 = ifelse(VKT_NonUrban_Class2 > 0, 1, 0),
           VKT_NonUrban_Class3 = ifelse(VKT_NonUrban_Class3 > 0, 1, 0),
           VKT_NonUrban_Class4 = ifelse(VKT_NonUrban_Class4 > 0, 1, 0),
           VKT_NonUrban_Local = ifelse(VKT_NonUrban_Local > 0, 1, 0),
           VKT_Class1_Total = ifelse(VKT_Class1_Total > 0, 1, 0),
           VKT_Class2_Total = ifelse(VKT_Class2_Total > 0, 1, 0),
           VKT_Class3_Total = ifelse(VKT_Class3_Total > 0, 1, 0),
           VKT_Class4_Total = ifelse(VKT_Class4_Total > 0, 1, 0),
           VKT_Local_Total = ifelse(VKT_Local_Total > 0, 1, 0),
           VKT_Urban_Total = ifelse(VKT_Urban_Total > 0, 1, 0),
           VKT_NonUrban_Total = ifelse(VKT_NonUrban_Total > 0, 1, 0),
           VKT_Total = ifelse(VKT_Total > 0, 1, 0),
           Lane_Km_Urban_Class1 = ifelse(Lane_Km_Urban_Class1 > 0, 1, 0),
           Lane_Km_Urban_Class2 = ifelse(Lane_Km_Urban_Class2 > 0, 1, 0),
           Lane_Km_Urban_Class3 = ifelse(Lane_Km_Urban_Class3 > 0, 1, 0),
           Lane_Km_Urban_Class4 = ifelse(Lane_Km_Urban_Class4 > 0, 1, 0),
           Lane_Km_Urban_Local = ifelse(Lane_Km_Urban_Local > 0, 1, 0),
           Lane_Km_NonUrban_Class1 = ifelse(Lane_Km_NonUrban_Class1 > 0, 1, 0),
           Lane_Km_NonUrban_Class2 = ifelse(Lane_Km_NonUrban_Class2 > 0, 1, 0),
           Lane_Km_NonUrban_Class3 = ifelse(Lane_Km_NonUrban_Class3 > 0, 1, 0),
           Lane_Km_NonUrban_Class4 = ifelse(Lane_Km_NonUrban_Class4 > 0, 1, 0),
           Lane_Km_NonUrban_Local = ifelse(Lane_Km_NonUrban_Local > 0, 1, 0),
           Lane_Km_Class1_Total = ifelse(Lane_Km_Class1_Total > 0, 1, 0),
           Lane_Km_Class2_Total = ifelse(Lane_Km_Class2_Total > 0, 1, 0),
           Lane_Km_Class3_Total = ifelse(Lane_Km_Class3_Total > 0, 1, 0),
           Lane_Km_Class4_Total = ifelse(Lane_Km_Class4_Total > 0, 1, 0),
           Lane_Km_Local_Total = ifelse(Lane_Km_Local_Total > 0, 1, 0),
           Lane_Km_Urban_Total = ifelse(Lane_Km_Urban_Total > 0, 1, 0),
           Lane_Km_NonUrban_Total = ifelse(Lane_Km_NonUrban_Total > 0, 1, 0),
           Lane_Km_Total = ifelse(Lane_Km_Total > 0, 1, 0),
           Total_Population = if_else(Total_Population > 0, 1, 0)) %>%
    group_by({{study_area}}) %>%
    summarise(VKT_Urban_Class1 = sum(VKT_Urban_Class1, na.rm = T),
              VKT_Urban_Class2 = sum(VKT_Urban_Class2, na.rm = T),
              VKT_Urban_Class3 = sum(VKT_Urban_Class3, na.rm = T),
              VKT_Urban_Class4 = sum(VKT_Urban_Class4, na.rm = T),
              VKT_Urban_Local = sum(VKT_Urban_Local, na.rm = T),
              VKT_NonUrban_Class1 = sum(VKT_NonUrban_Class1, na.rm = T),
              VKT_NonUrban_Class2 = sum(VKT_NonUrban_Class2, na.rm = T),
              VKT_NonUrban_Class3 = sum(VKT_NonUrban_Class3, na.rm = T),
              VKT_NonUrban_Class4 = sum(VKT_NonUrban_Class4, na.rm = T),
              VKT_NonUrban_Local = sum(VKT_NonUrban_Local, na.rm = T),
              VKT_Class1_Total = sum(VKT_Class1_Total, na.rm = T),
              VKT_Class2_Total = sum(VKT_Class2_Total, na.rm = T),
              VKT_Class3_Total = sum(VKT_Class3_Total, na.rm = T),
              VKT_Class4_Total = sum(VKT_Class4_Total, na.rm = T),
              VKT_Local_Total = sum(VKT_Local_Total, na.rm = T),
              VKT_Urban_Total = sum(VKT_Urban_Total, na.rm = T),
              VKT_NonUrban_Total = sum(VKT_NonUrban_Total , na.rm = T),
              VKT_Total = sum(VKT_Total, na.rm = T),
              
              Lane_Km_Urban_Class1 = sum(Lane_Km_Urban_Class1, na.rm = T),
              Lane_Km_Urban_Class2 = sum(Lane_Km_Urban_Class2 , na.rm = T),
              Lane_Km_Urban_Class3 = sum(Lane_Km_Urban_Class3, na.rm = T),
              Lane_Km_Urban_Class4 = sum(Lane_Km_Urban_Class4, na.rm = T),
              Lane_Km_Urban_Local = sum(Lane_Km_Urban_Local, na.rm = T),
              Lane_Km_NonUrban_Class1 = sum(Lane_Km_NonUrban_Class1 , na.rm = T),
              Lane_Km_NonUrban_Class2 = sum(Lane_Km_NonUrban_Class2 , na.rm = T),
              Lane_Km_NonUrban_Class3 = sum(Lane_Km_NonUrban_Class3, na.rm = T),
              Lane_Km_NonUrban_Class4 = sum(Lane_Km_NonUrban_Class4, na.rm = T),
              Lane_Km_NonUrban_Local = sum(Lane_Km_NonUrban_Local , na.rm = T),
              Lane_Km_Class1_Total = sum(Lane_Km_Class1_Total , na.rm = T),
              Lane_Km_Class2_Total = sum(Lane_Km_Class2_Total , na.rm = T),
              Lane_Km_Class3_Total = sum(Lane_Km_Class3_Total , na.rm = T),
              Lane_Km_Class4_Total = sum(Lane_Km_Class4_Total , na.rm = T),
              Lane_Km_Local_Total = sum(Lane_Km_Local_Total , na.rm = T),
              Lane_Km_Urban_Total = sum(Lane_Km_Urban_Total , na.rm = T),
              Lane_Km_NonUrban_Total = sum(Lane_Km_NonUrban_Total , na.rm = T),
              Lane_Km_Total = sum(Lane_Km_Total , na.rm = T),
              Total_Population = sum(Total_Population, na.rm = T)) %>%
    ungroup() %>%
    filter(if_all(c(!!!cols), ~ .x == number_of_periods)) %>%
    mutate(coverage_ind = "Yes") %>%
    select({{study_area}}, coverage_ind)
  
  coverage_df_out <- left_join(input_df, coverage_df, by = join_by({{study_area}})) %>%
    filter(coverage_ind == "Yes") %>%
    select(-coverage_ind)
  
  return(coverage_df_out)
  
}

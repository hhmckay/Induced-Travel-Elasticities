### VKT Elasticities - Analysis

# Run Function.R prior to running this script

################################################################################
### Data Cleaning
################################################################################

# Clean HPMS data for 1983, 1993, and 2003
hpms_1983 <- clean_hpms("HPMS/HPMS_1983.csv",
                        "State.FIPS.Code",
                        "County.FIPS.Code",
                        "Section.Identification",
                        "Rural.Urban",
                        "Functional.Class",
                        "Section.Grouped.Length",
                        "AADT",
                        "Number.of.Through.Lanes",
                        "HPMS/HPMS_SAMPLE_1983.csv",
                        "STATE.FIPS.CODE",
                        "COUNTY.FIPS.CODE",
                        "SECTION.IDENTIFICATION",
                        "RURAL.URBAN",
                        "FUNCTIONAL.CLASS",
                        "SECTION.GROUPED.LENGTH",
                        "AADT",
                        "NUMBER.OF.THROUGH.LANES",
                        "EXPANSION.FACTOR",
                        1983)

hpms_1993 <- clean_hpms("HPMS/HPMS_1993.csv",
                        "STATE",
                        "COUNTY",
                        "SECTION.ID",
                        "RURAL.URBAN",
                        "FUNCTIONAL.SYSTEM",
                        "SECTION.LENGTH",
                        "AADT",
                        "NUMBER.OF.LANES",
                        "HPMS/HPMS_SAMPLE_1993.csv",
                        "STATE",
                        "COUNTY",
                        "SECTION.ID",
                        "RURAL.URBAN",
                        "FUNCTIONAL.SYSTEM",
                        "SECTION.LENGTH",
                        "AADT",
                        "NUMBER.LANES",
                        "STANDARD.EXPANSION.FACTOR",
                        1993)

hpms_2003 <- clean_hpms("HPMS/HPMS_2003.csv",
                        "State_Code",
                        "County_Code",
                        "Section_ID",
                        "Rural.Urban",
                        "F_System",
                        "Section__Length",
                        "AADT",
                        "Through_Lanes",
                        "HPMS/HPMS_SAMPLE_2003.csv",
                        "STATE_CODE",
                        "COUNTY_CODE",
                        "SECTION_ID",
                        "RURAL_URBAN",
                        "F_SYSTEM",
                        "SECTION_LENGTH",
                        "AADT",
                        "THROUGH_LANES",
                        "STD_EXP_FACTOR",
                        2003)

### Read county population data
population <- read.csv("Population/POPULATION.csv") %>%
  mutate(County_ID = paste0(STATEFP, "_", COUNTYFP)) %>%
  rename(Total_Population = A00AA,
         Year = YEAR) %>%
  select(County_ID, Total_Population, Year) %>%
  mutate(Year = case_when(
    Year == 1980 ~ 1983,
    Year == 1990 ~ 1993,
    Year == 2000 ~ 2003,
    TRUE ~ NA
  )) %>%
  mutate(county_year_id = paste0(County_ID, Year)) %>%
  select(county_year_id, Total_Population)

# Read msa data
msa <- read.csv("MSAs/MSAs_2003.csv") %>%
  mutate(County_ID = paste0(State_Code, "_", County_Code)) %>%
  select(County_ID, MSA_Name)

# Combine HPMS FHs
hpms_df <- rbind(hpms_1983, hpms_1993, hpms_2003) %>%
  mutate(county_year_id = paste0(County_ID, Year)) %>%
  left_join(population, by = "county_year_id") %>%
  left_join(msa, by = "County_ID")

hpms_df_county <- hpms_df


################################################################################
### Models
################################################################################

### Table 1: MSA Elasticities by Functional Class

# Model 1: Class 1
hpms_df_msa_model <- full_coverage(hpms_df, County_ID, 3, VKT_Class1_Total, Lane_Km_Class1_Total, Total_Population) %>%
  group_by(MSA_Name, Year) %>%
  summarise(across(matches("^(VKT|Lane_Km|Total_Population)"), sum)) %>%
  ungroup() %>%
  filter(!is.na(MSA_Name))

df_panel <- pdata.frame(hpms_df_msa_model, index = c("MSA_Name", "Year"))

model1 <- run_model(df_panel, "log(VKT_Class1_Total + 1)", c("log(Lane_Km_Class1_Total + 1)", "log(Total_Population + 1)"),
                    "MSA_Name", "Year")

# Model 2: Class 2
hpms_df_msa_model <- full_coverage(hpms_df, County_ID, 3, VKT_Class2_Total, Lane_Km_Class2_Total, Total_Population) %>%
  group_by(MSA_Name, Year) %>%
  summarise(across(matches("^(VKT|Lane_Km|Total_Population)"), sum)) %>%
  ungroup() %>%
  filter(!is.na(MSA_Name))

df_panel <- pdata.frame(hpms_df_msa_model, index = c("MSA_Name", "Year"))

model2 <- run_model(df_panel, "log(VKT_Class2_Total + 1)", c("log(Lane_Km_Class2_Total + 1)", "log(Total_Population + 1)"),
                    "MSA_Name", "Year")

# Model 3: Class 3
hpms_df_msa_model <- full_coverage(hpms_df, County_ID, 3, VKT_Class3_Total, Lane_Km_Class3_Total, Total_Population) %>%
  group_by(MSA_Name, Year) %>%
  summarise(across(matches("^(VKT|Lane_Km|Total_Population)"), sum)) %>%
  ungroup() %>%
  filter(!is.na(MSA_Name))

df_panel <- pdata.frame(hpms_df_msa_model, index = c("MSA_Name", "Year"))

model3 <- run_model(df_panel, "log(VKT_Class3_Total + 1)", c("log(Lane_Km_Class3_Total + 1)", "log(Total_Population + 1)"),
                    "MSA_Name", "Year")

# Model 4: Class 1 - 3
hpms_df_msa_model <- full_coverage(hpms_df, County_ID, 3, VKT_Class3_Total, Lane_Km_Class3_Total, Total_Population) %>%
  group_by(MSA_Name, Year) %>%
  summarise(across(matches("^(VKT|Lane_Km|Total_Population)"), sum)) %>%
  ungroup() %>%
  filter(!is.na(MSA_Name)) %>%
  mutate(VKT_Class1_3_Total = VKT_Class1_Total + VKT_Class2_Total + VKT_Class3_Total,
         Lane_Km_Class1_3_Total = Lane_Km_Class1_Total + Lane_Km_Class2_Total + Lane_Km_Class3_Total)

df_panel <- pdata.frame(hpms_df_msa_model, index = c("MSA_Name", "Year"))

model4 <- run_model(df_panel, "log(VKT_Class1_3_Total + 1)", c("log(Lane_Km_Class1_3_Total + 1)", "log(Total_Population + 1)"),
                    "MSA_Name", "Year")

models <- list(model1, model2, model3, model4)

se_list <- lapply(models, function(mod) {
  se_vec <- mod$se
  names(se_vec) <- names(coef(mod))
  se_vec[names(coef(mod))]
})

se_p_list <- lapply(models, function(mod) {
  se_p_vec <- mod$se_p
  names(se_p_vec) <- names(coef(mod))
  se_p_vec[names(coef(mod))]
})

stargazer(models,
          type = "html",
          title = "Table 1: U.S. MSA VKT Elasticities by Roadway Class",
          se = se_list,
          p = se_p_list,
          dep.var.labels = c("VKT Class 1", "VKT Class 2", "VKT Class 3", "VKT Class 1 - 3"),
          covariate.labels = c("Lane Km Class 1", "Lane Km Class 2", "Lane Km Class 3", "Lane Km Class 1 - 3", "Population"),
          omit.stat = "all",
          add.lines = list(
            c("Observations", c(nobs(model1), nobs(model2), nobs(model3), nobs(model4))),
            c("Adjusted R2", c(model1$r2_adj, model2$r2_adj, model3$r2_adj, model4$r2_adj)),
            c("Adjusted R2 Within", c(model1$r2_within_adj, model2$r2_within_adj, model3$r2_within_adj, model4$r2_within_adj)),
            c("Robust F Statistic", c(model1$wald, model2$wald, model3$wald, model4$wald))))



### Table 2: CA MSA Elasticities by Functional Class

hpms_df_ca <- hpms_df %>%
  mutate(State_ID = sub("_.*", "", County_ID)) %>%
  filter(State_ID == "6")

# Model 1: Class 1
hpms_df_msa_model <- full_coverage(hpms_df_ca, County_ID, 3, VKT_Class1_Total, Lane_Km_Class1_Total, Total_Population) %>%
  group_by(MSA_Name, Year) %>%
  summarise(across(matches("^(VKT|Lane_Km|Total_Population)"), sum)) %>%
  ungroup() %>%
  filter(!is.na(MSA_Name))

df_panel <- pdata.frame(hpms_df_msa_model, index = c("MSA_Name", "Year"))

model1 <- run_model(df_panel, "log(VKT_Class1_Total + 1)", c("log(Lane_Km_Class1_Total + 1)", "log(Total_Population + 1)"),
                    "MSA_Name", "Year")

# Model 2: Class 2
hpms_df_msa_model <- full_coverage(hpms_df_ca, County_ID, 3, VKT_Class2_Total, Lane_Km_Class2_Total, Total_Population) %>%
  group_by(MSA_Name, Year) %>%
  summarise(across(matches("^(VKT|Lane_Km|Total_Population)"), sum)) %>%
  ungroup() %>%
  filter(!is.na(MSA_Name))

df_panel <- pdata.frame(hpms_df_msa_model, index = c("MSA_Name", "Year"))

model2 <- run_model(df_panel, "log(VKT_Class2_Total + 1)", c("log(Lane_Km_Class2_Total + 1)", "log(Total_Population + 1)"),
                    "MSA_Name", "Year")

# Model 3: Class 3
hpms_df_msa_model <- full_coverage(hpms_df_ca, County_ID, 3, VKT_Class3_Total, Lane_Km_Class3_Total, Total_Population) %>%
  group_by(MSA_Name, Year) %>%
  summarise(across(matches("^(VKT|Lane_Km|Total_Population)"), sum)) %>%
  ungroup() %>%
  filter(!is.na(MSA_Name))

df_panel <- pdata.frame(hpms_df_msa_model, index = c("MSA_Name", "Year"))

model3 <- run_model(df_panel, "log(VKT_Class3_Total + 1)", c("log(Lane_Km_Class3_Total + 1)", "log(Total_Population + 1)"),
                    "MSA_Name", "Year")

# Model 4: Class 1 - 3
hpms_df_msa_model <- full_coverage(hpms_df_ca, County_ID, 3, VKT_Class3_Total, Lane_Km_Class3_Total, Total_Population) %>%
  group_by(MSA_Name, Year) %>%
  summarise(across(matches("^(VKT|Lane_Km|Total_Population)"), sum)) %>%
  ungroup() %>%
  filter(!is.na(MSA_Name)) %>%
  mutate(VKT_Class1_3_Total = VKT_Class1_Total + VKT_Class2_Total + VKT_Class3_Total,
         Lane_Km_Class1_3_Total = Lane_Km_Class1_Total + Lane_Km_Class2_Total + Lane_Km_Class3_Total)

df_panel <- pdata.frame(hpms_df_msa_model, index = c("MSA_Name", "Year"))

model4 <- run_model(df_panel, "log(VKT_Class1_3_Total + 1)", c("log(Lane_Km_Class1_3_Total + 1)", "log(Total_Population + 1)"),
                    "MSA_Name", "Year")

models <- list(model1, model2, model3, model4)

se_list <- lapply(models, function(mod) {
  se_vec <- mod$se
  names(se_vec) <- names(coef(mod))
  se_vec[names(coef(mod))]
})

se_p_list <- lapply(models, function(mod) {
  se_p_vec <- mod$se_p
  names(se_p_vec) <- names(coef(mod))
  se_p_vec[names(coef(mod))]
})

stargazer(models,
          type = "html",
          title = "Table 2: CA MSA VKT Elasticities by Roadway Class",
          se = se_list,
          p = se_p_list,
          dep.var.labels = c("VKT Class 1", "VKT Class 2", "VKT Class 3", "VKT Class 1 - 3"),
          covariate.labels = c("Lane Km Class 1", "Lane Km Class 2", "Lane Km Class 3", "Lane Km Class 1 - 3", "Population"),
          omit.stat = "all",
          add.lines = list(
            c("Observations", c(nobs(model1), nobs(model2), nobs(model3), nobs(model4))),
            c("Adjusted R2", c(model1$r2_adj, model2$r2_adj, model3$r2_adj, model4$r2_adj)),
            c("Adjusted R2 Within", c(model1$r2_within_adj, model2$r2_within_adj, model3$r2_within_adj, model4$r2_within_adj)),
            c("Robust F Statistic", c(model1$wald, model2$wald, model3$wald, model4$wald))))



### Table 3: County Elasticities by Functional Class

# Model 1: Class 1
hpms_df_county_model <- full_coverage(hpms_df, County_ID, 3, VKT_Class1_Total, Lane_Km_Class1_Total, Total_Population) %>%
  filter(!is.na(County_ID))

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model1 <- run_model(df_panel, "log(VKT_Class1_Total + 1)", c("log(Lane_Km_Class1_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

# Model 2: Class 2
hpms_df_county_model <- full_coverage(hpms_df, County_ID, 3, VKT_Class2_Total, Lane_Km_Class2_Total, Total_Population) %>%
  filter(!is.na(County_ID))

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model2 <- run_model(df_panel, "log(VKT_Class2_Total + 1)", c("log(Lane_Km_Class2_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

# Model 3: Class 3
hpms_df_county_model <- full_coverage(hpms_df, County_ID, 3, VKT_Class3_Total, Lane_Km_Class3_Total, Total_Population) %>%
  filter(!is.na(County_ID))

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model3 <- run_model(df_panel, "log(VKT_Class3_Total + 1)", c("log(Lane_Km_Class3_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

# Model 4: Class 1 - 3
hpms_df_county_model <- full_coverage(hpms_df, County_ID, 3, VKT_Class3_Total, Lane_Km_Class3_Total, Total_Population) %>%
  filter(!is.na(County_ID)) %>%
  mutate(VKT_Class1_3_Total = VKT_Class1_Total + VKT_Class2_Total + VKT_Class3_Total,
         Lane_Km_Class1_3_Total = Lane_Km_Class1_Total + Lane_Km_Class2_Total + Lane_Km_Class3_Total)

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model4 <- run_model(df_panel, "log(VKT_Class1_3_Total + 1)", c("log(Lane_Km_Class1_3_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

models <- list(model1, model2, model3, model4)

se_list <- lapply(models, function(mod) {
  se_vec <- mod$se
  names(se_vec) <- names(coef(mod))
  se_vec[names(coef(mod))]
})

se_p_list <- lapply(models, function(mod) {
  se_p_vec <- mod$se_p
  names(se_p_vec) <- names(coef(mod))
  se_p_vec[names(coef(mod))]
})

stargazer(models,
          type = "html",
          title = "Table 3: U.S. County VKT Elasticities by Roadway Class",
          se = se_list,
          p = se_p_list,
          dep.var.labels = c("VKT Class 1", "VKT Class 2", "VKT Class 3", "VKT Class 1 - 3"),
          covariate.labels = c("Lane Km Class 1", "Lane Km Class 2", "Lane Km Class 3", "Lane Km Class 1 - 3", "Population"),
          omit.stat = "all",
          add.lines = list(
            c("Observations", c(nobs(model1), nobs(model2), nobs(model3), nobs(model4))),
            c("Adjusted R2", c(model1$r2_adj, model2$r2_adj, model3$r2_adj, model4$r2_adj)),
            c("Adjusted R2 Within", c(model1$r2_within_adj, model2$r2_within_adj, model3$r2_within_adj, model4$r2_within_adj)),
            c("Robust F Statistic", c(model1$wald, model2$wald, model3$wald, model4$wald))))


### Table 4: CA County Elasticities by Functional Class

hpms_df_ca <- hpms_df %>%
  mutate(State_ID = sub("_.*", "", County_ID)) %>%
  filter(State_ID == "6")

# Model 1: Class 1
hpms_df_county_model <- full_coverage(hpms_df_ca, County_ID, 3, VKT_Class1_Total, Lane_Km_Class1_Total, Total_Population) %>%
  filter(!is.na(County_ID))

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model1 <- run_model(df_panel, "log(VKT_Class1_Total + 1)", c("log(Lane_Km_Class1_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

# Model 2: Class 2
hpms_df_county_model <- full_coverage(hpms_df_ca, County_ID, 3, VKT_Class2_Total, Lane_Km_Class2_Total, Total_Population) %>%
  filter(!is.na(County_ID))

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model2 <- run_model(df_panel, "log(VKT_Class2_Total + 1)", c("log(Lane_Km_Class2_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

# Model 3: Class 3
hpms_df_county_model <- full_coverage(hpms_df_ca, County_ID, 3, VKT_Class3_Total, Lane_Km_Class3_Total, Total_Population) %>%
  filter(!is.na(County_ID))

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model3 <- run_model(df_panel, "log(VKT_Class3_Total + 1)", c("log(Lane_Km_Class3_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

# Model 4: Class 1 - 3
hpms_df_county_model <- full_coverage(hpms_df_ca, County_ID, 3, VKT_Class3_Total, Lane_Km_Class3_Total, Total_Population) %>%
  filter(!is.na(County_ID)) %>%
  mutate(VKT_Class1_3_Total = VKT_Class1_Total + VKT_Class2_Total + VKT_Class3_Total,
         Lane_Km_Class1_3_Total = Lane_Km_Class1_Total + Lane_Km_Class2_Total + Lane_Km_Class3_Total)

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model4 <- run_model(df_panel, "log(VKT_Class1_3_Total + 1)", c("log(Lane_Km_Class1_3_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

models <- list(model1, model2, model3, model4)

se_list <- lapply(models, function(mod) {
  se_vec <- mod$se
  names(se_vec) <- names(coef(mod))
  se_vec[names(coef(mod))]
})

se_p_list <- lapply(models, function(mod) {
  se_p_vec <- mod$se_p
  names(se_p_vec) <- names(coef(mod))
  se_p_vec[names(coef(mod))]
})

stargazer(models,
          type = "html",
          title = "Table 4: CA County VKT Elasticities by Roadway Class",
          se = se_list,
          p = se_p_list,
          dep.var.labels = c("VKT Class 1", "VKT Class 2", "VKT Class 3", "VKT Class 1 - 3"),
          covariate.labels = c("Lane Km Class 1", "Lane Km Class 2", "Lane Km Class 3", "Lane Km Class 1 - 3", "Population"),
          omit.stat = "all",
          add.lines = list(
            c("Observations", c(nobs(model1), nobs(model2), nobs(model3), nobs(model4))),
            c("Adjusted R2", c(model1$r2_adj, model2$r2_adj, model3$r2_adj, model4$r2_adj)),
            c("Adjusted R2 Within", c(model1$r2_within_adj, model2$r2_within_adj, model3$r2_within_adj, model4$r2_within_adj)),
            c("Robust F Statistic", c(model1$wald, model2$wald, model3$wald, model4$wald))))


### Table 5: Non-MSA County Elasticities by Functional Class

# Model 1: Class 1
hpms_df_county_model <- full_coverage(hpms_df, County_ID, 3, VKT_Class1_Total, Lane_Km_Class1_Total, Total_Population) %>%
  filter(!is.na(County_ID)) %>%
  filter(is.na(MSA_Name))

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model1 <- run_model(df_panel, "log(VKT_Class1_Total + 1)", c("log(Lane_Km_Class1_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

# Model 2: Class 2
hpms_df_county_model <- full_coverage(hpms_df, County_ID, 3, VKT_Class2_Total, Lane_Km_Class2_Total, Total_Population) %>%
  filter(!is.na(County_ID)) %>%
  filter(is.na(MSA_Name))

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model2 <- run_model(df_panel, "log(VKT_Class2_Total + 1)", c("log(Lane_Km_Class2_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

# Model 3: Class 3
hpms_df_county_model <- full_coverage(hpms_df, County_ID, 3, VKT_Class3_Total, Lane_Km_Class3_Total, Total_Population) %>%
  filter(!is.na(County_ID)) %>%
  filter(is.na(MSA_Name))

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model3 <- run_model(df_panel, "log(VKT_Class3_Total + 1)", c("log(Lane_Km_Class3_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

# Model 4: Class 1 - 3
hpms_df_county_model <- full_coverage(hpms_df, County_ID, 3, VKT_Class3_Total, Lane_Km_Class3_Total, Total_Population) %>%
  filter(!is.na(County_ID)) %>%
  filter(is.na(MSA_Name)) %>%
  mutate(VKT_Class1_3_Total = VKT_Class1_Total + VKT_Class2_Total + VKT_Class3_Total,
         Lane_Km_Class1_3_Total = Lane_Km_Class1_Total + Lane_Km_Class2_Total + Lane_Km_Class3_Total)

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model4 <- run_model(df_panel, "log(VKT_Class1_3_Total + 1)", c("log(Lane_Km_Class1_3_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

models <- list(model1, model2, model3, model4)

se_list <- lapply(models, function(mod) {
  se_vec <- mod$se
  names(se_vec) <- names(coef(mod))
  se_vec[names(coef(mod))]
})

se_p_list <- lapply(models, function(mod) {
  se_p_vec <- mod$se_p
  names(se_p_vec) <- names(coef(mod))
  se_p_vec[names(coef(mod))]
})

stargazer(models,
          type = "html",
          title = "Table 5: U.S. Non-MSA VKT Elasticities by Roadway Class",
          se = se_list,
          p = se_p_list,
          dep.var.labels = c("VKT Class 1", "VKT Class 2", "VKT Class 3", "VKT Class 1 - 3"),
          covariate.labels = c("Lane Km Class 1", "Lane Km Class 2", "Lane Km Class 3", "Lane Km Class 1 - 3", "Population"),
          omit.stat = "all",
          add.lines = list(
            c("Observations", c(nobs(model1), nobs(model2), nobs(model3), nobs(model4))),
            c("Adjusted R2", c(model1$r2_adj, model2$r2_adj, model3$r2_adj, model4$r2_adj)),
            c("Adjusted R2 Within", c(model1$r2_within_adj, model2$r2_within_adj, model3$r2_within_adj, model4$r2_within_adj)),
            c("Robust F Statistic", c(model1$wald, model2$wald, model3$wald, model4$wald))))



### Table 6: CA Non-MSA County Elasticities by Functional Class

hpms_df_ca <- hpms_df %>%
  mutate(State_ID = sub("_.*", "", County_ID)) %>%
  filter(State_ID == "6")

# Model 1: Class 1
hpms_df_county_model <- full_coverage(hpms_df_ca, County_ID, 3, VKT_Class1_Total, Lane_Km_Class1_Total, Total_Population) %>%
  filter(!is.na(County_ID)) %>%
  filter(is.na(MSA_Name))

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model1 <- run_model(df_panel, "log(VKT_Class1_Total + 1)", c("log(Lane_Km_Class1_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

# Model 2: Class 2
hpms_df_county_model <- full_coverage(hpms_df_ca, County_ID, 3, VKT_Class2_Total, Lane_Km_Class2_Total, Total_Population) %>%
  filter(!is.na(County_ID)) %>%
  filter(is.na(MSA_Name))

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model2 <- run_model(df_panel, "log(VKT_Class2_Total + 1)", c("log(Lane_Km_Class2_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

# Model 3: Class 3
hpms_df_county_model <- full_coverage(hpms_df_ca, County_ID, 3, VKT_Class3_Total, Lane_Km_Class3_Total, Total_Population) %>%
  filter(!is.na(County_ID)) %>%
  filter(is.na(MSA_Name))

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model3 <- run_model(df_panel, "log(VKT_Class3_Total + 1)", c("log(Lane_Km_Class3_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

# Model 4: Class 1 - 3
hpms_df_county_model <- full_coverage(hpms_df_ca, County_ID, 3, VKT_Class3_Total, Lane_Km_Class3_Total, Total_Population) %>%
  filter(!is.na(County_ID)) %>%
  filter(is.na(MSA_Name)) %>%
  mutate(VKT_Class1_3_Total = VKT_Class1_Total + VKT_Class2_Total + VKT_Class3_Total,
         Lane_Km_Class1_3_Total = Lane_Km_Class1_Total + Lane_Km_Class2_Total + Lane_Km_Class3_Total)

df_panel <- pdata.frame(hpms_df_county_model, index = c("County_ID", "Year"))

model4 <- run_model(df_panel, "log(VKT_Class1_3_Total + 1)", c("log(Lane_Km_Class1_3_Total + 1)", "log(Total_Population + 1)"),
                    "County_ID", "Year")

models <- list(model1, model2, model3, model4)

se_list <- lapply(models, function(mod) {
  se_vec <- mod$se
  names(se_vec) <- names(coef(mod))
  se_vec[names(coef(mod))]
})

se_p_list <- lapply(models, function(mod) {
  se_p_vec <- mod$se_p
  names(se_p_vec) <- names(coef(mod))
  se_p_vec[names(coef(mod))]
})

stargazer(models,
          type = "html",
          title = "Table 6: CA Non-MSA VKT Elasticities by Roadway Class",
          se = se_list,
          p = se_p_list,
          dep.var.labels = c("VKT Class 1", "VKT Class 2", "VKT Class 3", "VKT Class 1 - 3"),
          covariate.labels = c("Lane Km Class 1", "Lane Km Class 2", "Lane Km Class 3", "Lane Km Class 1 - 3", "Population"),
          omit.stat = "all",
          add.lines = list(
            c("Observations", c(nobs(model1), nobs(model2), nobs(model3), nobs(model4))),
            c("Adjusted R2", c(model1$r2_adj, model2$r2_adj, model3$r2_adj, model4$r2_adj)),
            c("Adjusted R2 Within", c(model1$r2_within_adj, model2$r2_within_adj, model3$r2_within_adj, model4$r2_within_adj)),
            c("Robust F Statistic", c(model1$wald, model2$wald, model3$wald, model4$wald))))

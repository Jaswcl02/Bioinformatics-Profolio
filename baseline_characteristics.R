#split metadata Sample_type -> sampletype and timepoint
meta_split <- do.call(rbind,strsplit(meta$Sample_type,"_"))
meta$Timepoint <- meta_split[,1]
meta$SampleType <- meta_split[,2]

#split sample_table Sample_type -> sampletype and timepoint
sampletable_split <- do.call(rbind,strsplit(dataset$sample_table$Sample_type,"_"))
dataset$sample_table$Timepoint <- sampletable_split[,1]
dataset$sample_table$SampleType <- sampletable_split[,2]

dataset$sample_table$GroupTime <- paste0(dataset$sample_table$Group, "_", dataset$sample_table$Timepoint)
dataset$sample_table$GroupSample <- paste0(dataset$sample_table$Group, "_", dataset$sample_table$SampleType)
dataset$meta$GroupSample <- paste0(dataset$meta$Group, "_", dataset$meta$SampleType)


#spilting dataset into baseline subset

meta_base <- subset(dataset$sample_table, 
                    grepl("W0", Timepoint, ignore.case = TRUE))

## Ensure OTU table matches new sample list
otu_base <- dataset$otu_table[, rownames(meta_base), drop = FALSE]

## Build  baseline dataset
dataset_base <- microtable$new(otu_table = otu_base,
                               tax_table = dataset$tax_table,
                               sample_table = meta_base)
dataset_base$tidy_dataset()
dataset_base


install.packages("tableone")


# Fix types
# Replace commas/spaces and coerce to numeric (safe if read as character)
numify <- function(x) as.numeric(gsub("[ ,]", "", as.character(x)))

meta_base <- meta_base %>%
  mutate(
    Age = numify(Age),
    BMI = numify(BMI),
    GroupSample = factor(Group)  # your 4 strata
  )

# Define variables
cont_vars <- c("Age", "BMI")  # continuous
cat_vars  <- c("Gender", "No_history_of_disease", "Allergy_History", "Asthma", 
               "Allergic_rhinitis", "Delivery_mode", "Feeding_mode", "Batch")  

# Keep only columns that exist in your data
cont_vars <- intersect(cont_vars, names(meta_base))
cat_vars  <- intersect(cat_vars,  names(meta_base))

# Make Table 1
tab1 <- CreateTableOne(vars   = c(cont_vars, cat_vars),
                       strata = "Group",   
                       data   = meta_base,  # Using original meta as in your code
                       test   = TRUE)

# Print with mean ± SD
print(tab1,
      showAllLevels = TRUE,
      quote = FALSE,
      noSpaces = TRUE,
      contDigits = 1)  # number of decimals

# Save to CSV
write.csv(print(tab1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE),
          "2708_baseline_characteristics.csv", row.names = TRUE)

# Alternative: Simple summary statistics for Age and BMI only
cat("\n=== Simple Age and BMI Statistics by Group ===\n")
age_bmi_summary <- meta_base %>%
  group_by(Group) %>%
  summarise(
    Age_mean = round(mean(Age, na.rm = TRUE), 1),
    Age_sd = round(sd(Age, na.rm = TRUE), 1),
    BMI_mean = round(mean(BMI, na.rm = TRUE), 1),
    BMI_sd = round(sd(BMI, na.rm = TRUE), 1),
    Age_mean_sd = paste0(Age_mean, " ± ", Age_sd),
    BMI_mean_sd = paste0(BMI_mean, " ± ", BMI_sd),
    n = n(),
    .groups = 'drop'
  ) %>%
  select(Group, Age_mean_sd, BMI_mean_sd, n)

print(age_bmi_summary)
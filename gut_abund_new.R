# Step 1: More thorough debugging of the NA issue
cat("=== DEBUGGING NA ISSUE ===\n")

# Check if NAs exist in Genus column after filtering
cat("NAs in Genus column after filtering:", sum(is.na(dataset_gut$tax_table$Genus)), "\n")
cat("Empty strings in Genus column:", sum(dataset_gut$tax_table$Genus == "", na.rm = TRUE), "\n")
cat("'NA' strings in Genus column:", sum(dataset_gut$tax_table$Genus == "NA", na.rm = TRUE), "\n")

# Check unique values in Genus column
unique_genera <- unique(dataset_gut$tax_table$Genus)
cat("First 10 unique genera:", paste(head(unique_genera, 10), collapse = ", "), "\n")
cat("Are there any actual NA values?", any(is.na(unique_genera)), "\n")

# Step 2: More aggressive filtering
# Remove any features where Genus is problematic
valid_genus_strict <- !is.na(dataset_gut$tax_table$Genus) & 
  dataset_gut$tax_table$Genus != "" & 
  dataset_gut$tax_table$Genus != "NA" &
  !grepl("^\\s*$", dataset_gut$tax_table$Genus) &  # Remove whitespace-only
  !grepl("unassigned|unclassified|unknown", dataset_gut$tax_table$Genus, ignore.case = TRUE)

cat("Features removed by strict filtering:", sum(!valid_genus_strict), "\n")
cat("Features remaining:", sum(valid_genus_strict), "\n")

# Apply strict filtering
dataset_gut$otu_table <- dataset_gut$otu_table[valid_genus_strict, ]
dataset_gut$tax_table <- dataset_gut$tax_table[valid_genus_strict, ]
dataset_gut$tax_table

# Step 3: Recreate the microtable object completely
dataset_gut_clean <- microtable$new(
  otu_table = dataset_gut$otu_table,
  tax_table = dataset_gut$tax_table,
  sample_table = dataset_gut$sample_table
)
dataset_gut_clean$tidy_dataset()


t1 <-  trans_abund$new(
  dataset = dataset_gut_clean,
  taxrank = "Genus",
  show = 0,
  ntaxa = 15,
  groupmean = NULL,
  group_morestats = FALSE,
  delete_taxonomy_lineage = TRUE,
  delete_taxonomy_prefix = TRUE,
  prefix = NULL,
  use_percentage = TRUE,
  input_taxaname = NULL,
  high_level = NULL,
  high_level_fix_nsub = NULL
)
dataset_gut_clean$sample_table
t1$plot_box(group ="GroupTime",xtext_angle = 60)

t1$data_abund 


# Choose genera you're interested in
genera_of_interest <- c("Faecalibacterium", "Fusobacterium")

# Filter for specific genera
selected_genera_data <- t1$data_abund %>%
  filter(Taxonomy %in% genera_of_interest, !is.na(Taxonomy))

# Display available genera (first 30)
cat("Available genera in dataset:\n")
available_genera <- unique(t1$data_abund$Taxonomy)
print(head(available_genera, 40))


# ===== STEP 2: VISUALIZATION OF RELATIVE ABUNDANCE =====

t1$data_abund <- selected_genera_data
View(t1$data_abund)
t1$plot_box(
    group = "GroupTime",
    xtext_angle = 45)

#=====

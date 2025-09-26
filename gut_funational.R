library(microeco)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(reshape2)

install.packages("file2meco", repos = BiocManager::repositories())
library("file2meco")
# Correct way to read your file
pathway_table <- read.delim("~/Downloads/OneDrive_1_25-08-2025/gut-picrust2_out_pipeline_split/pathways_out/path_abun_unstrat.tsv", 
                            row.names = 1, check.names = FALSE)

data("MetaCyc_pathway_map")


meta <- as.data.frame(meta_stool)

# Align rows of metadata to the columns (sample IDs) of your pathway table
if ("SampleID" %in% names(meta)) {
  meta <- meta %>% filter(SampleID %in% colnames(pathway_table)) %>%
    slice(match(colnames(pathway_table), SampleID))
  rownames(meta) <- meta$SampleID
} else {
  # If your SampleIDs are the rownames already:
  stopifnot(all(rownames(meta) %in% colnames(pathway_table)))
  meta <- meta[ match(colnames(pathway_table), rownames(meta)), , drop = FALSE]
}

stopifnot("GroupTime" %in% colnames(meta))  # <-- REQUIRED

## --- 2) Build microtable and compute abundances ---
tmp <- microtable$new(
  otu_table    = pathway_table,
  tax_table    = MetaCyc_pathway_map,
  sample_table = meta
)
tmp$tidy_dataset()
tmp$cal_abund()



# 1) Build the trans_abund object (returns a NEW object)
ta <- trans_abund$new(
  dataset = tmp,
  taxrank = "pathway",  # <- use "pathway" if you want individual pathways
  show = 0,
  ntaxa = 20,
  groupmean = "GroupTime",
  group_morestats = FALSE,
  delete_taxonomy_lineage = TRUE,
  delete_taxonomy_prefix = TRUE,
  prefix = NULL,
  use_percentage = TRUE
)
ta$data_abund
ta$data_abund$SE <-round(ta$data_abund$SE,4)
ta$data_abund$Abundance <-round(ta$data_abund$Abundance,4)
View(ta$data_abund)

write.csv(ta$data_abund,file = "2508_gut_pathway_abundance.csv",row.names = T)

t2 <- trans_diff$new(dataset = tmp,  # Use the full microtable object
                     method = 'lme', 
                     alpha = 0.05,
                     taxa_level = "pathway",  # pathways stored at Genus level
                     filter_thres = 0.001,
                     formula = "Group*Timepoint + (1|Subject_ID)",
                     p_adjust_method = "fdr")

View(t2$res_diff)
t2$res_diff
write.csv(t2$res_diff,file = "2508_gut_lme_pathway_2.csv",row.names = T)
#===ploting

library(dplyr)

t2$res_diff <- t2$res_diff[t2$res_diff$Factors != "c(Intercept,	
TimepointW12,GroupPlacebo)", ]

# now plot without intercept
t2$plot_diff_bar(
  heatmap_cell    = "Estimate",
  heatmap_sig     = "Significance",
  heatmap_lab_fill = "Coefficient"
)

#===ploting====
# Define the specific pathways you want to plot
target_pathways <- c(
  "superpathway of tetrahydrofolate biosynthesis and salvage",
  "6-hydroxymethyl-dihydropterin diphosphate biosynthesis III (Chlamydia)",
  "mixed acid fermentation",
  "superpathway of tetrahydrofolate biosynthesis",
  "pyrimidine deoxyribonucleotides de novo biosynthesis II",
  "superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis (E. coli)",
  "D-galacturonate degradation I",
  "sucrose degradation IV (sucrose phosphorylase)",
  "superpathway of N-acetylneuraminate degradation",
  "superpathway of glycolysis and the Entner-Doudoroff pathway"
)

# Remove duplicates (I noticed you have the E. coli pathway twice)
target_pathways <- unique(target_pathways)

print(paste("Target pathways:", length(target_pathways)))

# Filter the differential results to only include these pathways AND specific factors
target_factors <- c("Timepoint", "Group", "Group:Timepoint")

selected_data <- t2$res_diff %>%
  filter(Taxa %in% target_pathways) %>%
  filter(Factors %in% target_factors)

print(paste("Found", nrow(selected_data), "records for target pathways"))
print("Available pathways in data:")
print(unique(selected_data$Taxa))

# If some pathways are missing, let's try partial matching
if(nrow(selected_data) == 0) {
  print("No exact matches found. Trying partial matching...")
  
  # Try to find pathways that contain key terms
  key_terms <- c("tetrahydrofolate", "dihydropterin", "mixed acid", "pyrimidine", 
                 "galacturonate", "sucrose", "acetylneuraminate", "glycolysis", "Entner")
  
  selected_data <- t2$res_diff %>%
    filter(Factors %in% target_factors) %>%
    filter(sapply(Taxa, function(x) any(sapply(key_terms, function(term) grepl(term, x, ignore.case = TRUE)))))
  
  print(paste("Found", nrow(selected_data), "records with partial matching"))
  print("Found pathways:")
  print(unique(selected_data$Taxa))
}

# If we still have data, create the heatmap
if(nrow(selected_data) > 0) {
  
  # Check data structure
  print("Data structure:")
  print(selected_data[, c("Taxa", "Factors", "Estimate", "P.unadj", "Significance")])
  
  # Method 1: Simple ggplot2 heatmap (works with any data structure)
  library(ggplot2)
  
  # Create a cleaner factor name for plotting
  selected_data$Clean_Factor <- selected_data$Factors
  selected_data$Clean_Factor <- gsub(":", " × ", selected_data$Clean_Factor)  # Make interaction clearer
  
  # Create the heatmap
  p_heat <- ggplot(selected_data, aes(x = Clean_Factor, y = Taxa, fill = Estimate)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = Significance), size = 4, fontface = "bold") +
    scale_fill_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red", 
      midpoint = 0,
      name = "Coefficient"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    labs(
      title = "Differential Pathway Analysis",
      subtitle = "Selected Metabolic Pathways",
      x = "Treatment Contrast",
      y = "Metabolic Pathway"
    ) +
    # Wrap long pathway names
    scale_y_discrete(labels = function(x) {
      sapply(x, function(name) {
        if(nchar(name) > 50) {
          # Split at logical points
          words <- strsplit(name, " ")[[1]]
          if(length(words) > 6) {
            paste(c(paste(words[1:6], collapse = " "), 
                    paste(words[7:length(words)], collapse = " ")), 
                  collapse = "\n")
          } else {
            name
          }
        } else {
          name
        }
      })
    })
  
  # Display the plot
  print(p_heat)
  
  # Save the plot
  ggsave("selected_pathways_heatmap.png", p_heat, 
         width = 12, height = 8, dpi = 300, bg = "white")
  
  ggsave("selected_pathways_heatmap.pdf", p_heat, 
         width = 12, height = 8, bg = "white")
  
  print("Heatmap saved as 'selected_pathways_heatmap.png' and '.pdf'")
  
  # Method 2: If you want to use ComplexHeatmap and have matrix data
  if(length(unique(selected_data$Taxa)) > 1 && length(unique(selected_data$Factors)) > 1) {
    
    # Try ComplexHeatmap approach
    if(require(ComplexHeatmap, quietly = TRUE)) {
      library(circlize)
      
      # Create matrices
      heatmap_wide <- reshape(selected_data[, c("Taxa", "Factors", "Estimate")], 
                              idvar = "Taxa", 
                              timevar = "Factors", 
                              direction = "wide", 
                              v.names = "Estimate")
      
      rownames(heatmap_wide) <- heatmap_wide$Taxa
      heatmap_matrix <- as.matrix(heatmap_wide[, -1])
      heatmap_matrix[is.na(heatmap_matrix)] <- 0
      
      # Significance matrix
      sig_wide <- reshape(selected_data[, c("Taxa", "Factors", "Significance")], 
                          idvar = "Taxa", 
                          timevar = "Factors", 
                          direction = "wide", 
                          v.names = "Significance")
      
      rownames(sig_wide) <- sig_wide$Taxa
      sig_matrix <- as.matrix(sig_wide[, -1])
      sig_matrix[is.na(sig_matrix)] <- ""
      
      # Clean column names
      colnames(heatmap_matrix) <- gsub("Estimate\\.", "", colnames(heatmap_matrix))
      colnames(sig_matrix) <- gsub("Significance\\.", "", colnames(sig_matrix))
      
      # Color function
      col_fun <- colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)), 
                            c("blue", "white", "red"))
      
      # Create ComplexHeatmap
      ht <- Heatmap(
        heatmap_matrix,
        name = "Coefficient",
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 11),
        column_names_rot = 45,
        cluster_rows = FALSE,  # Keep pathway order
        cluster_columns = FALSE,
        # Add significance stars
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(sig_matrix[i, j] %in% c("*", "**", "***")) {
            grid.text(sig_matrix[i, j], x, y, gp = gpar(fontsize = 10, col = "black", fontface = "bold"))
          }
        },
        width = unit(8, "cm"),
        height = unit(12, "cm")
      )
      
      # Draw ComplexHeatmap
      pdf("selected_pathways_complexheatmap.pdf", width = 10, height = 8)
      draw(ht)
      dev.off()
      
      print("ComplexHeatmap also saved as 'selected_pathways_complexheatmap.pdf'")
    }
  }
  
} else {
  print("No data found for the specified pathways.")
  print("Available pathways in your data:")
  print(head(unique(t2$res_diff$Taxa), 20))
}

#======

# Assuming your sample_table has columns: Subject_ID, Treatment, Time

###spilting stool dataset into treatment subset 
meta_path_T_gut <- subset(tmp$sample_table, 
                             grepl("MHLJDD", Group, ignore.case = TRUE))

## Ensure OTU table matches new sample list
otu_path_T_gut <- tmp$otu_table[, rownames(meta_path_T_gut), drop = FALSE]

## 3) Build gut/stool dataset
dataset_path_T_gut <- microtable$new(otu_table = otu_path_T_gut,
                                        tax_table = MetaCyc_pathway_map,
                                        sample_table = meta_path_T_gut)
dataset_path_T_gut$tidy_dataset()



treat.obj <- trans_diff$new(dataset = dataset_path_T_gut,  # Use the full microtable object
                     method = 'wilcox', 
                     taxa_level = "pathway", 
                     group = "Timepoint",# pathways stored at Genus level
                     by_ID = "Subject_ID")
View(treat.obj$res_diff)

#======
# Assuming your sample_table has columns: Subject_ID, Treatment, Time

###spilting stool dataset into treatment subset 
meta_path_P_gut <- subset(tmp$sample_table, 
                          grepl("Placebo", Group, ignore.case = TRUE))

## Ensure OTU table matches new sample list
otu_path_P_gut <- tmp$otu_table[, rownames(meta_path_P_gut), drop = FALSE]

## 3) Build gut/stool dataset
dataset_path_P_gut <- microtable$new(otu_table = otu_path_P_gut,
                                     tax_table = MetaCyc_pathway_map,
                                     sample_table = meta_path_P_gut)
dataset_path_P_gut$tidy_dataset()



place.obj <- trans_diff$new(dataset = dataset_path_P_gut,  # Use the full microtable object
                            method = 'wilcox', 
                            taxa_level = "pathway", 
                            group = "Timepoint",# pathways stored at Genus level
                            by_ID = "Subject_ID")
View(place.obj$res_diff)
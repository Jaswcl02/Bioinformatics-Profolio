library(microeco)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(reshape2)

install.packages("file2meco", repos = BiocManager::repositories())
library("file2meco")
# Correct way to read your file
or_pathway_table <- read.delim("~/Downloads/OneDrive_1_25-08-2025/oral-picrust2_out_pipeline_split/pathways_out/path_abun_unstrat.tsv", 
                            row.names = 1, check.names = FALSE)

data("MetaCyc_pathway_map")


meta <- as.data.frame(meta_saliva)

# Align rows of metadata to the columns (sample IDs) of your pathway table
if ("SampleID" %in% names(meta)) {
  meta <- meta %>% filter(SampleID %in% colnames(or_pathway_table)) %>%
    slice(match(colnames(or_pathway_table), SampleID))
  rownames(meta) <- meta$SampleID
} else {
  # If your SampleIDs are the rownames already:
  stopifnot(all(rownames(meta) %in% colnames(or_pathway_table)))
  meta <- meta[ match(colnames(or_pathway_table), rownames(meta)), , drop = FALSE]
}

stopifnot("GroupTime" %in% colnames(meta))  # <-- REQUIRED

## --- 2) Build microtable and compute abundances ---
tmp_or <- microtable$new(
  otu_table    = or_pathway_table,
  tax_table    = MetaCyc_pathway_map,
  sample_table = meta
)
tmp_or$tidy_dataset()
tmp_or$cal_abund()



# 1) Build the trans_abund object (returns a NEW object)
to <- trans_abund$new(
  dataset = tmp_or,
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
to$data_abund
to$data_abund$SE <-round(to$data_abund$SE,4)
to$data_abund$Abundance <-round(to$data_abund$Abundance,4)
View(to$data_abund)

write.csv(to$data_abund,file = "2608_oral_pathway_abundance_2.csv",row.names = T)

t3 <- trans_diff$new(dataset = tmp_or,  # Use the full microtable object
                     method = 'lme', 
                     alpha = 0.05,
                     taxa_level = "pathway",  # pathways stored at Genus level
                     filter_thres = 0.001,
                     formula = "Group*Timepoint + (1|Subject_ID)",
                     p_adjust_method = "fdr")

View(t3$res_diff)
t3$res_diff
write.csv(t3$res_diff,file = "2608_oral_lme_pathway_2.csv",row.names = T)

#=====ploting

=# Required libraries
library(microeco)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)  # For combining plots

# Assuming your existing code up to t3$res_diff is already run...

# Function to create combined heatmap and boxplot visualization
create_pathway_visualization <- function(t3, tmp_or, top_n = 10) {
  
  # 1. Get top significant pathways
  top_pathways <- diff_results %>%
    filter(!is.na(Pr)) %>%
    arrange(Pr) %>%
    head(top_n) %>%
    pull(Taxa)
  
  # 2. Create heatmap data
  heatmap_data <- diff_results %>%
    filter(Taxa %in% top_pathways) %>%
    select(Taxa, Term, Estimate, Pr) %>%
    mutate(
      Significance = case_when(
        Pr < 0.001 ~ "***",
        Pr < 0.01 ~ "**", 
        Pr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      Taxa = factor(Taxa, levels = rev(top_pathways))  # Reverse for proper ordering
    )
  
  # 3. Create heatmap
  heatmap_plot <- ggplot(heatmap_data, aes(x = Term, y = Taxa)) +
    geom_tile(aes(fill = Estimate), color = "white", size = 0.5) +
    geom_text(aes(label = Significance), size = 4, fontface = "bold") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, name = "Coefficient"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 9),
      axis.title = element_blank(),
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    labs(title = "Top 10 Most Significant Pathway Changes")
  
  # 4. Prepare boxplot data
  # Extract abundance data for top pathways
  abundance_data <- microtable_obj$otu_table[top_pathways, , drop = FALSE]
  
  # Convert to relative abundance (%)
  abundance_data <- sweep(abundance_data, 2, colSums(microtable_obj$otu_table), "/") * 100
  
  # Reshape for plotting
  boxplot_data <- abundance_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Pathway") %>%
    reshape2::melt(id.vars = "Pathway", variable.name = "SampleID", value.name = "Abundance") %>%
    left_join(
      microtable_obj$sample_table %>% 
        tibble::rownames_to_column("SampleID") %>%
        select(SampleID, Group, Timepoint, GroupTime),
      by = "SampleID"
    ) %>%
    mutate(
      Pathway = factor(Pathway, levels = top_pathways),
      Group = factor(Group),
      Timepoint = factor(Timepoint)
    )
  
  # 5. Create boxplots for each term in your model
  create_boxplot <- function(grouping_var, title_suffix) {
    ggplot(boxplot_data, aes_string(x = grouping_var, y = "Abundance")) +
      geom_boxplot(aes_string(fill = grouping_var), alpha = 0.7, outlier.size = 1) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 0.8) +
      facet_wrap(~Pathway, scales = "free_y", nrow = 2) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8),
        legend.position = "bottom"
      ) +
      labs(
        y = "Relative Abundance (%)",
        title = paste("Pathway Abundance by", title_suffix)
      )
  }
  
  # Create individual boxplots
  group_boxplot <- create_boxplot("Group", "Treatment Group")
  timepoint_boxplot <- create_boxplot("Timepoint", "Timepoint")
  grouptime_boxplot <- create_boxplot("GroupTime", "Group × Timepoint")
  
  # 6. Combine all plots
  combined_plot <- heatmap_plot / 
    (group_boxplot | timepoint_boxplot) / 
    grouptime_boxplot +
    plot_layout(heights = c(1, 1.5, 1.5))
  
  return(list(
    combined_plot = combined_plot,
    heatmap = heatmap_plot,
    group_boxplot = group_boxplot,
    timepoint_boxplot = timepoint_boxplot,
    grouptime_boxplot = grouptime_boxplot,
    data = list(
      heatmap_data = heatmap_data,
      boxplot_data = boxplot_data
    )
  ))
}

# Create the visualization
pathway_viz <- create_pathway_visualization(t3$res_diff, tmp_or, top_n = 10)

# Display the combined plot
pathway_viz$combined_plot

# Save the plot
ggsave("pathway_lme_analysis.png", 
       pathway_viz$combined_plot, 
       width = 16, height = 12, 
       dpi = 300)

# Alternative: Create a more compact version similar to your example image
create_compact_visualization <- function(diff_results, microtable_obj, top_n = 10) {
  
  # Get top pathways
  top_pathways <- diff_results %>%
    filter(!is.na(Pr)) %>%
    arrange(Pr) %>%
    head(top_n) %>%
    pull(Taxa)
  
  # Heatmap data
  heatmap_data <- diff_results %>%
    filter(Taxa %in% top_pathways) %>%
    select(Taxa, Term, Estimate, Pr) %>%
    mutate(
      Significance = case_when(
        Pr < 0.001 ~ "***",
        Pr < 0.01 ~ "**", 
        Pr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      Taxa = factor(Taxa, levels = rev(top_pathways))
    )
  
  # Compact heatmap
  heatmap_compact <- ggplot(heatmap_data, aes(x = Term, y = Taxa)) +
    geom_tile(aes(fill = Estimate), color = "white") +
    geom_text(aes(label = Significance), size = 3, fontface = "bold") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#762A83",
      midpoint = 0, name = "Coefficient"
    ) +
    theme_void() +
    theme(
      axis.text.y = element_text(size = 8, hjust = 1),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      legend.key.height = unit(0.3, "cm")
    )
  
  # Boxplot data for the 4 main pathways that changed most
  top_4_pathways <- head(top_pathways, 4)
  
  abundance_data <- microtable_obj$otu_table[top_4_pathways, , drop = FALSE]
  abundance_data <- sweep(abundance_data, 2, colSums(microtable_obj$otu_table), "/") * 100
  
  boxplot_data <- abundance_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Pathway") %>%
    reshape2::melt(id.vars = "Pathway", variable.name = "SampleID", value.name = "Abundance") %>%
    left_join(
      microtable_obj$sample_table %>% 
        tibble::rownames_to_column("SampleID") %>%
        select(SampleID, Group, Timepoint),
      by = "SampleID"
    ) %>%
    mutate(
      Pathway = factor(Pathway, levels = top_4_pathways),
      Group = factor(Group),
      Timepoint = factor(Timepoint)
    )
  
  # Create compact boxplot
  boxplot_compact <- ggplot(boxplot_data, aes(x = Timepoint, y = Abundance)) +
    geom_boxplot(aes(fill = Group), alpha = 0.7, outlier.size = 0.5) +
    facet_wrap(~Pathway, scales = "free_y", nrow = 1) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 7),
      axis.text = element_text(size = 7),
      legend.position = "bottom"
    ) +
    labs(y = "Relative abundance (%)", x = "")
  
  # Combine
  combined_compact <- heatmap_compact / boxplot_compact +
    plot_layout(heights = c(2, 1))
  
  return(combined_compact)
}

# Create compact version
pathway_viz_compact <- create_compact_visualization(t3$res_diff, tmp_or, top_n = 10)
pathway_viz_compact

# Save compact version
ggsave("pathway_lme_compact.png", 
       pathway_viz_compact, 
       width = 12, height = 8, 
       dpi = 300)

# Print summary of top pathways
cat("Top 10 Most Significant Pathways:\n")
t3$res_diff %>%
  filter(!is.na(Pr)) %>%
  arrange(Pr) %>%
  head(10) %>%
  select(Taxa, Term, Estimate, Pr) %>%
  print()


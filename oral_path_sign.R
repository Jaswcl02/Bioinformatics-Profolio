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
View(to)

write.csv(to$data_abund,file = "2508_oral_pathway_abundance.csv",row.names = T)

t3 <- trans_diff$new(dataset = tmp_or,  # Use the full microtable object
                     method = 'lme', 
                     alpha = 0.05,
                     taxa_level = "pathway",  # pathways stored at Genus level
                     filter_thres = 0.001,
                     formula = "Group*Timepoint + (1|Subject_ID)",
                     p_adjust_method = "fdr")

View(t3$res_diff)
t3$res_diff
write.csv(t3$res_diff,file = "2508_oral_lme_pathway_2.csv",row.names = T)
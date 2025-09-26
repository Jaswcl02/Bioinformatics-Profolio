abund.obj <- trans_abund$new(
  dataset = dataset_saliva,
  taxrank = "Genus",
  show = 0,
  ntaxa = NULL,
  groupmean = "GroupTime",
  group_morestats = FALSE,
  delete_taxonomy_lineage = TRUE,
  delete_taxonomy_prefix = TRUE,
  prefix = NULL,
  use_percentage = TRUE,
  input_taxaname = NULL,
  high_level = NULL,
  high_level_fix_nsub = NULL
)



#=====
# factor hygiene
st$Group     <- factor(st$Group)
st$Timepoint <- factor(st$Timepoint, ordered = FALSE)
st$Batch     <- factor(st$Batch)
st$Subject_ID<- factor(st$Subject_ID)

# any empty cells in the Group×Timepoint design?
addmargins(table(st$Group, st$Timepoint))

# how Batch relates to Subject_ID (confounding check)
tab_SB <- table(st$Subject_ID, st$Batch)
summary(rowSums(tab_SB > 0))           # subjects across how many batches?
mean(rowSums(tab_SB > 0) == 1)         # proportion of subjects in exactly 1 batch

# is Batch basically unique per subject? (bad for crossed RE)
# If TRUE or near-TRUE, (1|Batch) + (1|Subject_ID) will be confounded/singular.

# any missing/constant covariate?
summary(st$BMI)
sd(as.numeric(st$BMI), na.rm = TRUE)
any(is.na(st$BMI))

tb <- table(st$Batch)
summary(tb);  mean(tb == 1)  # proportion of singleton batches

t1 <- trans_diff$new(dataset = dataset_saliva, method = 'lme', alpha = 0.05, 
                     taxa_level = "Genus",
                     filter_thres = 0.001,
                     formula = "Group*Timepoint+ (1 | Subject_ID) ")

View(t1$res_diff)

unique(t1$res_diff$term)

str(dataset_saliva$res_diff)

#===heatmap====
library(dplyr)

t1$res_diff <- t1$res_diff[t1$res_diff$Factors != "(Intercept)", ]

# now plot without intercept
t1$plot_diff_bar(
  heatmap_cell    = "Estimate",
  heatmap_sig     = "Significance",
  heatmap_lab_fill = "Coefficient"
)
#===line plot====
sel_genera <- c("Lautropia","Streptococcus","Corynebacterium","Mammalia")

# keep only those rows in the microeco object
abund.obj$data_abund <- abund.obj$data_abund %>%
  dplyr::filter(Taxonomy %in% sel_genera)

# now call the usual plot function
abund.obj$plot_line(position = position_dodge(0.3), xtext_angle = 0)


#====table plot 
library(dplyr)

sel_best_nonint <- sel_res %>%
  # remove intercepts and non-estimable rows
  filter(!Factors %in% c("(Intercept)", "(1 | Subject_ID)", "<none>", "Group", "Timepoint", "Group:Timepoint"),
         !is.na(P_unadj)) %>%
  # pick the lowest p-value per genus
  group_by(Taxon) %>%
  slice_min(order_by = P_unadj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Taxon)

sel_best_nonint

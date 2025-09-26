#===overall alpha indices====
dataset_saliva$alpha_diversity
dataset_saliva$alpha_diversity
alpha_dataset_saliva.obj <- trans_alpha$new(dataset = dataset_saliva, group = "GroupTime")
alpha_dataset_saliva.obj$cal_diff(method = "wilcox")
alpha_dataset_saliva.obj$data_stat


alpha_dataset_saliva.obj$plot_alpha(measure = "InvSimpson", add_sig = FALSE,  add = "point",
                                 point_size = 3,
                                 point_alpha = 0.8)

alpha_dataset_saliva.obj$plot_alpha(measure = "Shannon", add_sig = FALSE,  add = "point",
                                 point_size = 3,
                                 point_alpha = 0.8)


#paired alpha diversity for dataset_placebo_oral, comparing difference in W0 and w12
alpha_dataset_placebo_oral.obj <- trans_alpha$new(dataset = dataset_placebo_oral, group = "Timepoint",by_ID = 'Subject_ID')
alpha_dataset_placebo_oral.obj$cal_diff(method = "wilcox")

alpha_dataset_placebo_oral.obj$plot_alpha(measure = "Shannon", add_sig = TRUE,  add = "point",
                                         point_size = 3,
                                         point_alpha = 0.8)

head(dataset_placebo_oral$alpha_diversity)
write.csv(dataset_placebo_oral$alpha_diversity,file = "1408_dataset_placebo_oral_alpha_diversity_w0w12",row.names = TRUE)

alpha_dataset_placebo_oral.obj$res_diff
write.csv(dataset_placebo_oral$res_diff,file = "1408_dataset_placebo_oral_res_diff_w0w12",row.names = TRUE)

#paired alpha diversity for dataset_treatment_oral, comparing difference in W0 and w12
alpha_dataset_treatment_oral.obj <- trans_alpha$new(dataset = dataset_treatment_oral, group = "Timepoint",by_ID = 'Subject_ID')
alpha_dataset_treatment_oral.obj$cal_diff(method = "wilcox")
alpha_dataset_treatment_oral.obj$plot_alpha(measure = "Shannon", add_sig = TRUE,add = "point",
                                           point_size = 3,
                                           point_alpha = 0.8)
head(dataset_treatment_oral$alpha_diversity)
write.csv(dataset_treatment_gut$alpha_diversity,file = "1408_dataset_treatment_oral_alpha_diversity_w0w12",row.names = TRUE)
alpha_dataset_treatment_oral.obj$res_diff
write.csv(alpha_dataset_treatment_oral.obj$res_diff,file = "1408_dataset_treatment_oral_ares_diff_w0w12",row.names = TRUE)

#unpaired alpha diversity for dataset_w0_oral, comparing difference placebo and treatmenr in w0
alpha_dataset_w0_oral.obj <- trans_alpha$new(dataset = dataset_w0_oral, group = "Group")
alpha_dataset_w0_oral.obj$cal_diff(method = "wilcox")
alpha_dataset_w0_oral.obj$res_diff
alpha_dataset_w0_oral.obj$plot_alpha(measure = "Shannon", add_sig = TRUE,add = "point",
                                    point_size = 3,
                                    point_alpha = 0.8)
head(dataset_w0_oral$alpha_diversity)
write.csv(dataset_w0_oral$alpha_diversity,file = "dataset_w0_oral_alpha_diversity_PvsT",row.names = TRUE)

#unpaired alpha diversity for dataset_w12_oral, comparing difference placebo and treatmenr in w12
alpha_dataset_w12_oral.obj <- trans_alpha$new(dataset = dataset_w12_oral, group = "Group")
alpha_dataset_w12_oral.obj$cal_diff(method = "wilcox")
alpha_dataset_w12_oral.obj$res_diff 
alpha_dataset_w12_oral.obj$plot_alpha(measure = "Shannon", add_sig = TRUE,add = "point",
                                     point_size = 3,
                                     point_alpha = 0.8)
write.csv(dataset_w12_gut$alpha_diversity,file = "dataset_w12_oral_alpha_diversity_PvsT",row.names = TRUE)

###lme (method 1)
alpha_obj <- trans_alpha$new(dataset = dataset_saliva)


str(dataset_saliva$alpha_diversity)
rownames(dataset_saliva$alpha_diversity)
rownames(dataset_saliva$sample_table)

dataset_saliva$sample_table$Alpha_Diversity <- dataset_saliva$alpha_diversity$Shannon
str(dataset_saliva$sample_table$Alpha_Diversity)

library(lme4)

model <- lmer(Alpha_Diversity ~ Timepoint * Group + (1 | Subject_ID), data = dataset_saliva$sample_table)
summary(model)

library(broom.mixed)

fixed_effects <- tidy(model, effects = "fixed")
fixed_effects$Signif <- cut(
  fixed_effects$p.value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
  labels = c("***", "**", "*", ".", "")
)
library(knitr)

kable(
  fixed_effects[, c("term", "estimate", "std.error", "df", "statistic", "p.value", "Signif")],
  digits = 3,
  caption = "Fixed Effects from LME Model"
)

write.csv(
  fixed_effects,
  file = "saliva_lme_fixed_effects.csv",
  row.names = FALSE
)

#===lme table (working method)
## ---------------------------
## 1) Merge sample + alpha data

## ---------------------------
dalpha_dataset_oral.obj <- trans_alpha$new(dataset = dataset_saliva, group = "Timepoint")

sample_data <- dataset_saliva$sample_table
alpha_data  <- dataset_saliva$alpha_diversity

stopifnot(!is.null(rownames(sample_data)), !is.null(rownames(alpha_data)))

common_ids  <- intersect(rownames(sample_data), rownames(alpha_data))
if (length(common_ids) == 0) stop("No overlapping sample IDs between sample_table and alpha_diversity.")

# align by common sample IDs
sample_data <- sample_data[common_ids, , drop = FALSE]
alpha_data  <- alpha_data[common_ids, , drop = FALSE]

merged_data <- cbind(sample_data, alpha_data)

## ---------------------------
## 2) Model function (crossed RE: Subject_ID + Batch)
## ---------------------------
library(lme4)
library(lmerTest)

run_lmer_analysis <- function(data, response_var) {
  # ensure factors
  data$Group      <- as.factor(data$Group)
  data$Timepoint  <- as.factor(data$Timepoint)
  data$Subject_ID <- as.factor(data$Subject_ID)
  data$Batch      <- as.factor(data$Batch)
  
  # build formula: fixed + crossed random intercepts
  form <- as.formula(paste(response_var, "~ Group * Timepoint + (1|Subject_ID) + (1|Batch)"))
  
  # fit model
  fit <- lmer(form, data = data, REML = TRUE, na.action = na.omit)
  
  # Type III ANOVA with Satterthwaite df
  aov_tab <- anova(fit, type = 3)
  list(model = fit, anova = aov_tab, summary = summary(fit))
}

## ---------------------------
## 3) Run across metrics
## ---------------------------
# pick alpha metrics from alpha_data columns (auto-detect)
alpha_metrics <- colnames(alpha_data)

results_list <- list()
for (metric in alpha_metrics) {
  cat("\n=== Analysis for", metric, "===\n")
  # skip non-numeric columns defensively
  if (!is.numeric(merged_data[[metric]])) {
    warning(sprintf("Metric '%s' is not numeric; skipping.", metric))
    next
  }
  out <- try(run_lmer_analysis(merged_data, metric), silent = TRUE)
  if (inherits(out, "try-error")) {
    warning(sprintf("Model failed for '%s': %s", metric, as.character(out)))
  } else {
    print(out$anova)
    results_list[[metric]] <- out
  }
}

## ---------------------------
## 4) Build summary table of p-values
## ---------------------------
create_summary_table <- function(results_list) {
  rows <- lapply(names(results_list), function(metric) {
    aov_tab <- results_list[[metric]]$anova
    # p column can vary by method; pick the first that exists
    p_col <- intersect(c("Pr(>F)", "p-value", "Pr(>Chisq)"), colnames(aov_tab))[1]
    if (is.na(p_col)) return(NULL)
    
    # fetch safely; use NA if row not present
    get_p <- function(term) {
      r <- aov_tab[rownames(aov_tab) == term, p_col, drop = TRUE]
      if (length(r) == 0) NA_real_ else as.numeric(r)
    }
    
    data.frame(
      Alpha_diversity = metric,
      Group           = get_p("Group"),
      Timepoint       = get_p("Timepoint"),
      Group_Timepoint = get_p("Group:Timepoint"),
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, rows)
  if (!is.null(out)) {
    out$Group           <- round(out$Group, 3)
    out$Timepoint       <- round(out$Timepoint, 3)
    out$Group_Timepoint <- round(out$Group_Timepoint, 3)
  }
  out
}

summary_table <- create_summary_table(results_list)

cat("\n==== Summary table (p-values) ====\n")
print(summary_table, row.names = FALSE)

## (Optional) save to CSV
write.csv(summary_table, "2508_oral_alpha_div_lmer_pvals.csv", row.names = FALSE)

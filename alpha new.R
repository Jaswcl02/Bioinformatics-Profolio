
#===overall alpha indices====
dataset_gut$alpha_diversity
alpha_dataset_gut.obj <- trans_alpha$new(dataset = dataset_gut, group = "GroupTime")
alpha_dataset_gut.obj$cal_diff(method = "wilcox")
alpha_dataset_gut.obj$data_stat

alpha_dataset_gut.obj$plot_alpha(measure = "InvSimpson",, add_sig = T,  add = "point",
                                         point_size = 3,
                                         point_alpha = 0.8)

alpha_dataset_gut.obj$plot_alpha(measure = "Shannon", add_sig = T,  add = "point",
                                 point_size = 3,
                                 point_alpha = 0.8)



#paired alpha diversity for dataset_placebo_gut, comparing difference in W0 and w12
alpha_dataset_placebo_gut.obj <- trans_alpha$new(dataset = dataset_placebo_gut, group = "Timepoint",by_ID = 'Subject_ID')
alpha_dataset_placebo_gut.obj$cal_diff(method = "wilcox")

alpha_dataset_placebo_gut.obj$plot_alpha(measure = "Shannon", add_sig = TRUE,  add = "point",
                                         point_size = 3,
                                         point_alpha = 0.8)
alpha_dataset_placebo_gut.obj$data_stat
write.csv(dataset_placebo_gut$alpha_diversity,file = "1208_dataset_placebo_gut_alpha_diversity_w0w12",row.names = TRUE)

alpha_dataset_placebo_gut.obj$res_diff
write.csv(as.data.frame(dataset_placebo_gut$res_diff), file = "1208_dataset_placebo_gut_res_diff_w0w12.csv", row.names = TRUE)

#paired alpha diversity for dataset_treatment_gut, comparing difference in W0 and w12
alpha_dataset_treatment_gut.obj <- trans_alpha$new(dataset = dataset_treatment_gut, group = "Timepoint",by_ID = 'Subject_ID')
alpha_dataset_treatment_gut.obj$cal_diff(method = "wilcox")
alpha_dataset_treatment_gut.obj$plot_alpha(measure = "Shannon", add_sig = TRUE,add = "point",
                                           point_size = 3,
                                           point_alpha = 0.8)
head(dataset_treatment_gut$alpha_diversity)
write.csv(dataset_treatment_gut$alpha_diversity,file = "1208_dataset_treatment_gut_alpha_diversity_w0w12",row.names = TRUE)
alpha_dataset_treatment_gut.obj$res_diff
write.csv(alpha_dataset_treatment_gut.obj$res_diff,file = "1208_dataset_treatment_gut_ares_diff_w0w12",row.names = TRUE)

#unpaired alpha diversity for dataset_w0_gut, comparing difference placebo and treatmenr in w0
alpha_dataset_w0_gut.obj <- trans_alpha$new(dataset = dataset_w0_gut, group = "Group")
alpha_dataset_w0_gut.obj$cal_diff(method = "wilcox")
alpha_dataset_w0_gut.obj$plot_alpha(measure = "Shannon", add_sig = TRUE,add = "point",
                                    point_size = 3,
                                    point_alpha = 0.8)
alpha_dataset_w0_gut.obj$res_diff
head(dataset_w0_gut$alpha_diversity)
write.csv(dataset_w0_gut$alpha_diversity,file = "dataset_w0_gut_alpha_diversity_PvsT",row.names = TRUE)

#unpaired alpha diversity for dataset_w12_gut, comparing difference placebo and treatmenr in w12
alpha_dataset_w12_gut.obj <- trans_alpha$new(dataset = dataset_w12_gut, group = "Group")
alpha_dataset_w12_gut.obj$cal_diff(method = "wilcox")
alpha_dataset_w12_gut.obj$plot_alpha(measure = "Shannon", add_sig = TRUE,add = "point",
                                     point_size = 3,
                                     point_alpha = 0.8)
head(dataset_w12_gut$alpha_diversity)
alpha_dataset_w12_gut.obj$res_diff
write.csv(dataset_w12_gut$alpha_diversity,file = "dataset_w12_gut_alpha_diversity_PvsT",row.names = TRUE)

#t.test

dataset_placebo_gut$alpha_diversity
ad <- dataset_placebo_gut$alpha_diversity
ad$Timepoint <- dataset_placebo_gut$sample_table[rownames(ad),"Timepoint"]

tt<-t.test(Shannon ~ Timepoint, data = ad)
tt

res_df <- data.frame(
  group1    = levels(factor(ad$Timepoint))[1],
  group2    = levels(factor(ad$Timepoint))[2],
  mean1     = tt$estimate[1],
  mean2     = tt$estimate[2],
  t_value   = unname(tt$statistic),
  df        = unname(tt$parameter),
  p_value   = tt$p.value,
  conf_low  = tt$conf.int[1],
  conf_high = tt$conf.int[2]
)

# Save as CSV
write.csv(res_df, "dataset_placebo_gut_alpha_diversity_w0w12_t.csv", row.names = FALSE)

####
dataset_treatment_gut$alpha_diversity
b <- dataset_treatment_gut$alpha_diversity
b$Timepoint <- dataset_treatment_gut$sample_table[rownames(b),"Timepoint"]

treat_t <-t.test(Shannon ~ Timepoint, data =b)
treat_t

dataset_treatment_gutW0W12_df <- data.frame(
  group1    = levels(factor(ad$Timepoint))[1],
  group2    = levels(factor(ad$Timepoint))[2],
  mean1     = treat_t$estimate[1],
  mean2     = treat_t$estimate[2],
  t_value   = unname(treat_t$statistic),
  df        = unname(treat_t$parameter),
  p_value   = treat_t$p.value,
  conf_low  = treat_t$conf.int[1],
  conf_high = treat_t$conf.int[2]
)

# Save as CSV
write.csv(dataset_treatment_gutW0W12_df, "dataset_treatment_gut_alpha_diversity_w0w12_t.csv", row.names = FALSE)

###lme
alpha_obj <- trans_alpha$new(dataset = dataset_gut)
# Check the structure of your alpha diversity data
str(alpha_obj$data_alpha)


alpha_dataset_gut.obj$cal_diff(
  method  = "lme",
  formula = "Group*Timepoint+ (1 | Subject_ID)",alpha = 0.05,return_model = TRUE )


# Look for missing values
sum(is.na(alpha_obj$data_alpha))

# Check sample sizes per group
table(alpha_obj$data_alpha$Group, alpha_obj$data_alpha$Timepoint)

# Check number of observations per subject
table(alpha_obj$data_alpha$Subject_ID)

# Convert categorical variables to factors
alpha_obj$data_alpha$Subject_ID <- as.factor(alpha_obj$data_alpha$Subject_ID)
alpha_obj$data_alpha$Group <- as.factor(alpha_obj$data_alpha$Group)
alpha_obj$data_alpha$Timepoint <- as.factor(alpha_obj$data_alpha$Timepoint)
alpha_obj$data_alpha$Measure <- as.factor(alpha_obj$data_alpha$Measure)

# See what measures you have
cat("Alpha diversity measures:\n")
print(table(alpha_obj$data_alpha$Measure))

# Work with just one measure first - let's use "Observed" (richness)
observed_data <- alpha_obj$data_alpha[alpha_obj$data_alpha$Measure == "Observed", ]

# Check the structure
cat("\nObserved richness data structure:\n")
cat("Dimensions:", dim(observed_data), "\n")
print(table(observed_data$Group, observed_data$Timepoint))

# Check subjects per group and timepoint
library(dplyr)
cat("\nSubjects per group:\n")
subject_summary <- observed_data %>% 
  group_by(Group) %>% 
  summarise(n_subjects = n_distinct(Subject_ID), n_obs = n()) %>%
  print()

cat("\nSubject-timepoint breakdown:\n")
timepoint_summary <- observed_data %>%
  group_by(Subject_ID, Group) %>%
  summarise(
    n_timepoints = n_distinct(Timepoint),
    timepoints = paste(sort(unique(Timepoint)), collapse = ", "),
    .groups = "drop"
  ) %>%
  print()

cat("\nSubjects by number of timepoints:\n")
print(table(timepoint_summary$n_timepoints))

#=====
# Identify subjects with both timepoints
complete_subjects <- timepoint_summary %>%
  filter(n_timepoints == 2) %>%
  pull(Subject_ID)

cat("Complete subjects (both timepoints):", length(complete_subjects), "\n")
cat("Incomplete subjects (one timepoint):", sum(timepoint_summary$n_timepoints == 1), "\n")

# Filter to complete subjects only
complete_observed <- observed_data %>%
  filter(Subject_ID %in% complete_subjects)

cat("\nComplete data structure:\n")
cat("Dimensions:", dim(complete_observed), "\n")
print(table(complete_observed$Group, complete_observed$Timepoint))

# Check that we now have balanced data
complete_check <- complete_observed %>%
  group_by(Subject_ID, Group) %>%
  summarise(n_timepoints = n_distinct(Timepoint), .groups = "drop")

cat("All subjects have 2 timepoints:", all(complete_check$n_timepoints == 2), "\n")

# Try the mixed-effects model
cat("Trying mixed-effects model with complete subjects:\n")
tryCatch({
  obs_alpha$cal_diff(method = "lme", formula = "Group*Timepoint + (1|Subject_ID)")
  cat("✓ SUCCESS! Mixed-effects model worked\n")
  print(obs_alpha$res_diff)
}, error = function(e) {
  cat("Still failed:", e$message, "\n")
  
  # Try without interaction as backup
  tryCatch({
    obs_alpha$cal_diff(method = "lme", formula = "Group + Timepoint + (1|Subject_ID)")
    cat("✓ Mixed-effects without interaction worked\n")
    print(obs_alpha$res_diff)
  }, error = function(e2) {
    cat("Mixed-effects still failing, using linear model:\n")
    obs_alpha$cal_diff(method = "lm", formula = "Group*Timepoint")
    print(obs_alpha$res_diff)
  })
})

str(dataset_gut$alpha_diversity)
rownames(dataset_gut$alpha_diversity)
rownames(dataset_gut$sample_table)

dataset_gut$sample_table$Alpha_Diversity <- dataset_gut$alpha_diversity$Shannon
str(dataset_gut$sample_table$Alpha_Diversity)



library(lme4)

model <- lmer(Alpha_Diversity ~ Timepoint * Group + (1 | Subject_ID), data = dataset_gut$sample_table)
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
  file = "lme_fixed_effects.csv",
  row.names = FALSE
)



library(lme4)
library(broom.mixed)
library(knitr)
library(dplyr)

# Get your complete subjects (from the previous analysis)
complete_subjects <- timepoint_summary %>%
  filter(n_timepoints == 2) %>%
  pull(Subject_ID)

# Filter to complete subjects only
complete_alpha_data <- alpha_obj$data_alpha %>%
  filter(Subject_ID %in% complete_subjects)

# Get all alpha diversity measures
all_measures <- unique(complete_alpha_data$Measure)
print(all_measures)

# Function to analyze each measure
analyze_alpha_measure <- function(measure_name, data) {
  cat("\n=== Analyzing", measure_name, "===\n")
  
  # Filter for specific measure
  measure_data <- data %>%
    filter(Measure == measure_name) %>%
    mutate(
      Subject_ID = as.factor(Subject_ID),
      Group = as.factor(Group),
      Timepoint = as.factor(Timepoint)
    )
  
  cat("Sample sizes:\n")
  print(table(measure_data$Group, measure_data$Timepoint))
  
  tryCatch({
    # Fit the model
    model <- lmer(Value ~ Timepoint * Group + (1 | Subject_ID), 
                  data = measure_data)
    
    # Extract results
    fixed_effects <- tidy(model, effects = "fixed")
    fixed_effects$Measure <- measure_name
    fixed_effects$Signif <- cut(
      fixed_effects$p.value,
      breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
      labels = c("***", "**", "*", ".", "")
    )
    
    cat("✓ Success for", measure_name, "\n")
    return(fixed_effects)
    
  }, error = function(e) {
    cat("✗ Failed for", measure_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Analyze all measures
all_results <- list()
for(measure in all_measures) {
  result <- analyze_alpha_measure(as.character(measure), complete_alpha_data)
  if(!is.null(result)) {
    all_results[[as.character(measure)]] <- result
  }
}

# Combine all results
combined_results <- do.call(rbind, all_results)

# Display results
print(kable(
  combined_results[, c("Measure", "term", "estimate", "std.error", "df", "statistic", "p.value", "Signif")],
  digits = 3,
  caption = "Fixed Effects from LME Models - All Alpha Diversity Measures"
))

# Save results
write.csv(
  combined_results,
  file = "all_alpha_diversity_lme_results.csv",
  row.names = FALSE
)



#======lme 
## ---------------------------
## 1) Merge sample + alpha data
## ---------------------------
sample_data <- dataset_gut$sample_table
alpha_data  <- dataset_gut$alpha_diversity

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
write.csv(summary_table, "2008_gut_alpha_div_lmer_pvals.csv"", row.names = FALSE)


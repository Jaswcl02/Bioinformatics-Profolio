dataset_treatment_gut$cal_betadiv(method = "bray")
dataset_treatment_gut.obj <- trans_beta$new(dataset = dataset_treatment_gut, measure = "bray", group = "Timepoint")
dataset_treatment_gut.obj$cal_ordination(method = "PCoA")
dataset_treatment_gut.obj$plot_ordination(
  plot_color = "Timepoint",
  plot_shape = "Timepoint",
  plot_type = c("point", "ellipse"))
dataset_treatment_gut.obj$res_ordination
write.csv(dataset_treatment_gut$beta_diversity, file = "dataset_treatment_gut_timepointcompare_beta.csv",row.names = TRUE)

#check dispersion before Manova
dataset_treatment_gut.obj$cal_betadisper()
dataset_treatment_gut.obj$res_betadisper 

#paired manova
dataset_treatment_gut.obj$cal_manova(
  group        = "Timepoint",
  permutations = 9999,
  strata       = dataset_treatment_gut.obj$sample_table$Subject_ID  # <-- pairing
)

head(dataset_treatment_gut.obj$res_manova)
write.csv(dataset_treatment_gut.obj$res_manova, file = "1508_beta_treatment_gut_w0w12res_manova.csv",row.names = TRUE)

dataset_treatment_gut.obj$cal_group_distance(
  within_group  = FALSE,                # between W0 and W12
  by_group      = "Subject_ID",         # pair within each subject
  ordered_group = c("W0","W12")         # keep the order
)
write.csv(dataset_treatment_gut.obj$cal_group_distance, file = "1508_beta_treatment_gut_w0w12_group_distance_paired.csv",row.names = TRUE)


####
dataset_placebo_gut$cal_betadiv(method = "bray")
dataset_placebo_gut.obj <- trans_beta$new(dataset = dataset_placebo_gut, measure = "bray", group = "Timepoint")
dataset_placebo_gut.obj$cal_ordination(method = "PCoA")
dataset_placebo_gut.obj$plot_ordination(plot_color = "Timepoint", plot_shape = "Timepoint", plot_type = c("point", "ellipse"))
write.csv(dataset_placebo_gut$beta_diversity, file = "dataset_placebo_gut_timepointcompare_beta.csv",row.names = TRUE)


#check dispersion before Manova
dataset_placebo_gut.obj$cal_betadisper()
dataset_placebo_gut.obj$res_betadisper 

#paired manova
dataset_placebo_gut.obj$cal_manova(
  group        = "Timepoint",
  permutations = 9999,
  strata       = dataset_placebo_gut.obj$sample_table$Subject_ID  # <-- pairing
)

head(dataset_placebo_gut.obj$res_manova)
write.csv(dataset_placebo_gut.obj$res_manova, file = "1508_beta_placebo_gut_w0w12res_manova.csv",row.names = TRUE)

dataset_placebo_gut.obj$cal_group_distance(
  within_group  = FALSE,                # between W0 and W12
  by_group      = "Subject_ID",         # pair within each subject
  ordered_group = c("W0","W12")         # keep the order
)
write.csv(dataset_placebo_gut.obj$cal_group_distance, file = "1508_beta_placebo_gut_w0w12_group_distance_paired.csv",row.names = TRUE)


####
dataset_w0_gut$cal_betadiv(method = "bray")
dataset_w0_gut.obj <- trans_beta$new(dataset = dataset_w0_gut, measure = "bray", group = "Group")
dataset_w0_gut.obj $cal_ordination(method = "PCoA")
dataset_w0_gut.obj $plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type =  c("point", "ellipse"))
write.csv(dataset_w0_gut$beta_diversity, file = "dataset_w0_gut_groupcompare_beta.csv",row.names = TRUE)

dataset_w0_gut.obj$cal_manova(group = "Group", permutations = 9999)  
head(dataset_w0_gut.obj $res_manova)
write.csv(dataset_w0_gut.obj $res_manova, file = "1508_dataset_w0_gut_group_res_manova.csv",row.names = TRUE)
dataset_w0_gut.obj $res_manova
dataset_w0_gut.obj $cal_group_distance(within_group = TRUE)
dataset_w0_gut.obj$res_group_distance

write.csv(dataset_w0_gut.obj$res_group_distance, file = "1508_dataset_w0_gut_group_res_group_distance_.csv",row.names = TRUE)


dataset_w0_gut.obj $plot_group_distance()

####
dataset_w12_gut$cal_betadiv(method = "bray")
dataset_w12_gut.obj <- trans_beta$new(dataset = dataset_w12_gut, measure = "bray", group = "Group")
dataset_w12_gut.obj $cal_ordination(method = "PCoA")
dataset_w12_gut.obj $plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type =  c("point", "ellipse"))
write.csv(dataset_w12_gut.obj$beta_diversity, file = "dataset_w12_gut_groupcompare_beta.csv",row.names = TRUE)

dataset_w12_gut.obj$cal_manova(group = "Group", permutations = 9999)  
head(dataset_w12_gut.obj $res_manova)
write.csv(dataset_w12_gut.obj $res_manova, file = "1508_dataset_w12_gut_group_res_manova.csv",row.names = TRUE)
dataset_w12_gut.obj $res_manova
dataset_w12_gut.obj $cal_group_distance(within_group = TRUE)
dataset_w12_gut.obj$res_group_distance

write.csv(dataset_w12_gut.obj$res_group_distance, file = "1508_dataset_w12_gut_group_res_group_distance_.csv",row.names = TRUE)
dataset_w12_gut.obj $plot_group_distance()


###overall (bary)
library(microeco)
library(vegan)   # adonis2, vegdist
library(ape)     # pcoa
library(ggplot2)

# Compute beta once in microeco (optional; we’ll still use vegan for PERMANOVA)


# Group labels like your slide
st <- dataset_gut$sample_table
st$GroupTime <- paste0(st$Group, "_", st$Timepoint)
st$GroupTime <- factor(st$GroupTime,
                       levels = c("Placebo_W0","Placebo_W12","MHLJDD_W0","MHLJDD_W12"))
vars <- c("GroupTime","BMI","Age","Sex","Group","Timepoint","Subject_ID")
st <- st[complete.cases(st[, vars, drop = FALSE]), , drop = FALSE]

st$Delivery_mode <- droplevels(factor(trimws(as.character(st$Delivery_mode))))
st$Sex         <- droplevels(factor(trimws(as.character(st$Sex))))
st$Group          <- droplevels(factor(trimws(as.character(st$Group))))
st$Timepoint      <- droplevels(factor(trimws(as.character(st$Timepoint))))
st$Subject_ID     <- droplevels(factor(trimws(as.character(st$Subject_ID))))
st$BMI            <- as.numeric(st$BMI)
st$Age            <- as.numeric(st$Age)


dataset_gut$sample_table <- st
dataset_gut$sample_table
# Make sure the metadata and OTU matrix align (rows = samples)
otu <- t(dataset_gut$otu_table)
st  <- st[rownames(otu), , drop = FALSE]

###


#beta diversity and plot
dataset_gut$cal_betadiv(method = "bray",Group= "GroupTime")
t2 <-trans_beta$new(dataset =dataset_gut, measure = "bray", group = "GroupTime")
dataset_gut$sample_table
t2$cal_ordination(method = "PCoA")
t2$plot_ordination(plot_color = "GroupTime", plot_type = c("point", "ellipse"))

write.csv(dataset_gut$beta_diversity, file = "dataset__gut_GroupTime_beta.csv",row.names = TRUE)

####
t2$cal_manova(
  manova_set = "Delivery_mode + BMI + Age + Sex + Batch + Group+ GroupTime + Timepoint + Group:Timepoint",strata = t1$sample_table$Subject_ID 
)
t2$res_manova
t2$cal_manova(manova_set = "Sex+Age+Timepoint+Group+Group:Timepoint",strata = t1$sample_table$Subject_ID )
t2$res_manova

t2$cal_group_distance(
  within_group = TRUE,
  by_group = "Subject_ID"  # This pairs by subject
)
t2$cal_group_distance_diff(method = "wilcox")
t2$plot_group_distance(add = "mean",add_sig=F)

#==btw group distance

t2$cal_group_distance(within_group = FALSE)
t2$cal_group_distance_diff(method = "wilcox")
t2$plot_group_distance(plot_type="ggviolin", plot_group_order = c("MHLJDD_W0 vs MHLJDD_W12","Placebo_W0 vs Placebo_W12","MHLJDD_W0 vs Placebo_W0","MHLJDD_W12 vs Placebo_W12"),add_sig = F, 
                       add=c("mean","mean_se"),xtext_angle=60)
t2$plot_group_distance(add = "mean")

###overall (cal_group_distance)
## ================== SETUP ==================
library(dplyr)
library(tidyr)
library(vegan)      # vegdist
library(ggplot2)
library(ggpubr)     # compare_means, stat_pvalue_manual
library(patchwork)

## dataset_gut: microeco::microtable (already in your env)
## Uses: dataset_gut$otu_table (features x samples), dataset_gut$sample_table (samples x metadata)

## ================== 1) META + ORDER ==================
meta <- dataset_gut$sample_table

# Build GroupTime like "Placebo_W0" etc.
if (!"GroupTime" %in% names(meta)) {
  meta$GroupTime <- with(meta, paste0(Group, "_", Timepoint))
}
# Order exactly as you want them to show (edit if needed)
lev <- c("Placebo_W0","Placebo_W12","MHLJDD_W0","MHLJDD_W12")
meta$GroupTime <- factor(meta$GroupTime, levels = lev)
meta <- droplevels(meta)

# Align OTU table columns to meta rownames (samples)
otu <- dataset_gut$otu_table
otu <- otu[, rownames(meta), drop = FALSE]

## ================== 2) HELPERS TO GET DISTANCES ==================
# (A) Within-group: for each GroupTime, all pairwise distances among samples in that group
get_within_bc <- function(group_label){
  s <- rownames(meta)[meta$GroupTime == group_label]
  if (length(s) < 2) return(NULL)
  # vegdist expects rows = samples; our otu has rows = features, cols = samples -> transpose
  m <- as.matrix(vegdist(t(otu[, s, drop = FALSE]), method = "bray"))
  tibble(Group = group_label, Distance = m[upper.tri(m, diag = FALSE)])
}

# (B) Between-group: all cross-sample distances between two GroupTime levels
get_between_bc <- function(g1, g2, label = NULL){
  s1 <- rownames(meta)[meta$GroupTime == g1]
  s2 <- rownames(meta)[meta$GroupTime == g2]
  if (length(s1) == 0 || length(s2) == 0) return(NULL)
  m  <- as.matrix(vegdist(t(otu[, c(s1, s2), drop = FALSE]), method = "bray"))
  m12 <- m[s1, s2, drop = FALSE]
  tibble(Pair = ifelse(is.null(label), paste(g1, "vs", g2), label),
         Distance = as.vector(m12))
}

pair_label <- function(a,b) paste(a, "vs", b)

## ================== 3) BUILD dfA (WITHIN) & dfB (BETWEEN) ==================
# dfA: within-group distributions for each GroupTime (Panel A)
dfA <- do.call(bind_rows, lapply(lev, get_within_bc))
if (is.null(dfA) || nrow(dfA) == 0) stop("dfA is empty: not enough samples per GroupTime for within-group distances.")
dfA$Group <- factor(dfA$Group, levels = lev)

# dfB: your 4 between-group distributions (Panel B)
pairs_to_use <- c(
  pair_label(lev[1], lev[2]),  # Placebo_W0 vs Placebo_W12
  pair_label(lev[3], lev[4]),  # MHLJDD_W0 vs MHLJDD_W12
  pair_label(lev[3], lev[1]),  # MHLJDD_W0 vs Placebo_W0
  pair_label(lev[4], lev[2])   # MHLJDD_W12 vs Placebo_W12
)

dfB <- bind_rows(
  get_between_bc(lev[1], lev[2], pairs_to_use[1]),
  get_between_bc(lev[3], lev[4], pairs_to_use[2]),
  get_between_bc(lev[3], lev[1], pairs_to_use[3]),
  get_between_bc(lev[4], lev[2], pairs_to_use[4])
)
if (is.null(dfB) || nrow(dfB) == 0) stop("dfB is empty: check your groups/labels for between-group distances.")
dfB$Pair <- factor(dfB$Pair, levels = pairs_to_use)

## ================== 4) COLORS ==================
cols_A <- c(                      # Panel A (GroupTime)
  "Placebo_W0"  = "#1b9e77",      # green
  "Placebo_W12" = "#d95f02",      # orange
  "MHLJDD_W0"   = "#7570b3",      # purple
  "MHLJDD_W12"  = "#e7298a"       # pink
)
# FIXED: use setNames so names come from `pairs_to_use`
cols_B <- setNames(
  c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
  pairs_to_use
)

## ================== 5) SIGNIFICANCE (pre-compute & force show) ==================
# A: we only want (lev[1] vs lev[2]) and (lev[3] vs lev[4])
allA <- ggpubr::compare_means(Distance ~ Group, data = dfA, method = "wilcox.test")

testsA <- allA %>%
  dplyr::filter(
    (group1 %in% c(lev[1], lev[2]) & group2 %in% c(lev[1], lev[2])) |
      (group1 %in% c(lev[3], lev[4]) & group2 %in% c(lev[3], lev[4]))
  ) %>%
  # order brackets in a stable L->R fashion using factor order in `lev`
  dplyr::arrange(
    pmin(match(group1, lev), match(group2, lev)),
    pmax(match(group1, lev), match(group2, lev))
  )

# Put brackets inside the plot area; length must match nrow(testsA)
ymaxA <- max(dfA$Distance, na.rm = TRUE)
testsA$y.position <- ymaxA * (1.08 + 0.08 * (seq_len(nrow(testsA)) - 1))

# B: we only want the bracket between the two cross-arm pairs at T0 and T12
allB <- ggpubr::compare_means(Distance ~ Pair, data = dfB, method = "wilcox.test")

testsB <- allB %>%
  dplyr::filter(
    group1 %in% c(pairs_to_use[3], pairs_to_use[4]) &
      group2 %in% c(pairs_to_use[3], pairs_to_use[4])
  ) %>%
  dplyr::slice(1)  # keep just that single comparison

ymaxB <- max(dfB$Distance, na.rm = TRUE)
testsB$y.position <- ymaxB * 1.10

## ================== 6) PLOTS ==================
pA <- ggplot(dfA, aes(Group, Distance, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  scale_fill_manual(values = cols_A, guide = "none") +
  coord_cartesian(ylim = c(0, ymaxA*1.25), clip = "off") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.margin = unit(c(8,16,8,8), "pt")) +
  labs(x = NULL, y = "Bray–Curtis distance") +
  ggpubr::stat_pvalue_manual(
    testsA,
    label = "p.signif",        # "p.format" for numeric values
    xmin  = "group1",
    xmax  = "group2",
    y.position = "y.position",
    tip.length  = 0.01, bracket.size = 0.6
  )

pB <- ggplot(dfB, aes(Pair, Distance, fill = Pair)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  scale_fill_manual(values = cols_B, guide = "none") +
  coord_cartesian(ylim = c(0, ymaxB*1.25), clip = "off") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        plot.margin = unit(c(8,16,8,8), "pt")) +
  labs(x = NULL, y = "Bray–Curtis distance") +
  ggpubr::stat_pvalue_manual(
    testsB,
    label = "p.signif",
    xmin  = "group1",
    xmax  = "group2",
    y.position = "y.position",
    tip.length  = 0.01, bracket.size = 0.6
  )

pB
final_plot <- (pA | pB) + plot_annotation(
  title = "Within- and between-group differences in Bray–Curtis distance"
)
final_plot

## ================== 7) (OPTIONAL) GLOBAL PERMANOVA WITH PAIRING ==================
# Accounts for repeated measures by subject as a block/strata:
otuM  <- t(dataset_gut$otu_table[, rownames(dataset_gut$sample_table), drop = FALSE]) # samples x features
distM <- vegdist(otuM, method = "bray")
adon  <- adonis2(
  distM ~ Group * Timepoint,
  data = dataset_gut$sample_table,
  permutations = 9999,
  strata = dataset_gut$sample_table$Subject_ID
)
print(adon)

#=====

#import cleaned files and install library
install.packages("BiocManager")
install.packages("file2meco", repos = BiocManager::repositories())
install.packages("MicrobiomeStat", repos = BiocManager::repositories())
install.packages("WGCNA", repos = BiocManager::repositories())
BiocManager::install("ggtree")
BiocManager::install("metagenomeSeq")
BiocManager::install("ALDEx2")
BiocManager::install("ANCOMBC")


# allow more waiting time to download each package
options(timeout = 1000)
# If a package is not installed, it will be installed from CRAN
# First select the packages of interest
tmp <- c("microeco", "mecoturn", "MASS", "GUniFrac", "ggpubr", "randomForest", "ggdendro", "ggrepel", "agricolae", "igraph", "picante", "pheatmap", "rgexf", 
         "ggalluvial", "ggh4x", "rcompanion", "FSA", "gridExtra", "aplot", "NST", "GGally", "ggraph", "networkD3", "poweRlaw", "ggtern", "SRS", "performance")
# Now check or install
for(x in tmp){
  if(!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
  }
}

library(microeco)
otu <- read.table("~/Documents/biomedtech/OneDrive_1_15-08-2025/otu_cleaned.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE,quote = "",comment.char="",colClasses ="character")
tax <- read.table("~/Documents/biomedtech/OneDrive_1_15-08-2025/taxonomy_cleaned.tsv", sep="\t", header=TRUE,row.names = 1,check.names = FALSE,quote = "",comment.char="")


otu <- data.frame(lapply(otu,as.numeric),row.names=rownames(otu),check.names=FALSE)

meta <- read.table("~/Documents/biomedtech/OneDrive_1_15-08-2025/sample_metadata.tsv", sep="\t", header=TRUE, row.names=1,check.names = FALSE,quote="",comment.char = "")
meta <- meta[!(rownames(meta) %in%"categorical"), ]


colnames(tax) <- trimws(colnames(tax))#trim and clean tax and otu to make dataset
rownames(otu) <- trimws(rownames(otu))

all(colnames(otu) %in% rownames(tax))#confirm matching on col and row names tax and otu to make dataset
all(rownames(tax) %in% colnames(otu))

colnames(otu) <- trimws(colnames(otu))#trim and clean meta and otu to make dataset
rownames(meta) <- trimws(rownames(meta))

all(colnames(otu) %in% rownames(meta))#confirm matching on col and row names tax and otu to make dataset
all(rownames(meta) %in% colnames(otu))


dataset <- microtable$new(otu_table = otu, tax_table = tax, sample_table = meta) # make dataset
dataset$tidy_dataset()

# STEP 1: Remove NAs from OTU table
# Replace NAs with 0 (common in microbiome data)
otu[is.na(otu)] <- 0


#separate taxon of tax

tax_table_with_id <- cbind(FeatureID = rownames(tax),tax)
head(tax_table_with_id)

# Step 2: Function to handle splitting and missing assignments
split_tax <- function(taxon) {
  if (taxon == "Unassigned" | is.na(taxon) | taxon == "") {
    # Return all NAs for unassigned
    return(rep(NA, 7))
  }
  # Remove possible leading underscores, then split by ";"
  levels <- unlist(strsplit(gsub("^_+", "", taxon), ";"))
  # Pad with NA if not all 7 levels present
  if (length(levels) < 7) levels <- c(levels, rep(NA, 7 - length(levels)))
  return(levels)
}

# Step 3: Apply splitting to every row
tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_split <- t(sapply(tax_table_with_id$Taxon, split_tax))
colnames(tax_split) <- tax_levels

# Step 4: Combine back with the original table (remove original 'Taxon' column if you wish)
tax_table_new <- cbind(tax_table_with_id[ , !(names(tax_table_with_id) %in% "Taxon")], tax_split)
#replace tax with new tax
dataset$tax_table <- tax_table_new

#check alignment 
head(rownames(dataset$otu_table)) 
head(rownames(tax_table_new))
tax_table_new <- as.data.frame(tax_table_new,stringAsFactor = FALSE)
rownames(tax_table_new) <- tax_table_new$FeatureID
tax_table_new$FeatureID <- NULL

dataset$tax_table<- tax_table_new
identical(rownames(dataset$otu_table),rownames(dataset$tax_table))

tidy_taxonomy(
  dataset$otu_table,
  column = "all",
  pattern = c(".*unassigned.*", ".*uncultur.*", ".*unknown.*", ".*unidentif.*",
              ".*unclassified.*", ".*No blast hit.*", ".*Incertae.sedis.*"),
  replacement = "",
  ignore.case = TRUE,
  na_fill = ""
)

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

#spilting dataset into stool subset 
meta_stool <- subset(dataset$sample_table, 
                     grepl("Stool", SampleType, ignore.case = TRUE))
write.csv(meta_stool,file = "meta_stool.csv",row.names = T)

## Ensure OTU table matches new sample list
otu_stool <- dataset$otu_table[, rownames(meta_stool), drop = FALSE]

## Build gut/stool dataset
dataset_gut <- microtable$new(otu_table = otu_stool,
                              tax_table = dataset$tax_table,
                              sample_table = meta_stool)
dataset_gut$tidy_dataset()

#spilting dataset into saliva subset 
meta_saliva <- subset(dataset$sample_table, 
                     grepl("Saliva", SampleType, ignore.case = TRUE))

## Ensure OTU table matches new sample list
otu_saliva <- dataset$otu_table[, rownames(meta_saliva), drop = FALSE]

## 3) Build saliva dataset
dataset_saliva <- microtable$new(otu_table = otu_saliva,
                              tax_table = dataset$tax_table,
                              sample_table = meta_saliva)
dataset_saliva$tidy_dataset()


###spilting stool dataset into treatment subset 
meta_treatment_gut <- subset(dataset_gut$sample_table, 
                      grepl("MHLJDD", Group, ignore.case = TRUE))

## Ensure OTU table matches new sample list
otu_treatment_gut <- dataset_gut$otu_table[, rownames(meta_treatment_gut), drop = FALSE]

## 3) Build gut/stool dataset
dataset_treatment_gut <- microtable$new(otu_table = otu_treatment_gut,
                                 tax_table = dataset$tax_table,
                                 sample_table = meta_treatment_gut)
dataset_treatment_gut$tidy_dataset()

#spilting stool dataset Placebo-gut_ds
meta_placebo_gut <- subset(dataset_gut$sample_table, 
                             grepl("Placebo", Group, ignore.case = TRUE))

## Ensure OTU table matches new sample list
otu_placebo_gut <- dataset_gut$otu_table[, rownames(meta_placebo_gut), drop = FALSE]

## 3) Build gut/stool dataset
dataset_placebo_gut <- microtable$new(otu_table = otu_placebo_gut,
                                        tax_table = dataset$tax_table,
                                        sample_table = meta_placebo_gut)
dataset_placebo_gut$tidy_dataset()

###spilting stool dataset w0_ds 
meta_w0_gut <- subset(dataset_gut$sample_table, 
                           grepl("W0", Timepoint, ignore.case = TRUE))

## Ensure OTU table matches new sample list
otu_w0_gut <- dataset_gut$otu_table[, rownames(meta_w0_gut), drop = FALSE]

## 3) Build gut/stool dataset
dataset_w0_gut <- microtable$new(otu_table = otu_w0_gut,
                                      tax_table = dataset$tax_table,
                                      sample_table = meta_w0_gut)
dataset_w0_gut$tidy_dataset()

###spilting stool dataset w12 dataset
meta_w12_gut <- subset(dataset_gut$sample_table, 
                           grepl("W12", Timepoint, ignore.case = TRUE))

## Ensure OTU table matches new sample list
otu_w12_gut <- dataset_gut$otu_table[, rownames(meta_w12_gut), drop = FALSE]

## 3) Build gut/stool dataset
dataset_w12_gut <- microtable$new(otu_table = otu_w12_gut,
                                      tax_table = dataset$tax_table,
                                      sample_table = meta_w12_gut)
dataset_w12_gut$tidy_dataset()

a


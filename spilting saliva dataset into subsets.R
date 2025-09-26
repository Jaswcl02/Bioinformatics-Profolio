###spilting saliva dataset into treatment subset 
meta_treatment_oral <- subset(dataset_saliva$sample_table, 
                             grepl("MHLJDD", Group, ignore.case = TRUE))

## Ensure OTU table matches new sample list
otu_treatment_oral <- dataset_saliva$otu_table[, rownames(meta_treatment_oral), drop = FALSE]

## 3) Build gut/stool dataset
dataset_treatment_oral <- microtable$new(otu_table = otu_treatment_oral,
                                        tax_table = dataset$tax_table,
                                        sample_table = meta_treatment_oral)
dataset_treatment_oral$tidy_dataset()

#spilting saliva dataset Placebo-oral_ds
meta_placebo_oral <- subset(dataset_saliva$sample_table, 
                           grepl("Placebo", Group, ignore.case = TRUE))

## Ensure OTU table matches new sample list
otu_placebo_oral <- dataset_saliva$otu_table[, rownames(meta_placebo_oral), drop = FALSE]

## 3) Build gut/stool dataset
dataset_placebo_oral <- microtable$new(otu_table = otu_placebo_oral,
                                      tax_table = dataset$tax_table,
                                      sample_table = meta_placebo_oral)
dataset_placebo_oral$tidy_dataset()

###spilting saliva dataset w0_ds 
meta_w0_oral <- subset(dataset_saliva$sample_table, 
                      grepl("W0", Timepoint, ignore.case = TRUE))

## Ensure OTU table matches new sample list
otu_w0_oral <- dataset_saliva$otu_table[, rownames(meta_w0_oral), drop = FALSE]

## 3) Build gut/stool dataset
dataset_w0_oral <- microtable$new(otu_table = otu_w0_oral,
                                 tax_table = dataset$tax_table,
                                 sample_table = meta_w0_oral)
dataset_w0_oral$tidy_dataset()

###spilting saliva dataset w12 dataset
meta_w12_oral <- subset(dataset_saliva$sample_table, 
                       grepl("W12", Timepoint, ignore.case = TRUE))

## Ensure OTU table matches new sample list
otu_w12_oral <- dataset_saliva$otu_table[, rownames(meta_w12_oral), drop = FALSE]

## 3) Build gut/stool dataset
dataset_w12_oral <- microtable$new(otu_table = otu_w12_oral,
                                  tax_table = dataset$tax_table,
                                  sample_table = meta_w12_oral)
dataset_w12_oral$tidy_dataset()
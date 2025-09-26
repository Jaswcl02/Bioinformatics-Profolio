Bioinformatics Portfolio – Microeco Microbiome Analysis
Overview

Hi yeah this is my personal bioinformatics portfolio demonstrating how to process, analyze and present microbiome sequencing data using R and supporting tools. The project focuses on building clean datasets suitable for downstream analysis, performing exploratory and inferential statistics, and compiling results into professional presentations.

Data processing

Raw OTU tables, taxonomic assignments and sample metadata were cleaned and imported into R.

The microeco package was used to combine OTU, taxonomic and sample tables into a single microtable object.

Data cleaning steps included trimming whitespace from identifiers, ensuring matching row/column names across tables, replacing missing OTU values with zeros, splitting hierarchical taxonomy strings into separate columns (Kingdom → Species), and filtering out unassigned/uncultured taxa.

Sample metadata were tidied by splitting the Sample_type column into Timepoint and SampleType and constructing derived variables (GroupTime, GroupSample).

Subsets of the full dataset were generated for stool and saliva samples, for specific treatment groups (e.g., MHLJDD vs. placebo), and for baseline (W0) vs. post‑treatment (W12) time points.

Statistical analyses

A variety of statistical techniques were applied to explore and compare microbial communities:

Alpha diversity: Richness and diversity indices were computed to summarise within‑sample diversity. Differences between groups (treatment vs. control, stool vs. saliva, time points) were assessed with parametric (t‑tests, ANOVA) or non‑parametric (Wilcoxon rank‑sum, Kruskal–Wallis) tests depending on distributional assumptions.

Beta diversity: Dissimilarity matrices (e.g., Bray‑Curtis) were computed and visualised with ordination methods such as PCA and PCoA. Permutational multivariate analysis of variance (PERMANOVA) tested for overall community composition differences.

Differential abundance: Compositional differential abundance methods (ALDEx2, ANCOM‑BC, metagenomeSeq) were used to identify taxa with significant changes between groups while controlling for false discovery rates.

Correlation and network analysis: Pearson and Spearman correlations quantified relationships between taxa and metadata variables. Weighted gene co‑expression network analysis (WGCNA) explored co‑abundance modules and their association with phenotypes.

Random forest and biomarker discovery: Machine‑learning methods (random forest classifiers) were trained to identify microbial signatures that distinguish treatment groups or sample types.

Visualisation and presentation

Results were visualised using ggplot2 and companion packages (ggpubr, ggalluvial, pheatmap, etc.), producing box plots, violin plots, stacked bar charts, ordination biplots, dendrograms and heatmaps.

Figures and key findings were assembled into a slide deck via the PptxGenJS library. A Python script (pptx_to_img.py) converts PPTX files to images, and create_montage.py builds montages for summarising multiple slides.

A detailed write‑up of the workflow is provided in project_summary.txt, and the R script microeco_dataset_subsetting.R shows how the dataset was constructed and subsetted.

Getting started

Clone the repository and install the required R packages listed in microeco_dataset_subsetting.R.

Follow the steps in the script to build the microtable object and generate subsets.

Run the analysis scripts (to be added) to reproduce the statistics and visualisations.

View the project_summary.txt for a narrative overview of the project and the analyses performed.

This README summarises the major steps and analyses performed in the project and provides a roadmap for exploring and extending the portfolio.

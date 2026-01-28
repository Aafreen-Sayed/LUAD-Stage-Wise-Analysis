# -------------------------------------------------------------------------
# Script: 02_StageWise_DGE_Pipeline.R
# Purpose: Stage-wise evolutionary differential expression analysis (Adjacent Comparisons)
# Method: Linear Models for Microarray Data (limma)
# -------------------------------------------------------------------------

library(limma)
library(openxlsx)
library(dplyr)

# Load data from Step 1
load("cleaned_data.RData")

# FUNCTION: Perform Pairwise Stage Comparison
get_deg_contrast <- function(stage1, stage2, exprs_data, pheno) {
  
  message("Running Contrast: ", stage1, " vs ", stage2)
  
  # Subset data for these two stages
  samples1 <- rownames(pheno)[pheno$clean_stage == stage1]
  samples2 <- rownames(pheno)[pheno$clean_stage == stage2]
  keep_samples <- c(samples1, samples2)
  exp_subset <- exprs_data[, keep_samples]
  
  # Design Matrix
  group <- factor(ifelse(colnames(exp_subset) %in% samples1, stage1, stage2))
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  # Fit Linear Model
  fit <- lmFit(exp_subset, design)
  
  # -------------------------------------------------------
  # KEY EVOLUTIONARY LOGIC: Defining the Dynamic Contrast
  # Comparing Stage N vs Stage N-1 to capture transition events
  # -------------------------------------------------------
  contrast_formula <- paste(stage2, "-", stage1) # e.g., "Stage 2 - Stage 1"
  contrast <- makeContrasts(contrasts = contrast_formula, levels = design)
  
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # Extract Results
  results <- topTable(fit2, number = Inf, sort.by = "P")
  results$Stage_Transition <- paste(stage1, "to", stage2)
  results$Gene.symbol <- rownames(results) # Assuming rownames are gene symbols
  
  return(results)
}

# Run the Sequential "Evolutionary" Contrasts
deg_1_vs_2 <- get_deg_contrast("Stage 1", "Stage 2", exprs_data, valid_pheno)
deg_2_vs_3 <- get_deg_contrast("Stage 2", "Stage 3", exprs_data, valid_pheno)
deg_3_vs_4 <- get_deg_contrast("Stage 3", "Stage 4", exprs_data, valid_pheno)

# Save Results
output_file <- "LUAD_Evolutionary_DEGs.xlsx"
wb <- createWorkbook()

addWorksheet(wb, "Stage1_to_Stage2")
writeData(wb, "Stage1_to_Stage2", deg_1_vs_2)

addWorksheet(wb, "Stage2_to_Stage3")
writeData(wb, "Stage2_to_Stage3", deg_2_vs_3)

addWorksheet(wb, "Stage3_to_Stage4")
writeData(wb, "Stage3_to_Stage4", deg_3_vs_4)

saveWorkbook(wb, output_file, overwrite = TRUE)
message("✅ Stage-wise analysis complete. Results saved to: ", output_file)
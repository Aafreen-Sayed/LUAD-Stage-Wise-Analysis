# -------------------------------------------------------------------------
# Script: 01_Data_Preprocessing.R
# Purpose: Preprocessing of raw GEO microarray data and clinical stage standardization
# Author: Aafreen Sayed
# -------------------------------------------------------------------------

library(GEOquery)
library(limma)
library(dplyr)
library(stringr)

# 1. Load Data (Example using one of your GSE IDs)
# Note: In practice, this runs across the multiple datasets listed in the methodology
gse_id <- "GSE33745" 
message("Processing dataset: ", gse_id)

# Download and parse dataset
gset <- getGEO(gse_id, GSEMatrix = TRUE)[[1]]
exprs_data <- exprs(gset)
pheno <- pData(gset)

# 2. Stage Cleaning Logic (The "Evolutionary" Framework Setup)
# Identify the correct column for staging info dynamically
stage_col <- grep("stage|tumou?r_stage|clinical_stage|pathologic_stage|ajcc", 
                  colnames(pheno), ignore.case = TRUE, value = TRUE)[1]

if (is.na(stage_col)) stop("❌ No recognizable stage column found.")

# Clean the staging column (Standardizing Roman Numerals I, II, III, IV)
pheno$clean_stage <- toupper(pheno[[stage_col]])
pheno$clean_stage <- gsub("STAGE\\s*", "", pheno$clean_stage)
pheno$clean_stage <- gsub("STG\\s*", "", pheno$clean_stage)

# Convert Roman numerals to numeric (Crucial for sequential ordering)
pheno$clean_stage <- gsub("IV", "4", pheno$clean_stage)
pheno$clean_stage <- gsub("III", "3", pheno$clean_stage)
pheno$clean_stage <- gsub("II", "2", pheno$clean_stage)
pheno$clean_stage <- gsub("I", "1", pheno$clean_stage)

# Remove substage letters (merging IA/IB into Stage 1)
pheno$clean_stage <- gsub("[^1-4]", "", pheno$clean_stage)

# Format final stage label
pheno$clean_stage <- paste0("Stage ", pheno$clean_stage)

# Filter only valid stage samples (Stage 1 to 4)
valid_pheno <- pheno[pheno$clean_stage %in% paste("Stage", 1:4), ]
exprs_data <- exprs_data[, rownames(valid_pheno)]

message("✅ Preprocessing Complete. Samples retained: ", nrow(valid_pheno))

# Save processed objects for the next step
save(exprs_data, valid_pheno, file = "cleaned_data.RData")
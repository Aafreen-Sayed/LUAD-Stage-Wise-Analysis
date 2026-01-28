# -------------------------------------------------------------------------
# Script: 03_Functional_Enrichment.R
# Purpose: Functional enrichment analysis (GO & KEGG) of stage-specific drivers
# Tools: clusterProfiler, org.Hs.eg.db
# -------------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(openxlsx)
library(ggplot2)

# Load the DEGs generated in Step 2
# (In a real run, read the Excel file produced in Step 2)
# degs <- read.xlsx("LUAD_Evolutionary_DEGs.xlsx", sheet = "Stage1_to_Stage2")

# Example: Run enrichment on a list of significant genes
# Filter for significant upregulation (logFC > 1, p < 0.05)
sig_genes_symbol <- degs$Gene.symbol[degs$logFC > 1 & degs$adj.P.Val < 0.05]

# Convert Gene Symbols to ENTREZ IDs (Required for KEGG/GO)
gene_entrez <- bitr(sig_genes_symbol, 
                    fromType = "SYMBOL",
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# 1. GO Enrichment (Biological Process)
message("Running GO Enrichment...")
ego <- enrichGO(gene = gene_entrez$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP", # Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# 2. KEGG Pathway Enrichment
message("Running KEGG Enrichment...")
ekegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

# 3. Visualization (The "Money Plots")
# Dotplot for top pathways
p1 <- dotplot(ego, showCategory = 15, title = "GO Enrichment: Stage Transition Drivers")

# Network plot (Gene-Pathway linkages)
p2 <- cnetplot(ego, showCategory = 5, circular = TRUE, colorEdge = TRUE)

# Save Plots
pdf("Enrichment_Visualization.pdf", width = 12, height = 8)
print(p1)
print(p2)
dev.off()

message("✅ Enrichment analysis complete. PDF saved.")
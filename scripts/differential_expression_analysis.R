# ============================================================
# Identifying Druggable Targets in Alzheimer's Disease
# Differential expression, clustering, and functional enrichment
# Dataset: GSE5281 (GEO)
# ============================================================

# --- Setup ---
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "affy", "clusterProfiler",
                        "org.Hs.eg.db", "GO.db"), update = FALSE)
install.packages(c("ggplot2", "pheatmap"))

library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# --- Load dataset ---
gse <- getGEO("GSE5281", GSEMatrix = TRUE)
expr_data <- exprs(gse[[1]])
metadata <- pData(gse[[1]])

# --- Assign groups from sample title (control vs AD) ---
group <- ifelse(grepl("control", metadata$title), "control", "AD")
group <- factor(group)

# --- Preprocessing: log2-transform (raw data was on a linear scale) ---
expr_data <- log2(expr_data + 1)

# --- Differential expression (limma) ---
design <- model.matrix(~group)
fit <- lmFit(expr_data, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")

# --- Significant genes ---
sig_genes <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
nrow(sig_genes)  # 2972

# --- Map probe IDs to gene symbols via platform annotation ---
gpl <- getGEO("GPL570", AnnotGPL = TRUE)
platform_data <- Table(gpl)
id_to_symbol <- setNames(platform_data$`Gene symbol`, platform_data$ID)
sig_genes$gene_symbol <- id_to_symbol[rownames(sig_genes)]
write.csv(sig_genes, "sig_genes_with_symbols.csv")

# --- Figure 1: Volcano plot ---
results$significant <- ifelse(
  results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
  "Significant", "Not significant"
)
ggplot(results, aes(x = logFC, y = -log10(P.Value), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey")) +
  labs(title = "Volcano plot: AD vs Control (GSE5281)",
       x = "log2 fold change", y = "-log10 p-value", color = "") +
  theme_minimal()

# --- Figure 2: Heatmap of top 30 DEGs ---
top30 <- rownames(sig_genes)[order(sig_genes$adj.P.Val)][1:30]
heatmap_data <- expr_data[top30, ]
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(expr_data)
pheatmap(heatmap_data, annotation_col = annotation_col, show_colnames = FALSE,
         scale = "row", main = "Top 30 differentially expressed genes, AD vs control")

# --- Functional enrichment (GO Biological Process) ---
gene_list <- unique(sig_genes$gene_symbol)
gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]

ego <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                 ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05)

# --- Figure 3: GO enrichment dot plot ---
dotplot(ego, showCategory = 15) +
  ggtitle("Top enriched biological processes, AD DEGs")

# --- Export gene list for drug-target mining (Python) ---
write.csv(data.frame(gene = gene_list), "gene_list_for_drug_mining.csv",
          row.names = FALSE)

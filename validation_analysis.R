# ============================================================
# Independent validation of AD differential expression findings
# Dataset: GSE63060 (AddNeuroMed Cohort, peripheral blood)
# Compares against significant genes from GSE5281 (brain tissue)
# Requires: sig_genes (from differential_expression_analysis.R)
#           already loaded in the R session
# ============================================================

library(GEOquery)
library(limma)

# --- Load dataset ---
gse2 <- getGEO("GSE63060", GSEMatrix = TRUE)
expr_data2 <- exprs(gse2[[1]])
metadata2 <- pData(gse2[[1]])

# --- Filter to AD vs control only (drop MCI) ---
keep <- metadata2$`status:ch1` %in% c("AD", "CTL")
expr_data2 <- expr_data2[, keep]
metadata2 <- metadata2[keep, ]
group2 <- factor(ifelse(metadata2$`status:ch1` == "AD", "AD", "control"))
table(group2)  # AD: 145, control: 104

# --- Check transform: this dataset is already log2-scale, no transform needed ---
summary(as.vector(expr_data2))

# --- Differential expression ---
design2 <- model.matrix(~group2)
fit2 <- lmFit(expr_data2, design2)
fit2 <- eBayes(fit2)
results2 <- topTable(fit2, coef = 2, number = Inf, adjust.method = "BH")

sig_genes2 <- results2[results2$adj.P.Val < 0.05 & abs(results2$logFC) > 1, ]
nrow(sig_genes2)  # 8

# --- Map probe IDs to gene symbols (Illumina GPL6947 annotation) ---
gpl2 <- getGEO("GPL6947", AnnotGPL = TRUE)
platform_data2 <- Table(gpl2)
id_to_symbol2 <- setNames(platform_data2$`Gene symbol`, platform_data2$ID)
sig_genes2$gene_symbol <- id_to_symbol2[rownames(sig_genes2)]

write.csv(sig_genes2, "sig_genes_GSE63060_validation.csv")

# --- Compare against primary dataset (GSE5281) significant genes ---
# sig_genes must already be loaded from differential_expression_analysis.R
overlap <- intersect(sig_genes2$gene_symbol, sig_genes$gene_symbol)
overlap <- overlap[overlap != ""]
overlap        # "NDUFA1"
length(overlap) # 1

# NDUFA1 (mitochondrial Complex I subunit) replicates across brain (GSE5281)
# and blood (GSE63060), and is mechanistically consistent with the
# aerobic/cellular respiration pathways identified via GO enrichment
# in the primary analysis.

# --- Hypergeometric test: is the overlap more than expected by chance? ---
# Shared background: genes measured on both microarray platforms
background <- intersect(unique(platform_data$`Gene symbol`),
                         unique(platform_data2$`Gene symbol`))
background <- background[background != "" & !is.na(background)]
N <- length(background)                                              # 17,201
K <- length(intersect(unique(sig_genes$gene_symbol), background))     # 1,758 (GSE5281 sig genes on background)
n <- length(intersect(unique(sig_genes2$gene_symbol), background))    # 2 (GSE63060 sig genes on background)
x <- length(overlap[overlap != ""])                                   # 1 (NDUFA1)

p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
p_value  # 0.194 — not statistically significant

# Result: with only 2 significant blood genes landing on the shared
# background, this test is underpowered. p = 0.194 means the overlap
# cannot be distinguished from chance with this sample size — this is
# reported as an inconclusive validation, not a confirmed replication.
# The biological plausibility of NDUFA1 (mitochondrial Complex I,
# consistent with the respiration pathways from GO enrichment) stands
# independently of this statistical test.

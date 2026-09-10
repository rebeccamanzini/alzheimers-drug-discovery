# Identifying Druggable Targets in Alzheimer's Disease Through Independent Computational Transcriptomic Analysis

**Rebecca Manzini** · Independent project · September 2026

## Abstract

Alzheimer's disease (AD) remains a leading cause of dementia with limited disease-modifying treatment options. This project applies an independent, purely computational pipeline to public transcriptomic data to test whether unsupervised differential expression analysis can recover established AD biology and surface genes with existing pharmacological relevance. Using microarray data from 161 post-mortem brain samples ([GEO accession GSE5281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281)), differential expression analysis identified 2,972 significantly dysregulated probes between AD and control tissue. Gene Ontology enrichment of these genes converged on synaptic vesicle transport and mitochondrial/metabolic processes, both established mechanisms in AD pathology. Cross-referencing the top significant genes against the [Drug Gene Interaction Database (DGIdb)](https://dgidb.org) returned 449 gene-drug interactions, including SNAP25, a gene central to the synaptic vesicle cycle pathway identified by enrichment analysis and a known target of botulinum toxin. This convergence across independent analytical steps — statistical, pathway-level, and pharmacological — suggests the pipeline captured genuine disease-relevant signal rather than statistical noise.

## Research Question

Can differential gene expression analysis of AD brain tissue identify dysregulated genes with existing pharmacological relevance, connecting a purely data-driven transcriptomic signature to potential therapeutic targets?

## Dataset

| | |
|---|---|
| Source | Gene Expression Omnibus (GEO), accession **GSE5281** |
| Type | Affymetrix microarray, human post-mortem brain tissue |
| Platform | GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array) |
| Groups | Alzheimer's disease (AD) vs. healthy control |
| Sample size | 161 total (87 AD, 74 control) |
| Brain regions | Entorhinal cortex, hippocampus, middle temporal gyrus, posterior cingulate, superior frontal gyrus, primary visual cortex |

## Methods

1. **Preprocessing** — raw expression values log2-transformed (log2(x+1))
2. **Differential expression** — limma (R/Bioconductor), moderated t-statistics via eBayes
3. **Multiple testing correction** — Benjamini-Hochberg (BH)
4. **Significance thresholds** — adjusted p < 0.05, |log2 fold change| > 1
5. **Gene annotation** — probe-to-symbol mapping via GPL570 platform annotation
6. **Functional enrichment** — GO Biological Process via clusterProfiler
7. **Drug-target mining** — top 50 significant genes queried against DGIdb (GraphQL API)

**Tools:** R (limma, GEOquery, clusterProfiler, org.Hs.eg.db, ggplot2, pheatmap), Python (pandas, requests)

## Results

### Differential expression

2,972 probes met significance thresholds. The volcano plot below shows greater spread among downregulated genes than upregulated genes — consistent with the AD literature, where downregulation typically predominates, reflecting broad neuronal loss and reduced synaptic/metabolic activity.

![Volcano plot](figures/volcano_plot_AD_vs_control.png)

*Figure 1. Volcano plot of differential gene expression, AD vs. control (GSE5281). Red points met significance thresholds (adjusted p < 0.05, |log2FC| > 1); n = 2,972 significant probes.*

Unsupervised hierarchical clustering of the top 30 differentially expressed genes largely recovered the known AD/control grouping without being given group labels, indicating a strong underlying transcriptional signal.

![Heatmap](figures/heatmap_top30_genes.png)

*Figure 2. Heatmap of top 30 differentially expressed genes by adjusted p-value, row-scaled expression, with hierarchical clustering of genes and samples.*

### Functional enrichment

GO Biological Process enrichment identified vesicle-mediated transport in synapse, synaptic vesicle cycle, aerobic/cellular respiration, and regulation of neuron projection development among the top enriched terms — aligning with two of the most established mechanisms in AD pathology: synaptic dysfunction and mitochondrial/metabolic impairment.

![GO enrichment dot plot](figures/go_enrichment_dotplot.png)

*Figure 3. Dot plot of top enriched GO Biological Process terms. Dot size = gene count; colour = BH-adjusted p-value; x-axis = gene ratio.*

### Drug-target mining

Querying the top 50 significant genes against DGIdb returned **449 gene-drug interactions**, with 19/50 genes having at least one known interaction. The standout result was **SNAP25**, targeted by multiple botulinum toxin compounds via a cleavage mechanism — mechanistically consistent with botulinum toxin's established mode of action, and directly connecting back to the synaptic vesicle cycle pathway identified independently through GO enrichment. Two separate analytical steps converging on the same biology.

*Note: interaction scores reflect DGIdb's internal weighting, not clinical efficacy. Approval status of individual compounds should be verified separately before any therapeutic interpretation.*

## Discussion

This analysis demonstrates that an unsupervised, purely computational approach applied to public transcriptomic data can recover established AD biology and identify existing pharmacological entry points without disease-specific assumptions built into the pipeline. The convergence between differential expression, pathway enrichment, and drug-target mining around synaptic mechanisms — most notably the SNAP25/botulinum toxin result — strengthens confidence that the pipeline captured genuine biological signal.

## Limitations and Future Directions

- **Single dataset, single platform** — validation against an independent AD dataset (e.g. GSE63060) would strengthen confidence
- **Interaction scores ≠ clinical efficacy** — DGIdb scores reflect database weighting, not approval status or therapeutic validity
- **Correlative, not causal** — functional validation (e.g. in vitro knockdown/rescue) would be needed to establish causal gene relevance
- **Brain region not modelled as a covariate** — a region-stratified re-analysis could reveal signal currently averaged out

## Repository Structure

```
├── README.md
├── scripts/
│   ├── differential_expression_analysis.R   # full R pipeline
│   └── drug_target_mining.py                # DGIdb query script
└── figures/
    ├── volcano_plot_AD_vs_control.png
    ├── heatmap_top30_genes.png
    └── go_enrichment_dotplot.png
```

## References

- Gene Expression Omnibus, accession GSE5281. National Center for Biotechnology Information.
- Ritchie, M.E. et al. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*.
- Yu, G. et al. (2012) clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS*.
- Freshour, S.L. et al. (2021) Integration of the Drug-Gene Interaction Database (DGIdb 4.0) with open crowdsource efforts. *Nucleic Acids Research*.

# RNA-seq Analysis of Catuaí and Costa Rica 95

This project analyzes the transcriptional response of two *Coffea arabica* cultivars, Catuaí and Costa Rica 95 (CR95), to *Xylella fastidiosa*, using saline-treated samples as controls.

The analysis includes read preprocessing, differential expression, GO enrichment, KEGG annotation, co-expression networks, and promoter motif discovery.

## Project Files

| Script | Contents |
| --- | --- |
| `Analisis_diferencial_0.1(1).R` | DESeq2, PCA, heatmaps, KEGG annotation, and gene lists for promoter analysis. |
| `CAMERA_0.1(1).R` | GO enrichment with CAMERA, bar plots, UpSet plots, and GO term comparisons. |
| `WGCNA_bien(1).R` | Co-expression network, modules, and export of tables to Cytoscape with CR95 DEGs. |

## 1. Read Processing in Galaxy

Read quality was assessed using **FastQC and MultiQC**. Reads were processed with **Trimmomatic** to remove TruSeq3 paired-end adapters, the first 12 nucleotides, and trailing bases with quality scores below 30.

```text
FastQC / MultiQC → Trimmomatic → FastQC / MultiQC

HEADCROP:12
TRAILING:30
```

The processed reads were aligned to the *C. arabica* ET-39 HiFi reference genome (`GCF_036785885.1`) using **HISAT2** and quantified with **featureCounts** to obtain the count matrix.

**Related figures:**

- **Figure 1:** quality, base composition, and adapter content before and after processing.
- **Supplementary Figures S8–S9:** additional quality results.
- **Supplementary Figure S10a:** HISAT2 alignment.
- **Supplementary Figure S10b:** read assignment with featureCounts.

This stage was performed in Galaxy.

## 2. Differential Expression Analysis in R — DESeq2

**Script:** `Analisis_diferencial_0.1(1).R`  
**Section:** `8. RUN DESeq2 FOR EACH CULTIVAR`

The analysis was performed separately for Catuaí and CR95, comparing `xylella` against `saline`. Genes with zero counts across all samples of each cultivar were removed.

```r
metadata_cv$treatment <- factor(
  metadata_cv$treatment,
  levels = c("saline", "xylella")
)

counts_cv <- counts_cv[
  rowSums(counts_cv) > 0,
  ,
  drop = FALSE
]

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_cv),
  colData = metadata_cv,
  design = ~ treatment
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("treatment", "xylella", "saline"),
  alpha = 0.1
)
```

Genes were considered significant at **padj < 0.10**, with no minimum log2FoldChange magnitude threshold.

After converting `res` into a data frame, genes were classified by direction:

```r
res$Direction <- "Not significant"

res$Direction[
  !is.na(res$padj) &
    res$padj < 0.1 &
    res$log2FoldChange > 0
] <- "Up"

res$Direction[
  !is.na(res$padj) &
    res$padj < 0.1 &
    res$log2FoldChange < 0
] <- "Down"
```

**Related outputs:**

- **Table 1:** summary of analyzed, significant, upregulated, and downregulated genes.
- **Supplementary Figure S1:** differential expression analysis workflow.

The summary is saved as `DESeq2_summary_0.1.txt`.

## 3. PCA — Figure 2

**Script:** `Analisis_diferencial_0.1(1).R`  
**Block:** `PCA - ALL GENES`

VST was applied, and PCA was calculated using all retained genes from the 15 samples.

```r
dds_pca <- DESeqDataSetFromMatrix(
  countData = round(counts_pca),
  colData = metadata_pca,
  design = ~ cultivar + treatment
)

vsd_pca <- vst(
  dds_pca,
  blind = TRUE
)

pca_matrix <- t(assay(vsd_pca))

pca_all <- prcomp(
  pca_matrix,
  center = TRUE,
  scale. = FALSE
)
```

The plot was created with **ggplot2**, using color to represent cultivar and shape to represent treatment.

**Related figure:** Figure 2.  
**Output file:** `PCA_all_genes_Catuai_CR95_labeled_purple.pdf`.

## 4. Expression Heatmaps — Figure 3a-b

**Script:** `Analisis_diferencial_0.1(1).R`  
**Function:** `plot_DESeq2_heatmap()`

DEGs from each cultivar were plotted using **pheatmap**. Normalized counts were log2-transformed and scaled to Z-scores for each gene.

```r
heatmap_mat <- log2(heatmap_mat + 1)

# Después de eliminar las filas sin variación:
heatmap_z <- t(scale(t(heatmap_mat)))
```

The following calls generate the heatmaps for both cultivars:

```r
plot_DESeq2_heatmap(
  cv = "Catuai",
  deseq_results = deseq_results,
  norm_counts_list = norm_counts_list,
  metadata = metadata,
  outpathHeat = outpathHeat
)

plot_DESeq2_heatmap(
  cv = "CR95",
  deseq_results = deseq_results,
  norm_counts_list = norm_counts_list,
  metadata = metadata,
  outpathHeat = outpathHeat
)
```

**Related figures:**

- **Figure 3a:** Catuaí heatmap.
- **Figure 3b:** CR95 heatmap.

**Output files:**

- `Heatmap_Catuai_significant_genes_0.1.pdf`
- `Heatmap_CR95_significant_genes_0.1.pdf`

## 5. GO Enrichment — CAMERA

**Script:** `CAMERA_0.1(1).R`  
**Section:** `13. LOOP CAMERA`

**edgeR** was used to filter low-expression genes and normalize counts. **limma** was then used for voom transformation and CAMERA analysis.

The procedure was performed separately for each cultivar.

```r
dge <- DGEList(
  counts = counts_cv,
  group = treatment
)

keep <- filterByExpr(
  dge,
  group = treatment
)

dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]

dge <- calcNormFactors(dge)

design <- model.matrix(~ treatment)

v <- voom(
  dge,
  design,
  plot = FALSE
)

# idx contiene las posiciones de los genes de cada conjunto GO.
res <- camera(
  v,
  index = idx,
  design = design,
  contrast = 2
)
```

GO sets containing **at least five genes** after filtering were retained, and enrichment was considered significant at **FDR < 0.10**.

Contrast 2 corresponds to the effect of `xylella` relative to `saline`.

### Enrichment Bar Plots — Figure 4a-b

The `plot_camera_barplot()` function displays the 15 significant GO terms with the lowest FDR for each cultivar.

```r
plot_camera_barplot(
  sig_results[["Catuai"]],
  "Catuai",
  genes_analyzed[["Catuai"]],
  top_n = 15
)

plot_camera_barplot(
  sig_results[["CR95"]],
  "CR95",
  genes_analyzed[["CR95"]],
  top_n = 15
)
```

**Related figures:**

- **Figure 4a:** Catuaí GO terms.
- **Figure 4b:** CR95 GO terms.

### GO Term Comparison — Supplementary Figure S12

This figure was created by a teammate in **R using ggplot2**, based on the CAMERA results and the following code.

The 10 significant GO terms with the lowest FDR were selected for each cultivar. Their results were then retrieved for both cultivars for comparison, including terms that were not significant in one of them.

Point size indicates the number of genes in the GO set, color represents `-log10(FDR)`, and shape indicates enrichment direction: diamonds for Up and stars for Down.

```r
# 15. Comparative dotplot - GO terms across cultivars (CAMERA)
# ============================================================

top_terms_per_cv <- lapply(names(sig_results), function(cv) {
  res <- sig_results[[cv]]
  if (nrow(res) == 0) return(NULL)
  res_top <- head(res[order(res$FDR), ], 10)
  res_top$cultivar <- cv
  res_top
})

comparison_df <- do.call(rbind, top_terms_per_cv)

if (is.null(comparison_df) || nrow(comparison_df) == 0) {

  cat("No significant GO terms in either cultivar -- comparison plot skipped.\n")

} else {

  all_selected_terms <- unique(comparison_df$GOID)

  full_comparison <- do.call(rbind, lapply(names(camera_results), function(cv) {
    res <- camera_results[[cv]]
    res_sub <- res[res$GOID %in% all_selected_terms, ]
    if (nrow(res_sub) == 0) return(NULL)
    res_sub$cultivar <- cv
    res_sub
  }))

  term_labels <- unique(comparison_df[, c("GOID", "TERM")])

  full_comparison <- merge(
    full_comparison,
    term_labels,
    by = "GOID",
    suffixes = c("", "_label")
  )

  p_comparison <- ggplot(
    full_comparison,
    aes(
      x = cultivar,
      y = TERM_label,
      size = NGenes,
      color = -log10(FDR),
      shape = Direction
    )
  ) +
    geom_point() +
    scale_color_gradient(low = "#9B59B6", high = "#E67E22") +
    scale_shape_manual(values = c(Up = 18, Down = 8)) +
    guides(size = guide_legend(override.aes = list(shape = 8))) +
    labs(
      title = "GO terms compared across cultivars (CAMERA)",
      x = "Cultivar",
      y = "GO term",
      size = "Genes in set",
      color = "-log10(FDR)",
      shape = "Direction"
    ) +
    theme_bw()

  print(p_comparison)

  ggsave(
    paste0(outpathcount, "comparison_dotplot.pdf"),
    p_comparison,
    width = 9,
    height = 7
  )
}
```

**Output file:** `comparison_dotplot.pdf`.

### DEG Comparison Using an UpSet Plot — Supplementary Figure S11

**Script:** `Analisis_diferencial_0.1(1).R`  
**Block:** `UPSET - DESeq2 SIGNIFICANT GENES`

The UpSet plot was created in **R using UpSetR** to compare significant DEGs (**padj < 0.10**) across four groups: Catuaí UP, Catuaí DOWN, CR95 UP, and CR95 DOWN.

The horizontal bars show the total number of genes in each group. The vertical bars show the number of genes in each exclusive intersection, identified by the connected dots below.

```r
# ============================================================
# UPSET - DESeq2 SIGNIFICANT GENES
#
# ALL SIGNIFICANT GENES
# padj < 0.1
#
# Catuai UP
# Catuai DOWN
# CR95 UP
# CR95 DOWN
# ============================================================
# ============================================================
# UPSET PLOT
# Shared significant DEGs between cultivars and directions
# padj < 0.1
# ============================================================

library(UpSetR)

# Crear lista de conjuntos
upset_sets <- list(
  "Catuai UP"   = Catuai_UP,
  "Catuai DOWN" = Catuai_DOWN,
  "CR95 UP"     = CR95_UP,
  "CR95 DOWN"   = CR95_DOWN
)

# Convertir a formato UpSetR
upset_data <- fromList(upset_sets)

# Crear UpSet plot
upset(
  upset_data,
  sets = c(
    "Catuai UP",
    "Catuai DOWN",
    "CR95 UP",
    "CR95 DOWN"
  ),
  keep.order = TRUE,
  order.by = "freq",
  nintersects = NA,
  number.angles = 0,
  point.size = 3,
  line.size = 1,
  main.bar.color = "#7A35F0",
  sets.bar.color = "#1DB5A6",
  matrix.color = "#6D28D9",
  mainbar.y.label = "Number of shared genes",
  sets.x.label = "Number of significant genes"
)

pdf(
  file.path(
    outpathDE,
    "UpSet_DEGs_Catuai_CR95_0.1.pdf"
  ),
  width = 11,
  height = 8
)

upset(
  upset_data,
  sets = c(
    "Catuai UP",
    "Catuai DOWN",
    "CR95 UP",
    "CR95 DOWN"
  ),
  keep.order = TRUE,
  order.by = "freq",
  nintersects = NA,
  number.angles = 0,
  point.size = 3,
  line.size = 1,
  main.bar.color = "#7A35F0",
  sets.bar.color = "#1DB5A6",
  matrix.color = "#6D28D9",
  mainbar.y.label = "Number of shared genes",
  sets.x.label = "Number of significant genes"
)

dev.off()
```

**Output file:** `UpSet_DEGs_Catuai_CR95_0.1.pdf`.

## 6. KEGG Annotation and Pathway Heatmaps

**Script:** `Analisis_diferencial_0.1(1).R`  
**Blocks:** KEGG annotation and sections `20–21`.

DEGs were matched to their functional annotation to retrieve KO identifiers. Files were prepared for **KEGG Mapper Color**, using red for upregulated genes and blue for downregulated genes.

```r
kegg_color_CR95 <- kegg_CR95 %>%
  filter(
    !is.na(KEGG_ko),
    KEGG_ko != "-",
    KEGG_ko != ""
  ) %>%
  separate_rows(
    KEGG_ko,
    sep = "[,;]"
  ) %>%
  mutate(
    KEGG_ko = trimws(KEGG_ko),
    KO = sub("^ko:", "", KEGG_ko),
    Color = ifelse(log2FoldChange > 0, "red", "blue")
  )
```

The same procedure is applied to Catuaí.

**Supplementary Figures S2a-b to S7a-b:** pathway maps generated using KEGG Mapper Color.

### Pathway log2FoldChange Heatmaps

The `comp_heatmap()` function generates heatmaps using **pheatmap**, with one column per cultivar. Genes that are not significant in a cultivar appear in gray.

```r
for (ko in c(
  "ko04626",
  "ko04016",
  "ko04075",
  "ko00270",
  "ko00195",
  "ko03010"
)) {
  comp_heatmap(
    ko,
    file.path(
      outpathDE,
      paste0("heatmap_", ko, "_Catuai_vs_CR95.pdf")
    )
  )
}
```

**Related figures:**

- **Figure 5a:** MAPK signaling pathway – plant (`ko04016`).
- **Figure 5b:** plant-pathogen interaction (`ko04626`).
- **Supplementary Figures S13–S16:** additional pathway heatmaps.

Assignments were checked against NCBI annotations, and selected proteins were compared against *Arabidopsis thaliana* using **BLASTp**. These external searches are not included in the three scripts.

## 7. Co-expression Network — WGCNA and Cytoscape

**Script:** `WGCNA_bien(1).R`

Genes with fewer than 10 total counts were removed, VST was applied, and the 5,000 most variable genes across all samples were selected.

```r
keep <- rowSums(counts) >= 10

counts <- counts[
  keep,
  ,
  drop = FALSE
]

# Después de crear dds con design = ~ 1:
vsd <- vst(dds, blind = TRUE)

vst_mat <- assay(vsd)

gene_var <- apply(
  vst_mat,
  1,
  var
)

n_top_genes <- min(
  5000,
  length(gene_var)
)

top_genes <- names(
  sort(gene_var, decreasing = TRUE)
)[1:n_top_genes]
```

Sections `19–24` select the power, calculate the network, and group genes into modules.

```r
adj <- adjacency(
  datExpr,
  power = softPower,
  type = "unsigned"
)

TOM <- TOMsimilarity(
  adj,
  TOMType = "unsigned",
  verbose = 0
)

dissTOM <- 1 - TOM

geneTree <- hclust(
  as.dist(dissTOM),
  method = "average"
)
```

For the complete network, the 10 strongest TOM connections per gene were retained. Tables were exported to **Cytoscape**, where modules and log2FoldChange values were visualized.

### Greenyellow Module — Figure 6

The section `53. greenyellow MODULE + DESeq2` prepares the greenyellow module table with CR95 results.

```r
greenyellow_table <- cytoscape_DE[
  cytoscape_DE$module == "greenyellow",
  ,
  drop = FALSE
]
```

**Output files:**

- `Cytoscape_ALL_CLUSTERS_DESeq2_CR95_0.1.txt`
- `greenyellow_ALL_genes_DESeq2_CR95_0.1.txt`

**Related figures:**

- **Figure 6:** greenyellow module with CR95 DEGs.
- **Supplementary Figures S17a-b:** complete networks and DEG distributions for both cultivars.

The provided script includes the integration of CR95 DEGs; the integration of Catuaí DEGs is not included in this version.

## 8. Promoters and De Novo Motif Discovery

**Script:** `Analisis_diferencial_0.1(1).R`  
**Block:** `EXPORTAR LISTAS DE GENES PARA ANALISIS DE PROMOTORES`

Lists of significant UP and DOWN genes were prepared in R for each cultivar. For example:

```r
CR95_UP <- deseq_results[["CR95"]] %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.1,
    log2FoldChange > 0
  ) %>%
  dplyr::pull(Geneid) %>%
  unique()
```

The four lists were exported as text files to continue the analysis outside R.

The subsequent steps were:

1. **SAMtools:** index the reference genome.
2. **BEDTools:** extract 1,000 bp upstream of each gene, accounting for gene orientation.
3. **MEME:** search for up to five motifs, each 8–20 bp long.
4. **Tomtom and JASPAR:** compare motifs with known profiles, considering matches significant at **q-value < 0.05**.

**Related table:** Table 2, showing selected de novo motifs and their matches to transcription factor binding profiles.

The promoter extraction and MEME/Tomtom commands are not included in the attached scripts.

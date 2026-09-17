# ============================================================
# DESeq2 - DIFFERENTIAL EXPRESSION
# Catuai y CR95
# Xylella vs saline
#
# SIGNIFICANCE CRITERION:
# padj < 0.1
#
# UP   = log2FoldChange > 0
# DOWN = log2FoldChange < 0
#
# NO minimum log2FoldChange threshold
# ============================================================


# ============================================================
# 1. WORKING DIRECTORY
# ============================================================

setwd("~/Desktop/Changed padj 0.1")


# ============================================================
# 2. OUTPUT FOLDERS
# ============================================================

outpathDE <- file.path(
  "~/Desktop/Changed padj 0.1",
  "DESeq2_RESULTS"
)

dir.create(
  outpathDE,
  showWarnings = FALSE,
  recursive = TRUE
)


outpathHeat <- file.path(
  "~/Desktop/Changed padj 0.1",
  "HEATMAPS_DESeq2"
)

dir.create(
  outpathHeat,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 3. INSTALL / LOAD PACKAGES
# ============================================================

options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}


if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}


if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}


if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}


if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}


library(DESeq2)
library(pheatmap)
library(UpSetR)
library(dplyr)
library(tidyr)


# ============================================================
# 4. LOAD HISAT2 COUNTS
# ============================================================

counts <- read.table(
  "Matrix_hisat2 copy 2.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)


# ============================================================
# 5. LOAD METADATA
# ============================================================

metadata <- read.table(
  "metadata_copy copy.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


# ============================================================
# 6. FIX SAMPLE NAMES
# ============================================================

# Corregir Ca7x si aparece como Cax7

colnames(counts)[
  colnames(counts) == "Cax7"
] <- "Ca7x"


# Confirmar que todas las muestras del metadata existen
# en la matriz de counts

stopifnot(
  all(
    metadata$sample %in%
      colnames(counts)
  )
)


# Ordenar counts exactamente igual que metadata

counts <- counts[
  ,
  metadata$sample,
  drop = FALSE
]


stopifnot(
  identical(
    colnames(counts),
    metadata$sample
  )
)


# ============================================================
# 7. CHECK EXPERIMENTAL DESIGN
# ============================================================

print(
  table(
    metadata$cultivar,
    metadata$treatment
  )
)


# Esperado:
#
#          saline xylella
# Catuai       3       5
# CR95         4       3


# ============================================================
# 8. RUN DESeq2 FOR EACH CULTIVAR
# ============================================================

cultivars <- c(
  "Catuai",
  "CR95"
)


# Aquí guardaremos resultados DESeq2

deseq_results <- list()


# Aquí guardaremos normalized counts

norm_counts_list <- list()


# Aquí guardaremos objetos DESeq2

dds_list <- list()


# Tabla resumen

DE_summary <- data.frame()


for (cv in cultivars) {
  
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    "Running DESeq2 for:",
    cv,
    "\n"
  )
  
  cat(
    "====================================\n"
  )
  
  
  # ==========================================================
  # 8.1 SELECT CULTIVAR
  # ==========================================================
  
  metadata_cv <- metadata[
    metadata$cultivar == cv,
    ,
    drop = FALSE
  ]
  
  
  counts_cv <- counts[
    ,
    metadata_cv$sample,
    drop = FALSE
  ]
  
  
  # ==========================================================
  # 8.2 SET SALINE AS REFERENCE
  # ==========================================================
  
  metadata_cv$treatment <- factor(
    metadata_cv$treatment,
    levels = c(
      "saline",
      "xylella"
    )
  )
  
  
  rownames(metadata_cv) <-
    metadata_cv$sample
  
  
  print(
    table(
      metadata_cv$treatment
    )
  )
  
  
  # ==========================================================
  # 8.3 REMOVE GENES WITH ZERO COUNTS
  # ==========================================================
  
  counts_cv <- counts_cv[
    rowSums(counts_cv) > 0,
    ,
    drop = FALSE
  ]
  
  
  cat(
    "Genes analyzed:",
    nrow(counts_cv),
    "\n"
  )
  
  
  # ==========================================================
  # 8.4 CREATE DESeq2 OBJECT
  # ==========================================================
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_cv),
    colData = metadata_cv,
    design = ~ treatment
  )
  
  
  # ==========================================================
  # 8.5 RUN DESeq2
  # ==========================================================
  
  dds <- DESeq(
    dds
  )
  
  
  # Guardar objeto
  dds_list[[cv]] <- dds
  
  
  # ==========================================================
  # 8.6 XYELLA VS SALINE
  #
  # POSITIVE log2FC:
  # higher expression in Xylella
  #
  # NEGATIVE log2FC:
  # lower expression in Xylella
  # ==========================================================
  
  res <- results(
    dds,
    contrast = c(
      "treatment",
      "xylella",
      "saline"
    ),
    alpha = 0.1
  )
  
  
  # Convertir a dataframe
  
  res <- as.data.frame(
    res
  )
  
  
  # Agregar Gene ID
  
  res$Geneid <-
    rownames(res)
  
  
  # Reordenar columnas
  
  res <- res[
    ,
    c(
      "Geneid",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj"
    )
  ]
  
  
  # ==========================================================
  # 8.7 CLASSIFY DIRECTION
  #
  # padj < 0.1
  # ==========================================================
  
  res$Direction <-
    "Not significant"
  
  
  # UP
  
  res$Direction[
    !is.na(res$padj) &
      res$padj < 0.1 &
      res$log2FoldChange > 0
  ] <- "Up"
  
  
  # DOWN
  
  res$Direction[
    !is.na(res$padj) &
      res$padj < 0.1 &
      res$log2FoldChange < 0
  ] <- "Down"
  
  
  # ==========================================================
  # 8.8 SORT BY padj
  # ==========================================================
  
  res <- res[
    order(
      res$padj,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  
  # Guardar en lista
  
  deseq_results[[cv]] <-
    res
  
  
  # ==========================================================
  # 8.9 SAVE ALL GENES
  # ==========================================================
  
  write.table(
    res,
    file = file.path(
      outpathDE,
      paste0(
        "DESeq2_",
        cv,
        "_xylella_vs_saline_ALL_0.1.txt"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 8.10 SELECT SIGNIFICANT GENES
  #
  # padj < 0.1
  # ==========================================================
  
  sig <- res[
    !is.na(res$padj) &
      res$padj < 0.1,
    ,
    drop = FALSE
  ]
  
  
  # ==========================================================
  # 8.11 SAVE SIGNIFICANT GENES
  # ==========================================================
  
  write.table(
    sig,
    file = file.path(
      outpathDE,
      paste0(
        "DESeq2_",
        cv,
        "_xylella_vs_saline_SIGNIFICANT_0.1.txt"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 8.12 NORMALIZED COUNTS
  # ==========================================================
  
  norm_counts <- counts(
    dds,
    normalized = TRUE
  )
  
  
  norm_counts_list[[cv]] <-
    norm_counts
  
  
  write.table(
    norm_counts,
    file = file.path(
      outpathDE,
      paste0(
        "Normalized_counts_",
        cv,
        "_0.1.txt"
      )
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  
  # ==========================================================
  # 8.13 PRINT SUMMARY
  # ==========================================================
  
  cat(
    "\nDifferential expression summary for",
    cv,
    "\n"
  )
  
  
  print(
    table(
      res$Direction
    )
  )
  
  
  # ==========================================================
  # 8.14 ADD TO SUMMARY TABLE
  # ==========================================================
  
  DE_summary <- rbind(
    DE_summary,
    
    data.frame(
      Cultivar = cv,
      
      Genes_analyzed =
        nrow(res),
      
      Significant =
        sum(
          res$Direction !=
            "Not significant"
        ),
      
      Up =
        sum(
          res$Direction ==
            "Up"
        ),
      
      Down =
        sum(
          res$Direction ==
            "Down"
        )
    )
  )
}


# ============================================================
# 9. FINAL DESeq2 SUMMARY
# ============================================================

print(
  DE_summary
)


write.table(
  DE_summary,
  file = file.path(
    outpathDE,
    "DESeq2_summary_0.1.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 10. SHOW SIGNIFICANT COUNTS
# ============================================================

cat(
  "\n====================================\n"
)

cat(
  "FINAL SUMMARY padj < 0.1\n"
)

cat(
  "====================================\n"
)


cat(
  "\nCATUAI:\n"
)

print(
  table(
    deseq_results[["Catuai"]]$Direction
  )
)


cat(
  "\nCR95:\n"
)

print(
  table(
    deseq_results[["CR95"]]$Direction
  )
)

## heatmaps
# ============================================================
# HEATMAPS - SIGNIFICANT GENES
#
# Separate heatmaps:
# 1) Catuai
# 2) CR95
#
# Significant genes:
# padj < 0.1
#
# Values:
# log2(normalized counts + 1)
# Row-scaled Z-scores
#
# Colors:
# Blue  = lower relative expression
# White = average relative expression
# Red   = higher relative expression
# ============================================================


# ============================================================
# HEATMAP FUNCTION
# ============================================================

plot_DESeq2_heatmap <- function(
    cv,
    deseq_results,
    norm_counts_list,
    metadata,
    outpathHeat
) {
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    "Creating heatmap for:",
    cv,
    "\n"
  )
  
  cat(
    "====================================\n"
  )
  
  
  # ----------------------------------------------------------
  # 1. SELECT SIGNIFICANT GENES
  # padj < 0.1
  # ----------------------------------------------------------
  
  sig_genes <- deseq_results[[cv]]$Geneid[
    
    !is.na(
      deseq_results[[cv]]$padj
    ) &
      
      deseq_results[[cv]]$padj < 0.1
    
  ]
  
  
  cat(
    "Significant genes:",
    length(sig_genes),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Safety check
  # ----------------------------------------------------------
  
  if (length(sig_genes) == 0) {
    
    warning(
      paste(
        "No significant genes found for",
        cv
      )
    )
    
    return(
      invisible(NULL)
    )
  }
  
  
  # ----------------------------------------------------------
  # 2. GET NORMALIZED COUNTS
  # ----------------------------------------------------------
  
  norm_counts_cv <-
    norm_counts_list[[cv]]
  
  
  cat(
    "Normalized count matrix:",
    nrow(norm_counts_cv),
    "genes x",
    ncol(norm_counts_cv),
    "samples\n"
  )
  
  
  # ----------------------------------------------------------
  # 3. KEEP SIGNIFICANT GENES
  # ----------------------------------------------------------
  
  heatmap_mat <- norm_counts_cv[
    
    rownames(norm_counts_cv) %in%
      sig_genes,
    
    ,
    
    drop = FALSE
  ]
  
  
  cat(
    "Significant genes found in normalized counts:",
    nrow(heatmap_mat),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Safety check
  # ----------------------------------------------------------
  
  if (nrow(heatmap_mat) == 0) {
    
    warning(
      paste(
        "No significant genes could be matched for",
        cv
      )
    )
    
    return(
      invisible(NULL)
    )
  }
  
  
  # ----------------------------------------------------------
  # 4. LOG2 TRANSFORMATION
  # ----------------------------------------------------------
  
  heatmap_mat <-
    log2(
      heatmap_mat + 1
    )
  
  
  # ----------------------------------------------------------
  # 5. REMOVE ZERO-VARIANCE GENES
  #
  # Necessary before calculating Z-score
  # ----------------------------------------------------------
  
  gene_sd <- apply(
    heatmap_mat,
    1,
    sd
  )
  
  
  heatmap_mat <- heatmap_mat[
    
    !is.na(gene_sd) &
      gene_sd > 0,
    
    ,
    
    drop = FALSE
  ]
  
  
  cat(
    "Genes after removing zero-variance rows:",
    nrow(heatmap_mat),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Safety check
  # ----------------------------------------------------------
  
  if (nrow(heatmap_mat) < 2) {
    
    warning(
      paste(
        "Too few genes remain to create heatmap for",
        cv
      )
    )
    
    return(
      invisible(NULL)
    )
  }
  
  
  # ----------------------------------------------------------
  # 6. ROW Z-SCORE
  # ----------------------------------------------------------
  
  heatmap_z <- t(
    scale(
      t(heatmap_mat)
    )
  )
  
  
  # ----------------------------------------------------------
  # Remove any row that became NA after scaling
  # ----------------------------------------------------------
  
  heatmap_z <- heatmap_z[
    
    apply(
      heatmap_z,
      1,
      function(x) all(is.finite(x))
    ),
    
    ,
    
    drop = FALSE
  ]
  
  
  cat(
    "Genes plotted:",
    nrow(heatmap_z),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 7. METADATA FOR THIS CULTIVAR
  # ----------------------------------------------------------
  
  metadata_cv <- metadata[
    metadata$cultivar == cv,
    ,
    drop = FALSE
  ]
  
  
  rownames(metadata_cv) <-
    metadata_cv$sample
  
  
  # ----------------------------------------------------------
  # Make sure samples appear in correct order
  # ----------------------------------------------------------
  
  metadata_cv <- metadata_cv[
    colnames(heatmap_z),
    ,
    drop = FALSE
  ]
  
  
  stopifnot(
    identical(
      rownames(metadata_cv),
      colnames(heatmap_z)
    )
  )
  
  
  # ----------------------------------------------------------
  # 8. SAMPLE ANNOTATION
  # ----------------------------------------------------------
  
  annotation_col <- data.frame(
    
    Treatment =
      metadata_cv$treatment
    
  )
  
  
  rownames(annotation_col) <-
    metadata_cv$sample
  
  
  # ----------------------------------------------------------
  # 9. TREATMENT COLORS
  #
  # Saline  = grey
  # Xylella = red
  # ----------------------------------------------------------
  
  annotation_colors <- list(
    
    Treatment = c(
      
      saline = "#BDBDBD",
      
      xylella = "#D73027"
      
    )
  )
  
  
  # ----------------------------------------------------------
  # 10. HEATMAP COLORS
  #
  # Blue -> White -> Red
  # ----------------------------------------------------------
  
  heat_colors <- colorRampPalette(
    
    c(
      "#2166AC",
      "#67A9CF",
      "#D1E5F0",
      "#FFFFFF",
      "#FDDBC7",
      "#EF8A62",
      "#B2182B"
    )
    
  )(100)
  
  
  # ----------------------------------------------------------
  # 11. OUTPUT FILE
  # ----------------------------------------------------------
  
  heatmap_file <- file.path(
    
    outpathHeat,
    
    paste0(
      "Heatmap_",
      cv,
      "_significant_genes_0.1.pdf"
    )
  )
  
  
  cat(
    "Saving heatmap to:\n",
    heatmap_file,
    "\n"
  )
  
  
  # ==========================================================
  # 12. CREATE AND SAVE HEATMAP
  #
  # pheatmap saves directly to the PDF
  # ==========================================================
  
  pheatmap(
    
    heatmap_z,
    
    color =
      heat_colors,
    
    annotation_col =
      annotation_col,
    
    annotation_colors =
      annotation_colors,
    
    cluster_rows =
      TRUE,
    
    cluster_cols =
      TRUE,
    
    show_rownames =
      FALSE,
    
    show_colnames =
      TRUE,
    
    border_color =
      NA,
    
    fontsize =
      10,
    
    fontsize_col =
      10,
    
    main = paste0(
      cv,
      " - Significant genes"
    ),
    
    filename =
      heatmap_file,
    
    width =
      8,
    
    height =
      10
  )
  
  
  # ----------------------------------------------------------
  # Confirm that file was created
  # ----------------------------------------------------------
  
  if (file.exists(heatmap_file)) {
    
    cat(
      "SUCCESS:",
      cv,
      "heatmap created.\n"
    )
    
    cat(
      "File size:",
      file.info(heatmap_file)$size,
      "bytes\n"
    )
    
  } else {
    
    warning(
      paste(
        "Heatmap file was NOT created for",
        cv
      )
    )
    
  }
}


# ============================================================
# CATUAI HEATMAP
# ============================================================

plot_DESeq2_heatmap(
  
  cv = "Catuai",
  
  deseq_results =
    deseq_results,
  
  norm_counts_list =
    norm_counts_list,
  
  metadata =
    metadata,
  
  outpathHeat =
    outpathHeat
)


# ============================================================
# CR95 HEATMAP
# ============================================================

plot_DESeq2_heatmap(
  
  cv = "CR95",
  
  deseq_results =
    deseq_results,
  
  norm_counts_list =
    norm_counts_list,
  
  metadata =
    metadata,
  
  outpathHeat =
    outpathHeat
)


#PCA
# ============================================================
# PCA - ALL GENES
# Catuai + CR95
# Color = cultivar
# Shape = treatment
# Sample names shown above each point
# ============================================================


# ============================================================
# 1. INSTALL / LOAD GGPLOT2
# ============================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)


# ============================================================
# 2. PREPARE METADATA
# ============================================================

metadata_pca <- metadata

rownames(metadata_pca) <- metadata_pca$sample


# ============================================================
# 3. CHECK SAMPLE ORDER
# ============================================================

counts_pca <- counts[
  ,
  metadata_pca$sample,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(counts_pca),
    rownames(metadata_pca)
  )
)


# ============================================================
# 4. REMOVE GENES WITH ZERO COUNTS IN ALL SAMPLES
# ============================================================

counts_pca <- counts_pca[
  rowSums(counts_pca) > 0,
  ,
  drop = FALSE
]

cat(
  "Genes used before VST:",
  nrow(counts_pca),
  "\n"
)


# ============================================================
# 5. CREATE DESEQ2 OBJECT
# ============================================================

dds_pca <- DESeqDataSetFromMatrix(
  countData = round(counts_pca),
  colData = metadata_pca,
  design = ~ cultivar + treatment
)


# ============================================================
# 6. VARIANCE STABILIZING TRANSFORMATION
# ============================================================

vsd_pca <- vst(
  dds_pca,
  blind = TRUE
)


# ============================================================
# 7. PCA USING ALL GENES
# ============================================================

pca_matrix <- t(
  assay(vsd_pca)
)

pca_all <- prcomp(
  pca_matrix,
  center = TRUE,
  scale. = FALSE
)


# ============================================================
# 8. PERCENTAGE OF VARIANCE
# ============================================================

percentVar_all <- round(
  100 *
    pca_all$sdev^2 /
    sum(pca_all$sdev^2)
)

cat(
  "PC1:",
  percentVar_all[1],
  "%\n"
)

cat(
  "PC2:",
  percentVar_all[2],
  "%\n"
)


# ============================================================
# 9. CREATE PCA DATAFRAME
# ============================================================

pcaData_all <- data.frame(
  PC1 = pca_all$x[, 1],
  PC2 = pca_all$x[, 2],
  sample = rownames(pca_all$x),
  stringsAsFactors = FALSE
)

pcaData_all$cultivar <-
  metadata_pca[
    pcaData_all$sample,
    "cultivar"
  ]

pcaData_all$treatment <-
  metadata_pca[
    pcaData_all$sample,
    "treatment"
  ]


# ============================================================
# 10. CREATE PCA PLOT
# ============================================================

p_pca_all <- ggplot(
  pcaData_all,
  aes(
    x = PC1,
    y = PC2,
    color = cultivar,
    shape = treatment
  )
) +
  
  geom_point(
    size = 4
  ) +
  
  geom_text(
    aes(
      label = sample
    ),
    vjust = -1.7,
    size = 3.5,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "Catuai" = "#2563EB",
      "CR95" = "#8A2BE2"
    )
  ) +
  
  labs(
    title = "PCA of RNA-seq samples",
    
    x = paste0(
      "PC1: ",
      percentVar_all[1],
      "% variance"
    ),
    
    y = paste0(
      "PC2: ",
      percentVar_all[2],
      "% variance"
    ),
    
    color = "Cultivar",
    shape = "Treatment"
  ) +
  
  theme_bw(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )


# ============================================================
# 11. SHOW PCA
# ============================================================

print(
  p_pca_all
)


# ============================================================
# 12. SAVE PCA AS PDF
# ============================================================

ggsave(
  filename = file.path(
    outpathDE,
    "PCA_all_genes_Catuai_CR95_labeled_purple.pdf"
  ),
  plot = p_pca_all,
  width = 8,
  height = 7
)


# ============================================================
# 13. SAVE PCA COORDINATES
# ============================================================

write.table(
  pcaData_all,
  file = file.path(
    outpathDE,
    "PCA_all_genes_coordinates.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
### ANOVA 
# ============================================================
# 14. ANOVA ON PC SCORES  (cultivar + treatment)
# ============================================================

if (!requireNamespace("car", quietly = TRUE)) {
  install.packages("car")
}

library(car)


n_pc <- 4

pc_scores <- as.data.frame(
  pca_all$x[, 1:n_pc, drop = FALSE]
)

pc_scores$cultivar  <- metadata_pca[rownames(pc_scores), "cultivar"]
pc_scores$treatment <- metadata_pca[rownames(pc_scores), "treatment"]


pc_anova <- data.frame()

for (i in 1:n_pc) {
  
  pc <- paste0("PC", i)
  
  fit <- lm(
    as.formula(paste(pc, "~ cultivar + treatment")),
    data = pc_scores
  )
  
  a <- Anova(fit, type = "II")
  
  pc_anova <- rbind(
    pc_anova,
    data.frame(
      PC           = pc,
      PercentVar   = percentVar_all[i],
      F_cultivar   = round(a["cultivar",  "F value"], 2),
      p_cultivar   = signif(a["cultivar",  "Pr(>F)"], 3),
      F_treatment  = round(a["treatment", "F value"], 2),
      p_treatment  = signif(a["treatment", "Pr(>F)"], 3),
      df_resid     = a["Residuals", "Df"],
      R2           = round(summary(fit)$r.squared, 3)
    )
  )
}


print(pc_anova)


write.table(
  pc_anova,
  file = file.path(outpathDE, "PCA_score_ANOVA.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 11. LOAD ANNOTATION FOR KEGG
# ============================================================

annot_df <- read.delim(
  "fullAnnotation_clean copy.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Revisar nombres

colnames(
  annot_df
)


# ============================================================
# 12. AUTOMATICALLY IDENTIFY GENE ID COLUMN
# ============================================================

genes_deseq <-
  deseq_results[["CR95"]]$Geneid


overlap <- sapply(
  annot_df,
  
  function(x) {
    
    sum(
      as.character(x) %in%
        genes_deseq,
      na.rm = TRUE
    )
    
  }
)


print(
  sort(
    overlap,
    decreasing = TRUE
  )[1:min(10, length(overlap))]
)


gene_column <-
  names(
    which.max(
      overlap
    )
  )


cat(
  "\nGene ID column in annotation:",
  gene_column,
  "\n"
)


# ============================================================
# 13. JOIN DESeq2 + ANNOTATION
#
# ONLY SIGNIFICANT GENES
# padj < 0.1
# ============================================================


# ------------------------------------------------------------
# CR95
# ------------------------------------------------------------

kegg_CR95 <-
  deseq_results[["CR95"]] %>%
  
  filter(
    !is.na(padj),
    padj < 0.1
  ) %>%
  
  left_join(
    annot_df,
    by = setNames(
      gene_column,
      "Geneid"
    )
  )


# ------------------------------------------------------------
# CATUAI
# ------------------------------------------------------------

kegg_Catuai <-
  deseq_results[["Catuai"]] %>%
  
  filter(
    !is.na(padj),
    padj < 0.1
  ) %>%
  
  left_join(
    annot_df,
    by = setNames(
      gene_column,
      "Geneid"
    )
  )


# ============================================================
# 14. NUMBER OF SIGNIFICANT GENES WITH KEGG KO
# ============================================================

cat(
  "\nCR95 genes with KEGG KO:\n"
)


print(
  sum(
    !is.na(kegg_CR95$KEGG_ko) &
      kegg_CR95$KEGG_ko != "-" &
      kegg_CR95$KEGG_ko != ""
  )
)


cat(
  "\nCatuai genes with KEGG KO:\n"
)


print(
  sum(
    !is.na(kegg_Catuai$KEGG_ko) &
      kegg_Catuai$KEGG_ko != "-" &
      kegg_Catuai$KEGG_ko != ""
  )
)


# ============================================================
# 15. PREPARE CR95 FOR KEGG COLORING
# ============================================================

kegg_color_CR95 <-
  kegg_CR95 %>%
  
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
    
    KEGG_ko =
      trimws(
        KEGG_ko
      ),
    
    KO =
      sub(
        "^ko:",
        "",
        KEGG_ko
      ),
    
    Color =
      ifelse(
        log2FoldChange > 0,
        "red",
        "blue"
      )
  )


# ============================================================
# 16. PREPARE CATUAI FOR KEGG COLORING
# ============================================================

kegg_color_Catuai <-
  kegg_Catuai %>%
  
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
    
    KEGG_ko =
      trimws(
        KEGG_ko
      ),
    
    KO =
      sub(
        "^ko:",
        "",
        KEGG_ko
      ),
    
    Color =
      ifelse(
        log2FoldChange > 0,
        "red",
        "blue"
      )
  )


# ============================================================
# 17. SAVE KEGG COLOR FILES
# ============================================================


# ------------------------------------------------------------
# CR95
# ------------------------------------------------------------

KEGG_CR95_file <-
  kegg_color_CR95 %>%
  dplyr::select(
    KO,
    Color
  )


write.table(
  KEGG_CR95_file,
  file.path(
    outpathDE,
    "KEGG_Color_CR95_0.1.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


# ------------------------------------------------------------
# CATUAI
# ------------------------------------------------------------

KEGG_Catuai_file <-
  kegg_color_Catuai %>%
  dplyr::select(
    KO,
    Color
  )


write.table(
  KEGG_Catuai_file,
  file.path(
    outpathDE,
    "KEGG_Color_Catuai_0.1.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


# ============================================================
# 18. CLEAN KEGG FILES
# Remove duplicated KO + Color combinations
# ============================================================

KEGG_Catuai_clean <-
  kegg_color_Catuai %>%
  dplyr::select(
    KO,
    Color
  ) %>%
  dplyr::distinct()


KEGG_CR95_clean <-
  kegg_color_CR95 %>%
  dplyr::select(
    KO,
    Color
  ) %>%
  dplyr::distinct()


write.table(
  KEGG_Catuai_clean,
  file.path(
    outpathDE,
    "KEGG_Color_Catuai_clean_0.1.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


write.table(
  KEGG_CR95_clean,
  file.path(
    outpathDE,
    "KEGG_Color_CR95_clean_0.1.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
# ============================================================
# CR95 - KEGG ROJOS + GENE ID
# ============================================================

CR95_KEGG_RED_genes <-
  kegg_color_CR95 %>%
  
  dplyr::filter(
    Color == "red"
  ) %>%
  
  dplyr::select(
    Geneid,
    KO,
    Color,
    log2FoldChange,
    padj,
    Description,
    KEGG_Pathway
  ) %>%
  
  dplyr::distinct()


# Ver tabla
print(
  CR95_KEGG_RED_genes
)
write.table(
  CR95_KEGG_RED_genes,
  file = file.path(
    outpathDE,
    "CR95_KEGG_RED_with_GeneID_0.1.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
# ============================================================
# CATUAI - KEGG AZULES + GENE ID
# ============================================================

Catuai_KEGG_BLUE_genes <-
  kegg_color_Catuai %>%
  
  dplyr::filter(
    Color == "blue"
  ) %>%
  
  dplyr::select(
    Geneid,
    KO,
    Color,
    log2FoldChange,
    padj,
    Description,
    KEGG_Pathway
  ) %>%
  
  dplyr::distinct()


# Ver tabla
print(
  Catuai_KEGG_BLUE_genes
)
write.table(
  Catuai_KEGG_BLUE_genes,
  file = file.path(
    outpathDE,
    "Catuai_KEGG_BLUE_with_GeneID_0.1.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# unique KOs per cultivar
length(unique(kegg_color_Catuai$KO))
length(unique(kegg_color_CR95$KO))

# DEGs per KEGG pathway
pathway_counts <- function(df) {
  df %>%
    dplyr::filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-") %>%
    tidyr::separate_rows(KEGG_Pathway, sep = ",") %>%
    dplyr::mutate(KEGG_Pathway = trimws(KEGG_Pathway)) %>%
    dplyr::filter(grepl("^ko", KEGG_Pathway)) %>%
    dplyr::distinct(Geneid, KEGG_Pathway, Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
    dplyr::count(KEGG_Pathway, Direction) %>%
    tidyr::pivot_wider(names_from = Direction, values_from = n, values_fill = 0) %>%
    dplyr::mutate(Total = Up + Down) %>%
    dplyr::arrange(desc(Total))
}

head(pathway_counts(kegg_color_Catuai), 20)
head(pathway_counts(kegg_color_CR95), 20)





# ============================================================
# GENES CR95 CON MAYOR log2FoldChange
# padj < 0.1
# ============================================================

top_CR95_FC <- deseq_results[["CR95"]] %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.1,
    log2FoldChange > 0
  ) %>%
  dplyr::arrange(
    dplyr::desc(log2FoldChange)
  ) %>%
  dplyr::select(
    Geneid,
    log2FoldChange,
    padj,
    baseMean,
    Direction
  )

# Ver los primeros 20
head(top_CR95_FC, 150)
top_CR95_FC_annot <- kegg_CR95 %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.1,
    log2FoldChange > 0
  ) %>%
  dplyr::arrange(
    dplyr::desc(log2FoldChange)
  ) %>%
  dplyr::select(
    Geneid,
    log2FoldChange,
    padj,
    Direction,
    Description,
    KEGG_ko,
    KEGG_Pathway
  )

head(top_CR95_FC_annot, 150)

# ============================================================
# CREAR LISTAS UP / DOWN DIRECTAMENTE DESDE DESeq2 ACTUAL
# padj < 0.1
# ============================================================

CR95_UP <- deseq_results[["CR95"]] %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.1,
    log2FoldChange > 0
  ) %>%
  dplyr::pull(Geneid) %>%
  unique()

CR95_DOWN <- deseq_results[["CR95"]] %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.1,
    log2FoldChange < 0
  ) %>%
  dplyr::pull(Geneid) %>%
  unique()


Catuai_UP <- deseq_results[["Catuai"]] %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.1,
    log2FoldChange > 0
  ) %>%
  dplyr::pull(Geneid) %>%
  unique()

Catuai_DOWN <- deseq_results[["Catuai"]] %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.1,
    log2FoldChange < 0
  ) %>%
  dplyr::pull(Geneid) %>%
  unique()

cat("CR95 UP:", length(CR95_UP), "\n")
cat("CR95 DOWN:", length(CR95_DOWN), "\n")
cat("CR95 TOTAL:", length(CR95_UP) + length(CR95_DOWN), "\n")

cat("Catuai UP:", length(Catuai_UP), "\n")
cat("Catuai DOWN:", length(Catuai_DOWN), "\n")
cat("Catuai TOTAL:", length(Catuai_UP) + length(Catuai_DOWN), "\n")





# ============================================================
# EXPORTAR LISTAS DE GENES PARA ANALISIS DE PROMOTORES
# padj < 0.1
# ============================================================

writeLines(
  CR95_UP,
  "~/Desktop/CR95_UP_genes_0.1.txt"
)

writeLines(
  CR95_DOWN,
  "~/Desktop/CR95_DOWN_genes_0.1.txt"
)

writeLines(
  Catuai_UP,
  "~/Desktop/Catuai_UP_genes_0.1.txt"
)

writeLines(
  Catuai_DOWN,
  "~/Desktop/Catuai_DOWN_genes_0.1.txt"
)


# ============================================================
# VERIFICAR NUMERO DE GENES
# ============================================================

cat("CR95 UP:", length(CR95_UP), "\n")
cat("CR95 DOWN:", length(CR95_DOWN), "\n")
cat("Catuai UP:", length(Catuai_UP), "\n")
cat("Catuai DOWN:", length(Catuai_DOWN), "\n")


# Ver algunos Gene IDs
head(CR95_UP)
head(Catuai_DOWN)

# ============================================================
# GENES EXCLUSIVOS PARA ANALISIS DE MOTIVOS
# ============================================================

CR95_UP_only <- setdiff(
  CR95_UP,
  Catuai_DOWN
)

Catuai_DOWN_only <- setdiff(
  Catuai_DOWN,
  CR95_UP
)

# Los 17 con comportamiento opuesto
Shared_opposite <- intersect(
  CR95_UP,
  Catuai_DOWN
)

cat("CR95 UP exclusivos:", length(CR95_UP_only), "\n")
cat("Catuai DOWN exclusivos:", length(Catuai_DOWN_only), "\n")
cat("Compartidos CR95 UP / Catuai DOWN:", length(Shared_opposite), "\n")
writeLines(
  CR95_UP_only,
  "~/Desktop/CR95_UP_only_0.1.txt"
)

writeLines(
  Catuai_DOWN_only,
  "~/Desktop/Catuai_DOWN_only_0.1.txt"
)

writeLines(
  Shared_opposite,
  "~/Desktop/CR95_UP_Catuai_DOWN_shared17_0.1.txt"
)
# ============================================================
# 19. DEGs MAPPING TO A GIVEN KEGG PATHWAY
# ============================================================

pathway_genes <- function(df, ko_id) {
  df %>%
    dplyr::filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-") %>%
    tidyr::separate_rows(KEGG_Pathway, sep = ",") %>%
    dplyr::mutate(KEGG_Pathway = trimws(KEGG_Pathway)) %>%
    dplyr::filter(KEGG_Pathway == ko_id) %>%
    dplyr::select(Geneid, KO, log2FoldChange, padj, Description) %>%
    dplyr::distinct() %>%
    dplyr::arrange(desc(log2FoldChange))
}


# unique gene count (rows inflate when a gene carries several KOs)
n_genes <- function(df, ko_id) {
  length(unique(pathway_genes(df, ko_id)$Geneid))
}


for (ko in c("ko04626", "ko04016", "ko04075", "ko00270", "ko00195", "ko03010")) {
  
  cat("\n==============================\n")
  cat(ko, "\n")
  cat("==============================\n")
  
  cat("\nCatuai  — unique genes:", n_genes(kegg_color_Catuai, ko), "\n")
  print(as.data.frame(pathway_genes(kegg_color_Catuai, ko)))
  
  cat("\nCR95    — unique genes:", n_genes(kegg_color_CR95, ko), "\n")
  print(as.data.frame(pathway_genes(kegg_color_CR95, ko)))
}
for (ko in c("ko04626", "ko04016", "ko04075", "ko00270", "ko00195", "ko03010")) {
  for (cv in c("Catuai", "CR95")) {
    
    df <- if (cv == "Catuai") kegg_color_Catuai else kegg_color_CR95
    
    write.table(
      pathway_genes(df, ko),
      file = file.path(outpathDE, paste0("DEGs_", ko, "_", cv, "_0.1.txt")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}
##KEGG MAPPING
if (!requireNamespace("pathview", quietly = TRUE)) BiocManager::install("pathview")
library(pathview)

ko_vector <- function(df) {
  v <- df %>%
    dplyr::distinct(Geneid, KO, log2FoldChange) %>%
    dplyr::group_by(KO) %>%
    dplyr::summarise(lfc = mean(log2FoldChange), .groups = "drop")
  setNames(v$lfc, v$KO)
}

cat_v <- ko_vector(kegg_color_Catuai)
cr_v  <- ko_vector(kegg_color_CR95)

all_ko <- union(names(cat_v), names(cr_v))
mat <- cbind(Catuai = cat_v[all_ko], CR95 = cr_v[all_ko])
rownames(mat) <- all_ko          # keep NA for genes not DE in that cultivar

setwd(outpathDE)                 # pathview writes to the working directory

for (p in c("04626", "04016", "04075", "00270", "00195", "03010")) {
  pathview(
    gene.data    = mat,
    pathway.id   = p,
    species      = "ko",
    gene.idtype  = "KEGG",
    multi.state  = TRUE,
    same.layer   = FALSE,
    limit        = list(gene = 2.5),
    bins         = list(gene = 10),
    low  = list(gene = "blue"),
    mid  = list(gene = "gray95"),
    high = list(gene = "red"),
    na.col       = "transparent",
    out.suffix   = "Catuai_CR95",
    kegg.native  = TRUE
  )
}
#Heatmap para KEGG para comparar, que termino siendo como una verificación de los KEGG Color Pathways que no estaban del todo correctos
# ============================================================
# 20. COMPARATIVE HEATMAPS (Catuai vs CR95) PER KEGG PATHWAY
# ============================================================

library(pheatmap)


# ------------------------------------------------------------
# 20.1 PATHWAY TITLES
# ------------------------------------------------------------

pathway_titles <- c(
  ko04626 = "Plant-pathogen interaction",
  ko04016 = "MAPK signaling pathway - plant",
  ko04075 = "Plant hormone signal transduction",
  ko00270 = "Cysteine and methionine metabolism",
  ko00195 = "Photosynthesis",
  ko03010 = "Ribosome"
)


# ------------------------------------------------------------
# 20.2 SHORT GENE NAMES BY KO
# Genes without a match keep the eggNOG description.
# ------------------------------------------------------------

short_names <- c(
  K03781 = "Catalase",
  K13447 = "Rboh (NADPH oxidase)",
  K01183 = "Chitinase class I",
  K20547 = "Chitinase class I",
  K04730 = "LRR receptor-like kinase",
  K04733 = "LRR receptor-like kinase",
  K13420 = "LRR receptor-like kinase",
  K13416 = "SERK/BAK1",
  K02183 = "Calmodulin / CML",
  K13448 = "Calmodulin-like (CML)",
  K16465 = "Calmodulin-like (CML)",
  K10840 = "Calmodulin-like (CML)",
  K13412 = "Ca-dependent protein kinase",
  K05391 = "Cyclic nucleotide-gated channel",
  K13423 = "WRKY transcription factor",
  K13424 = "WRKY transcription factor",
  K13425 = "WRKY transcription factor",
  K18835 = "WRKY transcription factor",
  K13457 = "NB-LRR resistance protein",
  K13459 = "NB-ARC resistance protein",
  K14509 = "Ethylene receptor (ETR/ERS)",
  K14484 = "Aux/IAA",
  K14487 = "GH3 auxin-responsive",
  K14506 = "GH3 / JAR1",
  K14488 = "SAUR auxin-responsive",
  K14489 = "Cytokinin receptor (CRE1)",
  K14508 = "NPR1/NPR3-like",
  K14504 = "Xyloglucan transferase (TCH4)",
  K08235 = "Xyloglucan transferase (TCH4)",
  K05933 = "ACC oxidase",
  K00815 = "Tyrosine aminotransferase",
  K01738 = "Cysteine synthase",
  K00026 = "Malate dehydrogenase",
  K17686 = "Copper-transporting ATPase",
  K12619 = "5'-3' exoribonuclease",
  K20553 = "5'-3' exoribonuclease",
  K02639 = "Ferredoxin"
)


# ------------------------------------------------------------
# 20.3 GENES EXCLUDED FROM SPECIFIC PATHWAYS
# KEGG ortholog assignment inconsistent with the NCBI
# annotation of the Coffea arabica gene.
# ------------------------------------------------------------

exclude_genes <- list(
  ko00195 = c("LOC113723029", "LOC113698731")   # ABC transporter B15
)

# ------------------------------------------------------------
# 20.4 HEATMAP FUNCTION
# ------------------------------------------------------------

comp_heatmap <- function(ko_id, file, use_short = TRUE, desc_chars = 70) {
  
  
  g <- dplyr::bind_rows(
    pathway_genes(kegg_color_Catuai, ko_id) %>%
      dplyr::mutate(cultivar = "Catuai"),
    
    pathway_genes(kegg_color_CR95, ko_id) %>%
      dplyr::mutate(cultivar = "CR95")
  )
  
  
  if (!is.null(exclude_genes[[ko_id]])) {
    g <- g %>%
      dplyr::filter(!Geneid %in% exclude_genes[[ko_id]])
  }
  
  
  if (nrow(g) == 0) {
    message("No DEGs mapping to ", ko_id)
    return(invisible(NULL))
  }
  
  
  lab <- g %>%
    dplyr::mutate(short = unname(short_names[KO])) %>%
    dplyr::group_by(Geneid) %>%
    dplyr::summarise(
      short       = dplyr::first(stats::na.omit(short)),
      Description = dplyr::first(Description),
      .groups     = "drop"
    ) %>%
    dplyr::mutate(
      name = dplyr::case_when(
        use_short & !is.na(short) ~ short,
        TRUE ~ substr(Description, 1, desc_chars)
      ),
      label = paste0(Geneid, " - ", name)
    )
  
  
  m <- g %>%
    dplyr::distinct(Geneid, cultivar, log2FoldChange) %>%
    tidyr::pivot_wider(
      names_from  = cultivar,
      values_from = log2FoldChange
    ) %>%
    as.data.frame()
  
  rownames(m) <- lab$label[match(m$Geneid, lab$Geneid)]
  
  m <- as.matrix(m[, setdiff(colnames(m), "Geneid"), drop = FALSE])
  
  
  for (cv in c("Catuai", "CR95")) {
    if (!cv %in% colnames(m)) {
      m <- cbind(m, NA_real_)
      colnames(m)[ncol(m)] <- cv
    }
  }
  
  m <- m[, c("Catuai", "CR95"), drop = FALSE]
  
  
  ord <- order(
    !complete.cases(m),
    -m[, "Catuai"],
    -m[, "CR95"],
    na.last = TRUE
  )
  
  m <- m[ord, , drop = FALSE]
  
  
  lim <- max(abs(m), na.rm = TRUE)
  
  
  pheatmap(
    m,
    color        = colorRampPalette(c("blue", "white", "red"))(51),
    breaks       = seq(-lim, lim, length.out = 52),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    na_col       = "grey90",
    border_color = "grey70",
    fontsize_row = 8,
    fontsize_col = 11,
    angle_col    = 0,
    labels_col   = c("Catuaí", "CR95"),
    main         = ifelse(
      is.na(pathway_titles[ko_id]),
      ko_id,
      paste0(pathway_titles[ko_id], " (", ko_id, ")")
    ),
    filename     = file,
    width        = 11,
    height       = max(6, 0.16 * nrow(m) + 2)
  )
  
  invisible(m)
}


# ============================================================
# 21. GENERATE HEATMAPS
# ============================================================

for (ko in c("ko04626", "ko04016", "ko04075", "ko00270", "ko00195", "ko03010")) {
  
  comp_heatmap(
    ko,
    file.path(
      outpathDE,
      paste0("heatmap_", ko, "_Catuai_vs_CR95.pdf")
    )
  )
}


list.files(outpathDE, pattern = "^heatmap_")

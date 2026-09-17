# ============================================================
# CAMERA GO ENRICHMENT
# HISAT2 counts
# Catuai y CR95
# Xylella vs saline
#
# SIGNIFICANCE CRITERION:
# FDR < 0.1
#
# CAMERA:
# Up   = GO term enriched toward genes higher in Xylella
# Down = GO term enriched toward genes lower in Xylella
#
# VISUAL STYLE:
# Lilac / lavender palette
# ============================================================


# ============================================================
# 1. WORKING DIRECTORY
# ============================================================

setwd("~/Desktop/Changed padj 0.1")


# ============================================================
# 2. OUTPUT DIRECTORY
# ============================================================

outpathcount <- file.path(
  "~/Desktop/Changed padj 0.1",
  "CAMERA RESULTS"
)

dir.create(
  outpathcount,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 3. INSTALL / LOAD LIBRARIES
# ============================================================

options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)


# ------------------------------------------------------------
# CRAN packages
# ------------------------------------------------------------

cran_pkgs <- c(
  "ggplot2",
  "BiocManager"
)


for (pkg in cran_pkgs) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(pkg)
    
  }
}


# ------------------------------------------------------------
# Bioconductor packages
# ------------------------------------------------------------

bioc_pkgs <- c(
  "edgeR",
  "limma",
  "GO.db",
  "AnnotationDbi"
)


for (pkg in bioc_pkgs) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    BiocManager::install(pkg)
    
  }
}


# ------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------

library(edgeR)
library(limma)
library(ggplot2)
library(GO.db)
library(AnnotationDbi)


# ============================================================
# 4. LOAD HISAT2 COUNT MATRIX
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
# 6. FIX / ALIGN SAMPLE NAMES
# ============================================================

# Matrix_hisat2 puede tener Cax7
# mientras que metadata tiene Ca7x

colnames(counts)[
  colnames(counts) == "Cax7"
] <- "Ca7x"


# ------------------------------------------------------------
# Confirmar que todas las muestras existen
# ------------------------------------------------------------

stopifnot(
  all(
    metadata$sample %in%
      colnames(counts)
  )
)


# ------------------------------------------------------------
# Ordenar counts igual que metadata
# ------------------------------------------------------------

counts <- counts[
  ,
  metadata$sample,
  drop = FALSE
]


# ------------------------------------------------------------
# Confirmar mismo orden
# ------------------------------------------------------------

stopifnot(
  identical(
    colnames(counts),
    metadata$sample
  )
)


cat(
  "Samples loaded:",
  ncol(counts),
  "\n"
)


# ============================================================
# 7. LOAD FUNCTIONAL ANNOTATION
# ============================================================

annotation <- read.table(
  "fullAnnotation_clean copy.txt",
  header = TRUE,
  sep = "\t",
  quote = "\"",
  comment.char = "",
  fill = TRUE,
  stringsAsFactors = FALSE,
  na.strings = c(
    "-",
    "NA",
    ""
  )
)


# ============================================================
# 8. CREATE GENE -> GO RELATIONSHIP
# ============================================================

has_go <-
  !is.na(annotation$GOs) &
  !is.na(annotation$gene_id)


gene2go <- strsplit(
  annotation$GOs[has_go],
  ",\\s*"
)


names(gene2go) <-
  annotation$gene_id[has_go]


cat(
  "Genes with GO annotation:",
  length(gene2go),
  "\n"
)


# ============================================================
# 9. CONVERT GENE -> GO INTO GO -> GENES
# ============================================================

term2gene <- split(
  
  rep(
    names(gene2go),
    lengths(gene2go)
  ),
  
  unlist(gene2go)
)


# ------------------------------------------------------------
# Remove duplicated LOC IDs within each GO term
# ------------------------------------------------------------

term2gene <- lapply(
  term2gene,
  unique
)


cat(
  "Total GO terms:",
  length(term2gene),
  "\n"
)


# ============================================================
# 10. MINIMUM GO SET SIZE
# ============================================================

min_set_size <- 5


# ============================================================
# 11. IMPORTANT CHECKS
# ============================================================

cat(
  "\nGenes shared between counts and annotation:",
  
  sum(
    rownames(counts) %in%
      annotation$gene_id
  ),
  
  "of",
  nrow(counts),
  "\n"
)


cat(
  "Total GO terms:",
  length(term2gene),
  "\n"
)


cat(
  "GO terms with >= 5 genes:",
  
  sum(
    lengths(term2gene) >=
      min_set_size
  ),
  
  "\n"
)


cat(
  "\nSamples by cultivar and treatment:\n"
)


print(
  table(
    metadata$cultivar,
    metadata$treatment
  )
)


# ============================================================
# 12. CAMERA BY CULTIVAR
# ============================================================

cultivars <- c(
  "Catuai",
  "CR95"
)


# ------------------------------------------------------------
# Lists to save results
# ------------------------------------------------------------

camera_results <- list()

sig_results <- list()

voom_objects <- list()

genes_analyzed <- list()


# ------------------------------------------------------------
# Summary table
# ------------------------------------------------------------

CAMERA_summary <- data.frame()


# ============================================================
# 13. LOOP CAMERA
# ============================================================

for (cv in cultivars) {
  
  
  cat(
    "\n====================================\n"
  )
  
  
  cat(
    "Running CAMERA for:",
    cv,
    "\n"
  )
  
  
  cat(
    "====================================\n"
  )
  
  
  # ==========================================================
  # 13.1 SELECT CULTIVAR
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
  
  
  stopifnot(
    ncol(counts_cv) ==
      nrow(metadata_cv)
  )
  
  
  stopifnot(
    identical(
      colnames(counts_cv),
      metadata_cv$sample
    )
  )
  
  
  # ==========================================================
  # 13.2 TREATMENT
  #
  # saline = reference
  # xylella = condition
  #
  # UP:
  # enrichment toward genes higher in Xylella
  #
  # DOWN:
  # enrichment toward genes lower in Xylella
  # ==========================================================
  
  treatment <- factor(
    
    metadata_cv$treatment,
    
    levels = c(
      "saline",
      "xylella"
    )
  )
  
  
  cat(
    "\nTreatment samples:\n"
  )
  
  
  print(
    table(treatment)
  )
  
  
  # ==========================================================
  # 13.3 EDGER
  # ==========================================================
  
  dge <- DGEList(
    counts = counts_cv,
    group = treatment
  )
  
  
  # ----------------------------------------------------------
  # Filter low expressed genes
  # ----------------------------------------------------------
  
  keep <- filterByExpr(
    dge,
    group = treatment
  )
  
  
  dge <- dge[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # TMM normalization
  # ----------------------------------------------------------
  
  dge <- calcNormFactors(
    dge
  )
  
  
  cat(
    "Genes after edgeR filtering:",
    nrow(dge),
    "\n"
  )
  
  
  # ==========================================================
  # 13.4 DESIGN MATRIX
  # ==========================================================
  
  design <- model.matrix(
    ~ treatment
  )
  
  
  cat(
    "\nDesign matrix columns:\n"
  )
  
  
  print(
    colnames(design)
  )
  
  
  # ==========================================================
  # 13.5 VOOM
  # ==========================================================
  
  v <- voom(
    dge,
    design,
    plot = FALSE
  )
  
  
  voom_objects[[cv]] <- v
  
  
  # ----------------------------------------------------------
  # Number of genes actually used in CAMERA
  # ----------------------------------------------------------
  
  genes_analyzed[[cv]] <- nrow(v$E)
  
  
  cat(
    "Genes analyzed in CAMERA:",
    genes_analyzed[[cv]],
    "\n"
  )
  
  
  # ==========================================================
  # 13.6 PREPARE GO SETS
  # ==========================================================
  
  term2gene_cv <- lapply(
    
    term2gene,
    
    function(g) {
      
      intersect(
        g,
        rownames(v$E)
      )
      
    }
  )
  
  
  # ----------------------------------------------------------
  # Keep GO terms with >= 5 genes after filtering
  # ----------------------------------------------------------
  
  term2gene_cv <-
    term2gene_cv[
      
      lengths(term2gene_cv) >=
        min_set_size
      
    ]
  
  
  cat(
    "GO terms available for CAMERA:",
    length(term2gene_cv),
    "\n"
  )
  
  
  # ==========================================================
  # 13.7 CONVERT GENES INTO MATRIX POSITIONS
  # ==========================================================
  
  idx <- lapply(
    
    term2gene_cv,
    
    function(g) {
      
      match(
        g,
        rownames(v$E)
      )
      
    }
  )
  
  
  # ==========================================================
  # 13.8 RUN CAMERA
  # ==========================================================
  
  res <- camera(
    
    v,
    
    index = idx,
    
    design = design,
    
    contrast = 2
  )
  
  
  # ==========================================================
  # 13.9 ADD GO ID
  # ==========================================================
  
  res$GOID <-
    rownames(res)
  
  
  rownames(res) <- NULL
  
  
  # ==========================================================
  # 13.10 ADD GO TERM NAME + ONTOLOGY
  # ==========================================================
  
  term_info <- suppressMessages(
    
    AnnotationDbi::select(
      
      GO.db,
      
      keys = unique(
        res$GOID
      ),
      
      columns = c(
        "TERM",
        "ONTOLOGY"
      ),
      
      keytype = "GOID"
    )
  )
  
  
  res <- merge(
    
    res,
    
    term_info,
    
    by = "GOID",
    
    all.x = TRUE
  )
  
  
  # ----------------------------------------------------------
  # If GO.db cannot find term name, show GO ID
  # ----------------------------------------------------------
  
  res$TERM[
    is.na(res$TERM)
  ] <-
    res$GOID[
      is.na(res$TERM)
    ]
  
  
  res$ONTOLOGY[
    is.na(res$ONTOLOGY)
  ] <- "Unknown"
  
  
  # ==========================================================
  # 13.11 ADD SIGNIFICANCE COLUMN
  #
  # FDR < 0.1
  # ==========================================================
  
  res$Significance <-
    "Not significant"
  
  
  res$Significance[
    
    !is.na(res$FDR) &
      res$FDR < 0.1
    
  ] <- "Significant"
  
  
  # ==========================================================
  # 13.12 SORT BY FDR
  # ==========================================================
  
  res <- res[
    
    order(
      res$FDR,
      na.last = TRUE
    ),
    
    ,
    
    drop = FALSE
  ]
  
  
  # ==========================================================
  # 13.13 SAVE COMPLETE CAMERA RESULTS
  # ==========================================================
  
  camera_results[[cv]] <-
    res
  
  
  write.table(
    
    res,
    
    file = file.path(
      
      outpathcount,
      
      paste0(
        "CAMERA_",
        cv,
        "_xylella_vs_saline_ALL_0.1.txt"
      )
    ),
    
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
  )
  
  
  # ==========================================================
  # 13.14 SIGNIFICANT GO TERMS
  #
  # FDR < 0.1
  # ==========================================================
  
  sig <- res[
    
    !is.na(res$FDR) &
      res$FDR < 0.1,
    
    ,
    
    drop = FALSE
  ]
  
  
  sig_results[[cv]] <-
    sig
  
  
  # ==========================================================
  # 13.15 SAVE SIGNIFICANT GO TERMS
  # ==========================================================
  
  write.table(
    
    sig,
    
    file = file.path(
      
      outpathcount,
      
      paste0(
        "CAMERA_",
        cv,
        "_xylella_vs_saline_SIGNIFICANT_0.1.txt"
      )
    ),
    
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
  )
  
  
  # ==========================================================
  # 13.16 PRINT SUMMARY
  # ==========================================================
  
  cat(
    "\nCAMERA summary for",
    cv,
    "\n"
  )
  
  
  cat(
    "Genes analyzed:",
    genes_analyzed[[cv]],
    "\n"
  )
  
  
  cat(
    "GO terms analyzed:",
    nrow(res),
    "\n"
  )
  
  
  cat(
    "Significant GO terms:",
    nrow(sig),
    "\n"
  )
  
  
  cat(
    "\nSignificant directions:\n"
  )
  
  
  print(
    table(
      sig$Direction
    )
  )
  
  
  # ==========================================================
  # 13.17 ADD TO SUMMARY TABLE
  # ==========================================================
  
  CAMERA_summary <- rbind(
    
    CAMERA_summary,
    
    data.frame(
      
      Cultivar = cv,
      
      Genes_analyzed =
        genes_analyzed[[cv]],
      
      GO_terms_analyzed =
        nrow(res),
      
      Significant_GO =
        nrow(sig),
      
      Up =
        sum(
          sig$Direction == "Up",
          na.rm = TRUE
        ),
      
      Down =
        sum(
          sig$Direction == "Down",
          na.rm = TRUE
        )
    )
  )
}


# ============================================================
# 14. FINAL CAMERA SUMMARY
# ============================================================

cat(
  "\n====================================\n"
)

cat(
  "FINAL CAMERA SUMMARY\n"
)

cat(
  "FDR < 0.1\n"
)

cat(
  "====================================\n"
)


print(
  CAMERA_summary
)


# ============================================================
# 15. SAVE CAMERA SUMMARY
# ============================================================

write.table(
  
  CAMERA_summary,
  
  file = file.path(
    outpathcount,
    "CAMERA_summary_0.1.txt"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 16. NUMBER OF SIGNIFICANT GO TERMS
# ============================================================

cat(
  "\nNumber of significant GO terms:\n"
)


print(
  sapply(
    sig_results,
    nrow
  )
)


# ============================================================
# 17. CHECK UP / DOWN
# ============================================================

cat(
  "\nCatuai significant directions:\n"
)


print(
  table(
    sig_results[["Catuai"]]$Direction
  )
)


cat(
  "\nCR95 significant directions:\n"
)


print(
  table(
    sig_results[["CR95"]]$Direction
  )
)


# ============================================================
# 18. CAMERA BARPLOT FUNCTION
#
# MODERN LILAC DESIGN
# Top 15 significant GO terms
# ============================================================

plot_camera_barplot <- function(
    res,
    cv,
    n_genes_analyzed,
    top_n = 15
) {
  
  
  # ----------------------------------------------------------
  # If no significant GO terms
  # ----------------------------------------------------------
  
  if (nrow(res) == 0) {
    
    cat(
      "No significant GO terms for",
      cv,
      "-- skipping barplot.\n"
    )
    
    return(
      invisible(NULL)
    )
  }
  
  
  # ----------------------------------------------------------
  # Select Top GO terms by lowest FDR
  # ----------------------------------------------------------
  
  res_top <- head(
    
    res[
      order(res$FDR),
      ,
      drop = FALSE
    ],
    
    top_n
  )
  
  
  # ----------------------------------------------------------
  # Order GO terms in plot
  # ----------------------------------------------------------
  
  res_top$TERM <- factor(
    
    res_top$TERM,
    
    levels = rev(
      res_top$TERM
    )
  )
  
  
  # ----------------------------------------------------------
  # Create plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    
    res_top,
    
    aes(
      x = TERM,
      y = -log10(FDR),
      fill = Direction
    )
    
  ) +
    
    # --------------------------------------------------------
  # Bars
  # --------------------------------------------------------
  
  geom_col(
    width = 0.70,
    alpha = 0.90
  ) +
    
    # --------------------------------------------------------
  # Dot at end of bar
  # --------------------------------------------------------
  
  geom_point(
    aes(
      color = Direction
    ),
    size = 4
  ) +
    
    # --------------------------------------------------------
  # FDR = 0.1 threshold
  # --------------------------------------------------------
  
  geom_hline(
    yintercept = -log10(0.1),
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey45"
  ) +
    
    coord_flip() +
    
    # --------------------------------------------------------
  # LILAC / LAVENDER COLORS
  # --------------------------------------------------------
  
  scale_fill_manual(
    
    values = c(
      
      Up = "#C084FC",
      
      Down = "#818CF8"
      
    )
    
  ) +
    
    scale_color_manual(
      
      values = c(
        
        Up = "#9333EA",
        
        Down = "#6366F1"
        
      ),
      
      guide = "none"
    ) +
    
    scale_y_continuous(
      
      expand = expansion(
        mult = c(
          0,
          0.10
        )
      )
    ) +
    
    labs(
      
      title = paste0(
        cv,
        " | CAMERA GO enrichment"
      ),
      
      subtitle = paste0(
        "Xylella vs saline  •   Genes analyzed: ",
        format(
          n_genes_analyzed,
          big.mark = ","
        )
      ),
      
      x = NULL,
      
      y = expression(
        -log[10](FDR)
      ),
      
      fill =
        "Enrichment direction",
      
      caption =
        "Dashed line indicates the FDR = 0.1 significance threshold"
    ) +
    
    theme_minimal(
      base_size = 12
    ) +
    
    theme(
      
      # ------------------------------------------------------
      # Main title
      # ------------------------------------------------------
      
      plot.title = element_text(
        size = 18,
        face = "bold",
        color = "#6D28D9",
        margin = margin(
          b = 5
        )
      ),
      
      # ------------------------------------------------------
      # Subtitle
      # ------------------------------------------------------
      
      plot.subtitle = element_text(
        size = 11,
        color = "grey35",
        margin = margin(
          b = 17
        )
      ),
      
      # ------------------------------------------------------
      # GO term names
      # ------------------------------------------------------
      
      axis.text.y = element_text(
        size = 10,
        color = "grey15"
      ),
      
      axis.text.x = element_text(
        size = 9.5,
        color = "grey35"
      ),
      
      axis.title.x = element_text(
        size = 11,
        face = "bold",
        margin = margin(
          t = 10
        )
      ),
      
      # ------------------------------------------------------
      # Grid
      # ------------------------------------------------------
      
      panel.grid.major.y =
        element_blank(),
      
      panel.grid.minor =
        element_blank(),
      
      panel.grid.major.x =
        element_line(
          linewidth = 0.35,
          color = "grey88"
        ),
      
      # ------------------------------------------------------
      # Legend
      # ------------------------------------------------------
      
      legend.position =
        "top",
      
      legend.title = element_text(
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 10
      ),
      
      # ------------------------------------------------------
      # Caption
      # ------------------------------------------------------
      
      plot.caption = element_text(
        size = 9,
        color = "grey50",
        hjust = 0,
        margin = margin(
          t = 12
        )
      ),
      
      plot.margin = margin(
        15,
        25,
        15,
        15
      )
    )
  
  
  # ----------------------------------------------------------
  # Show plot
  # ----------------------------------------------------------
  
  print(p)
  
  
  # ----------------------------------------------------------
  # Save PDF
  # ----------------------------------------------------------
  
  ggsave(
    
    filename = file.path(
      
      outpathcount,
      
      paste0(
        "Barplot_CAMERA_",
        cv,
        "_xylella_vs_saline_LILAC_0.1.pdf"
      )
    ),
    
    plot = p,
    
    width = 10,
    height = 7
  )
}


# ============================================================
# 19. CREATE CATUAI BARPLOT
# ============================================================

plot_camera_barplot(
  sig_results[["Catuai"]],
  "Catuai",
  genes_analyzed[["Catuai"]],
  top_n = 15
)


# ============================================================
# 20. CREATE CR95 BARPLOT
# ============================================================

plot_camera_barplot(
  sig_results[["CR95"]],
  "CR95",
  genes_analyzed[["CR95"]],
  top_n = 15
)

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


# ------------------------------------------------------------
# Install / load UpSetR
# ------------------------------------------------------------

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}

library(UpSetR)


# ============================================================
# 1. CREATE UP / DOWN GENE LISTS
# ============================================================


# ------------------------------------------------------------
# CATUAI UP
# padj < 0.1
# log2FoldChange > 0
# ------------------------------------------------------------

Catuai_UP <-
  deseq_results[["Catuai"]]$Geneid[
    
    !is.na(
      deseq_results[["Catuai"]]$padj
    ) &
      
      deseq_results[["Catuai"]]$padj < 0.1 &
      
      deseq_results[["Catuai"]]$log2FoldChange > 0
  ]


# ------------------------------------------------------------
# CATUAI DOWN
# padj < 0.1
# log2FoldChange < 0
# ------------------------------------------------------------

Catuai_DOWN <-
  deseq_results[["Catuai"]]$Geneid[
    
    !is.na(
      deseq_results[["Catuai"]]$padj
    ) &
      
      deseq_results[["Catuai"]]$padj < 0.1 &
      
      deseq_results[["Catuai"]]$log2FoldChange < 0
  ]


# ------------------------------------------------------------
# CR95 UP
# ------------------------------------------------------------

CR95_UP <-
  deseq_results[["CR95"]]$Geneid[
    
    !is.na(
      deseq_results[["CR95"]]$padj
    ) &
      
      deseq_results[["CR95"]]$padj < 0.1 &
      
      deseq_results[["CR95"]]$log2FoldChange > 0
  ]


# ------------------------------------------------------------
# CR95 DOWN
# ------------------------------------------------------------

CR95_DOWN <-
  deseq_results[["CR95"]]$Geneid[
    
    !is.na(
      deseq_results[["CR95"]]$padj
    ) &
      
      deseq_results[["CR95"]]$padj < 0.1 &
      
      deseq_results[["CR95"]]$log2FoldChange < 0
  ]


# ============================================================
# 2. UPSET LIST
# ============================================================

strategy_list <- list(
  
  "Catuai UP" =
    Catuai_UP,
  
  "Catuai DOWN" =
    Catuai_DOWN,
  
  "CR95 UP" =
    CR95_UP,
  
  "CR95 DOWN" =
    CR95_DOWN
)


# ------------------------------------------------------------
# Check number of genes
# ------------------------------------------------------------

cat(
  "\nGenes in each UpSet group:\n"
)


print(
  sapply(
    strategy_list,
    length
  )
)


# ============================================================
# 3. SHARED SIGNIFICANT GENES
# ============================================================

Catuai_sig <- union(
  Catuai_UP,
  Catuai_DOWN
)

CR95_sig <- union(
  CR95_UP,
  CR95_DOWN
)


shared_genes <- intersect(
  Catuai_sig,
  CR95_sig
)


cat(
  "\nTotal significant genes shared between Catuai and CR95:",
  length(shared_genes),
  "\n"
)


# ============================================================
# 4. SHOW UPSET IN RSTUDIO
#
# Purple / turquoise design
# ============================================================

UpSetR::upset(
  
  UpSetR::fromList(
    strategy_list
  ),
  
  sets = c(
    "Catuai UP",
    "Catuai DOWN",
    "CR95 UP",
    "CR95 DOWN"
  ),
  
  keep.order = TRUE,
  
  # ----------------------------------------------------------
  # COLORS
  # ----------------------------------------------------------
  
  # Top intersection bars
  main.bar.color =
    "#7C3AED",
  
  # Left set-size bars
  sets.bar.color =
    "#14B8A6",
  
  # Intersection dots / connecting lines
  matrix.color =
    "#6D28D9",
  
  # Alternating matrix background
  shade.color =
    "#EDE9FE",
  
  shade.alpha =
    0.35,
  
  # ----------------------------------------------------------
  # SIZE
  # ----------------------------------------------------------
  
  point.size =
    3.5,
  
  line.size =
    1,
  
  text.scale = c(
    1.4,
    1.2,
    1.2,
    1.1,
    1.2,
    1
  ),
  
  nsets = 4,
  
  order.by =
    "freq",
  
  mainbar.y.label =
    "Number of shared genes",
  
  sets.x.label =
    "Number of significant genes"
)


# ============================================================
# 5. SAVE UPSET PDF
# ============================================================

pdf(
  
  file.path(
    outpathDE,
    "UpSet_Catuai_CR95_UP_DOWN_COLOR_0.1.pdf"
  ),
  
  width = 11,
  height = 7
)


UpSetR::upset(
  
  UpSetR::fromList(
    strategy_list
  ),
  
  sets = c(
    "Catuai UP",
    "Catuai DOWN",
    "CR95 UP",
    "CR95 DOWN"
  ),
  
  keep.order = TRUE,
  
  main.bar.color =
    "#7C3AED",
  
  sets.bar.color =
    "#14B8A6",
  
  matrix.color =
    "#6D28D9",
  
  shade.color =
    "#EDE9FE",
  
  shade.alpha =
    0.35,
  
  point.size =
    3.5,
  
  line.size =
    1,
  
  text.scale = c(
    1.4,
    1.2,
    1.2,
    1.1,
    1.2,
    1
  ),
  
  nsets = 4,
  
  order.by =
    "freq",
  
  mainbar.y.label =
    "Number of shared genes",
  
  sets.x.label =
    "Number of significant genes"
)


dev.off()


#====================================================
# 29. COMPARISON DOT PLOT
# ORIGINAL DESIGN
# ============================================================

p_comparison <- ggplot(
  
  comparison_df,
  
  aes(
    
    x = cultivar,
    
    y = TERM,
    
    size = NGenes,
    
    color = -log10(FDR),
    
    shape = Direction
  )
  
) +
  
  geom_point(
    alpha = 0.85
  ) +
  
  scale_color_viridis_c() +
  
  labs(
    
    title =
      "CAMERA GO terms compared across cultivars",
    
    subtitle =
      "Xylella vs saline | significant selection FDR < 0.1",
    
    x =
      "Cultivar",
    
    y =
      "GO term",
    
    size =
      "Genes in set",
    
    color =
      "-log10(FDR)",
    
    shape =
      "Direction"
  ) +
  
  theme_bw()


# ============================================================
# 30. SHOW DOT PLOT
# ============================================================

print(
  p_comparison
)


# ============================================================
# 31. SAVE DOT PLOT
# ============================================================

ggsave(
  
  filename = file.path(
    
    outpathcount,
    
    "Comparison_dotplot_Catuai_vs_CR95_0.1.pdf"
  ),
  
  plot = p_comparison,
  
  width = 10,
  
  height = 10
)

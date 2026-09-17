# ============================================================
# WGCNA - ANALISIS COMPLETO
# + DESeq2 CR95 Xylella vs saline
#
# DESeq2 significance:
# padj < 0.1
#
# TODOS LOS ARCHIVOS GENERADOS TERMINAN EN _0.1
# ============================================================


# ============================================================
# 1. WORKING DIRECTORY
# ============================================================

setwd(
  "~/Desktop/Changed padj 0.1"
)


outpathcount <- file.path(
  "~/Desktop/Changed padj 0.1",
  "CYTOSCAPE"
)


dir.create(
  outpathcount,
  showWarnings = FALSE,
  recursive = TRUE
)


# Carpeta donde están los resultados DESeq2

deseq_path <- file.path(
  "~/Desktop/Changed padj 0.1",
  "DESeq2_RESULTS"
)


options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)

options(
  stringsAsFactors = FALSE
)


# ============================================================
# 2. INSTALL PACKAGES
# ============================================================

cran_pkgs <- c(
  "WGCNA",
  "dynamicTreeCut",
  "fastcluster",
  "foreach",
  "doParallel",
  "BiocManager"
)


for (pkg in cran_pkgs) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(pkg)
    
  }
}


bioc_pkgs <- c(
  "impute",
  "preprocessCore",
  "GO.db",
  "AnnotationDbi",
  "DESeq2"
)


for (pkg in bioc_pkgs) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    BiocManager::install(pkg)
    
  }
}


# ============================================================
# 3. LOAD LIBRARIES
# ============================================================

library(WGCNA)
library(DESeq2)
library(AnnotationDbi)
library(GO.db)


# Evitar problemas de paralelización

disableWGCNAThreads()


# ============================================================
# 4. LOAD COUNT MATRIX
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
# 6. FIX SAMPLE NAME
# ============================================================

colnames(counts)[
  colnames(counts) == "Cax7"
] <- "Ca7x"


# ============================================================
# 7. CHECK SAMPLES
# ============================================================

missing_samples <- setdiff(
  metadata$sample,
  colnames(counts)
)


if (length(missing_samples) > 0) {
  
  stop(
    "These metadata samples do not exist in counts: ",
    paste(
      missing_samples,
      collapse = ", "
    )
  )
}


# Ordenar counts igual que metadata

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


cat(
  "Samples:",
  ncol(counts),
  "\n"
)


# ============================================================
# 8. FILTER LOW COUNT GENES
# ============================================================

counts <- round(
  counts
)


keep <- rowSums(
  counts
) >= 10


counts <- counts[
  keep,
  ,
  drop = FALSE
]


cat(
  "Genes after count filtering:",
  nrow(counts),
  "\n"
)


# ============================================================
# 9. VST TRANSFORMATION
# ============================================================

col_data <- data.frame(
  
  dummy = rep(
    1,
    nrow(metadata)
  ),
  
  row.names =
    metadata$sample
)


dds <- DESeqDataSetFromMatrix(
  
  countData =
    counts,
  
  colData =
    col_data,
  
  design =
    ~ 1
)


vsd <- vst(
  dds,
  blind = TRUE
)


vst_mat <- assay(
  vsd
)


cat(
  "VST matrix:",
  nrow(vst_mat),
  "genes x",
  ncol(vst_mat),
  "samples\n"
)


# ============================================================
# 10. SELECT TOP 5000 MOST VARIABLE GENES
# ============================================================

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
  
  sort(
    gene_var,
    decreasing = TRUE
  )
  
)[1:n_top_genes]


# WGCNA:
# rows    = samples
# columns = genes

datExpr <- as.data.frame(
  
  t(
    vst_mat[
      top_genes,
      ,
      drop = FALSE
    ]
  )
)


cat(
  "WGCNA expression matrix:",
  nrow(datExpr),
  "samples x",
  ncol(datExpr),
  "genes\n"
)


# ============================================================
# 11. DATA QUALITY
# ============================================================

gsg <- goodSamplesGenes(
  datExpr,
  verbose = 3
)


cat(
  "All genes/samples OK:",
  gsg$allOK,
  "\n"
)


if (!gsg$allOK) {
  
  datExpr <- datExpr[
    gsg$goodSamples,
    gsg$goodGenes,
    drop = FALSE
  ]
  
  
  cat(
    "Removed",
    sum(!gsg$goodGenes),
    "genes and",
    sum(!gsg$goodSamples),
    "samples.\n"
  )
}


cat(
  "NA values:",
  sum(
    is.na(datExpr)
  ),
  "\n"
)


cat(
  "Infinite values:",
  sum(
    is.infinite(
      as.matrix(datExpr)
    )
  ),
  "\n"
)


cat(
  "Final samples:",
  nrow(datExpr),
  "\n"
)


cat(
  "Final genes:",
  ncol(datExpr),
  "\n"
)


# ============================================================
# 12. SAMPLE CLUSTERING
# ============================================================

sampleTree <- hclust(
  dist(datExpr),
  method = "average"
)


# ------------------------------------------------------------
# Save PDF
# ------------------------------------------------------------

pdf(
  file.path(
    outpathcount,
    "sample_clustering_0.1.pdf"
  ),
  width = 10,
  height = 6
)


plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  cex.lab = 1,
  cex.axis = 1,
  cex.main = 1.2
)


dev.off()


# ------------------------------------------------------------
# Show in RStudio
# ------------------------------------------------------------

plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  cex.lab = 1,
  cex.axis = 1,
  cex.main = 1.2
)


# ============================================================
# 13. CREATE TRAITS
# ============================================================

traitData <- data.frame(
  
  xylella = as.numeric(
    metadata$treatment == "xylella"
  ),
  
  Catuai = as.numeric(
    metadata$cultivar == "Catuai"
  )
)


rownames(traitData) <-
  metadata$sample


# Ordenar igual que datExpr

traitData <- traitData[
  rownames(datExpr),
  ,
  drop = FALSE
]


stopifnot(
  identical(
    rownames(datExpr),
    rownames(traitData)
  )
)


print(
  traitData
)


# ============================================================
# 14. TRAIT COLORS
# ============================================================

traitColors <- numbers2colors(
  traitData,
  signed = FALSE
)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

pdf(
  file.path(
    outpathcount,
    "sample_dendrogram_traits_0.1.pdf"
  ),
  width = 10,
  height = 6
)


plotDendroAndColors(
  sampleTree,
  traitColors,
  groupLabels = names(traitData),
  main = "Sample dendrogram and trait heatmap"
)


dev.off()


# ------------------------------------------------------------
# Show
# ------------------------------------------------------------

plotDendroAndColors(
  sampleTree,
  traitColors,
  groupLabels = names(traitData),
  main = "Sample dendrogram and trait heatmap"
)


# ============================================================
# 15. SOFT THRESHOLD
# ============================================================

disableWGCNAThreads()


powers <- c(
  1:10,
  seq(
    12,
    20,
    2
  )
)


sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "unsigned",
  verbose = 5
)


print(
  sft$fitIndices[
    ,
    c(
      "Power",
      "SFT.R.sq",
      "mean.k."
    )
  ]
)


# ============================================================
# 16. SAVE SOFT THRESHOLD RESULTS
# ============================================================

write.table(
  
  sft$fitIndices,
  
  file = file.path(
    outpathcount,
    "soft_threshold_results_0.1.txt"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 17. SOFT THRESHOLD PLOT
# ============================================================

pdf(
  file.path(
    outpathcount,
    "soft_threshold_selection_0.1.pdf"
  ),
  width = 10,
  height = 5
)


par(
  mfrow = c(
    1,
    2
  )
)


# ------------------------------------------------------------
# Scale-free topology
# ------------------------------------------------------------

plot(
  
  sft$fitIndices$Power,
  
  -sign(
    sft$fitIndices$slope
  ) *
    sft$fitIndices$SFT.R.sq,
  
  xlab =
    "Soft threshold (power)",
  
  ylab =
    "Scale-free topology fit (R^2)",
  
  type = "n",
  
  main =
    "Scale independence"
)


text(
  
  sft$fitIndices$Power,
  
  -sign(
    sft$fitIndices$slope
  ) *
    sft$fitIndices$SFT.R.sq,
  
  labels =
    sft$fitIndices$Power,
  
  col =
    "red"
)


abline(
  h = 0.85,
  col = "blue",
  lty = 2
)


# ------------------------------------------------------------
# Mean connectivity
# ------------------------------------------------------------

plot(
  
  sft$fitIndices$Power,
  
  sft$fitIndices$mean.k.,
  
  xlab =
    "Soft threshold (power)",
  
  ylab =
    "Mean connectivity",
  
  type =
    "n",
  
  main =
    "Mean connectivity"
)


text(
  
  sft$fitIndices$Power,
  
  sft$fitIndices$mean.k.,
  
  labels =
    sft$fitIndices$Power,
  
  col =
    "red"
)


dev.off()


# ============================================================
# 18. SHOW SOFT THRESHOLD PLOT
# ============================================================

par(
  mfrow = c(
    1,
    2
  )
)


plot(
  
  sft$fitIndices$Power,
  
  -sign(
    sft$fitIndices$slope
  ) *
    sft$fitIndices$SFT.R.sq,
  
  xlab =
    "Soft threshold (power)",
  
  ylab =
    "Scale-free topology fit (R^2)",
  
  type =
    "n",
  
  main =
    "Scale independence"
)


text(
  
  sft$fitIndices$Power,
  
  -sign(
    sft$fitIndices$slope
  ) *
    sft$fitIndices$SFT.R.sq,
  
  labels =
    sft$fitIndices$Power,
  
  col =
    "red"
)


abline(
  h = 0.85,
  col = "blue",
  lty = 2
)


plot(
  
  sft$fitIndices$Power,
  
  sft$fitIndices$mean.k.,
  
  xlab =
    "Soft threshold (power)",
  
  ylab =
    "Mean connectivity",
  
  type =
    "n",
  
  main =
    "Mean connectivity"
)


text(
  
  sft$fitIndices$Power,
  
  sft$fitIndices$mean.k.,
  
  labels =
    sft$fitIndices$Power,
  
  col =
    "red"
)


# Regresar a una sola gráfica

par(
  mfrow = c(
    1,
    1
  )
)


# ============================================================
# 19. SELECT SOFT POWER
# ============================================================

r_sq_threshold <- 0.85


signed_R2 <-
  -sign(
    sft$fitIndices$slope
  ) *
  sft$fitIndices$SFT.R.sq


candidates <-
  sft$fitIndices$Power[
    signed_R2 >=
      r_sq_threshold
  ]


if (length(candidates) > 0) {
  
  softPower <- min(
    candidates
  )
  
} else {
  
  warning(
    "No power reached R^2 >= ",
    r_sq_threshold,
    "; using the highest power tested."
  )
  
  
  softPower <- max(
    sft$fitIndices$Power
  )
}


cat(
  "Selected soft-thresholding power:",
  softPower,
  "\n"
)


# ============================================================
# 20. ADJACENCY MATRIX
# ============================================================

adj <- adjacency(
  datExpr,
  power = softPower,
  type = "unsigned"
)


cat(
  "Adjacency matrix:",
  dim(adj),
  "\n"
)


# ============================================================
# 21. TOPOLOGICAL OVERLAP MATRIX
# ============================================================

TOM <- TOMsimilarity(
  adj,
  TOMType = "unsigned",
  verbose = 0
)


dissTOM <- 1 - TOM


dimnames(TOM) <- list(
  colnames(datExpr),
  colnames(datExpr)
)


# ============================================================
# 22. GENE CLUSTERING
# ============================================================

geneTree <- hclust(
  as.dist(dissTOM),
  method = "average"
)


# ============================================================
# 23. DYNAMIC TREE CUT
# ============================================================

min_module_size <- 30


dynamicMods <- cutreeDynamic(
  
  dendro = geneTree,
  
  distM = dissTOM,
  
  deepSplit = 2,
  
  pamRespectsDendro = FALSE,
  
  minClusterSize = min_module_size
)


dynamicColors <- labels2colors(
  dynamicMods
)


cat(
  "\nModules before merging:\n"
)


print(
  table(dynamicColors)
)


# ============================================================
# 24. MERGE SIMILAR MODULES
# ============================================================

merge_cut_height <- 0.25


merged <- mergeCloseModules(
  
  datExpr,
  
  dynamicColors,
  
  cutHeight =
    merge_cut_height,
  
  verbose = 0
)


mergedColors <-
  merged$colors


mergedMEs <-
  merged$newMEs


cat(
  "\nModules after merging:\n"
)


print(
  table(mergedColors)
)


# ============================================================
# 25. GENE DENDROGRAM + MODULE COLORS
# ============================================================

pdf(
  file.path(
    outpathcount,
    "gene_dendrogram_modules_0.1.pdf"
  ),
  width = 12,
  height = 6
)


plotDendroAndColors(
  
  geneTree,
  
  cbind(
    dynamicColors,
    mergedColors
  ),
  
  c(
    "Before merging",
    "After merging"
  ),
  
  dendroLabels = FALSE,
  
  hang = 0.03,
  
  addGuide = TRUE,
  
  guideHang = 0.05,
  
  main =
    "Gene dendrogram and module colors"
)


dev.off()


# Show in RStudio

plotDendroAndColors(
  
  geneTree,
  
  cbind(
    dynamicColors,
    mergedColors
  ),
  
  c(
    "Before merging",
    "After merging"
  ),
  
  dendroLabels = FALSE,
  
  hang = 0.03,
  
  addGuide = TRUE,
  
  guideHang = 0.05,
  
  main =
    "Gene dendrogram and module colors"
)


# ============================================================
# 26. FINAL MODULE ASSIGNMENTS
# ============================================================

moduleColors <-
  mergedColors


MEs <-
  mergedMEs


gene_modules <- data.frame(
  
  gene =
    colnames(datExpr),
  
  module =
    moduleColors,
  
  stringsAsFactors = FALSE
)


write.table(
  
  gene_modules,
  
  file = file.path(
    outpathcount,
    "gene_module_assignment_0.1.txt"
  ),
  
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)


# ============================================================
# 27. MODULE-TRAIT CORRELATIONS
# ============================================================

moduleTraitCor <- cor(
  MEs,
  traitData,
  use = "p"
)


moduleTraitP <- corPvalueStudent(
  moduleTraitCor,
  nrow(datExpr)
)


cat(
  "\nModule-trait correlations:\n"
)


print(
  moduleTraitCor
)


cat(
  "\nModule-trait p-values:\n"
)


print(
  moduleTraitP
)


# ============================================================
# 28. SAVE MODULE-TRAIT TABLES
# ============================================================

write.table(
  
  moduleTraitCor,
  
  file = file.path(
    outpathcount,
    "module_trait_correlations_0.1.txt"
  ),
  
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


write.table(
  
  moduleTraitP,
  
  file = file.path(
    outpathcount,
    "module_trait_pvalues_0.1.txt"
  ),
  
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


# ============================================================
# 29. MODULE-TRAIT HEATMAP
# ============================================================

textMatrix <- paste0(
  
  signif(
    moduleTraitCor,
    2
  ),
  
  "\n(",
  
  signif(
    moduleTraitP,
    1
  ),
  
  ")"
)


dim(textMatrix) <-
  dim(moduleTraitCor)


pdf(
  file.path(
    outpathcount,
    "module_trait_heatmap_0.1.pdf"
  ),
  width = 6,
  height = 9
)


par(
  mar = c(
    6,
    8.5,
    3,
    2
  )
)


labeledHeatmap(
  
  Matrix =
    moduleTraitCor,
  
  xLabels =
    names(traitData),
  
  yLabels =
    names(MEs),
  
  ySymbols =
    names(MEs),
  
  colorLabels =
    FALSE,
  
  colors =
    blueWhiteRed(50),
  
  textMatrix =
    textMatrix,
  
  setStdMargins =
    FALSE,
  
  cex.text =
    0.7,
  
  zlim =
    c(-1, 1),
  
  main =
    "Module-trait relationships"
)


dev.off()


# Show

par(
  mar = c(
    6,
    8.5,
    3,
    2
  )
)


labeledHeatmap(
  
  Matrix =
    moduleTraitCor,
  
  xLabels =
    names(traitData),
  
  yLabels =
    names(MEs),
  
  ySymbols =
    names(MEs),
  
  colorLabels =
    FALSE,
  
  colors =
    blueWhiteRed(50),
  
  textMatrix =
    textMatrix,
  
  setStdMargins =
    FALSE,
  
  cex.text =
    0.7,
  
  zlim =
    c(-1, 1),
  
  main =
    "Module-trait relationships"
)


# ============================================================
# 30. MODULE MOST CORRELATED WITH XYLELLA
# ============================================================

xylella_cor <-
  moduleTraitCor[
    ,
    "xylella"
  ]


best_module_ME <-
  names(
    which.max(
      abs(
        xylella_cor
      )
    )
  )


best_module <-
  sub(
    "^ME",
    "",
    best_module_ME
  )


cat(
  "\nModule most correlated with xylella:",
  best_module,
  "\n"
)


cat(
  "Correlation:",
  round(
    xylella_cor[
      best_module_ME
    ],
    2
  ),
  "\n"
)


cat(
  "P-value:",
  signif(
    moduleTraitP[
      best_module_ME,
      "xylella"
    ],
    2
  ),
  "\n"
)


cat(
  "Number of genes:",
  sum(
    moduleColors ==
      best_module
  ),
  "\n"
)


# ============================================================
# 31. MODULE MEMBERSHIP AND GENE SIGNIFICANCE
# ============================================================

geneModuleMembership <-
  as.data.frame(
    
    cor(
      datExpr,
      MEs,
      use = "p"
    )
    
  )


geneTraitSignificance <-
  as.data.frame(
    
    cor(
      datExpr,
      traitData$xylella,
      use = "p"
    )
    
  )


colnames(
  geneTraitSignificance
) <- "GS.xylella"


# ============================================================
# 32. HUB GENES OF BEST MODULE
# ============================================================

inModule <-
  moduleColors ==
  best_module


hub_table <- data.frame(
  
  gene =
    colnames(datExpr)[inModule],
  
  MM =
    geneModuleMembership[
      inModule,
      best_module_ME
    ],
  
  GS =
    geneTraitSignificance[
      inModule,
      "GS.xylella"
    ]
)


hub_table <- hub_table[
  
  order(
    -abs(
      hub_table$MM
    )
  ),
  
  ,
  
  drop = FALSE
]


write.table(
  
  hub_table,
  
  file = file.path(
    outpathcount,
    paste0(
      "hub_genes_",
      best_module,
      "_0.1.txt"
    )
  ),
  
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)


print(
  head(
    hub_table,
    15
  )
)


# ============================================================
# 33. MM vs GS PLOT
# ============================================================

pdf(
  
  file.path(
    outpathcount,
    paste0(
      "MM_vs_GS_",
      best_module,
      "_0.1.pdf"
    )
  )
)


plot(
  
  abs(
    hub_table$MM
  ),
  
  abs(
    hub_table$GS
  ),
  
  xlab =
    paste(
      "Module Membership in",
      best_module
    ),
  
  ylab =
    "Gene Significance for xylella",
  
  main =
    paste(
      "MM vs. GS -",
      best_module,
      "module"
    ),
  
  col =
    best_module,
  
  pch =
    19
)


dev.off()


# Show

plot(
  
  abs(
    hub_table$MM
  ),
  
  abs(
    hub_table$GS
  ),
  
  xlab =
    paste(
      "Module Membership in",
      best_module
    ),
  
  ylab =
    "Gene Significance for xylella",
  
  main =
    paste(
      "MM vs. GS -",
      best_module,
      "module"
    ),
  
  col =
    best_module,
  
  pch =
    19
)


# ============================================================
# 34. LOAD ANNOTATION
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
# 35. ONE ANNOTATION ROW PER GENE
# ============================================================

gene_annotation <- annotation[
  
  !duplicated(
    annotation$gene_id
  ),
  
  c(
    "gene_id",
    "Description",
    "Preferred_name"
  )
]


# ============================================================
# 36. ANNOTATE HUB GENES
# ============================================================

hub_table_annot <- merge(
  
  hub_table,
  
  gene_annotation,
  
  by.x =
    "gene",
  
  by.y =
    "gene_id",
  
  all.x =
    TRUE
)


hub_table_annot <- hub_table_annot[
  
  order(
    -abs(
      hub_table_annot$MM
    )
  ),
  
  ,
  
  drop = FALSE
]


write.table(
  
  hub_table_annot,
  
  file = file.path(
    
    outpathcount,
    
    paste0(
      "hub_genes_",
      best_module,
      "_annotated_0.1.txt"
    )
  ),
  
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)


print(
  head(
    hub_table_annot,
    10
  )
)


# ============================================================
# 37. EXPORT BEST MODULE TO CYTOSCAPE
# ============================================================

cytoscape_threshold <- 0.02


modProbes <- colnames(datExpr)[
  inModule
]


modTOM <- TOM[
  inModule,
  inModule
]


dimnames(modTOM) <- list(
  modProbes,
  modProbes
)


node_annot <-
  gene_annotation$Description[
    
    match(
      modProbes,
      gene_annotation$gene_id
    )
    
  ]


node_annot[
  is.na(node_annot)
] <- "No annotation"


cyt <- exportNetworkToCytoscape(
  
  modTOM,
  
  edgeFile = file.path(
    
    outpathcount,
    
    paste0(
      "CytoscapeInput-edges-",
      best_module,
      "_0.1.txt"
    )
  ),
  
  nodeFile = file.path(
    
    outpathcount,
    
    paste0(
      "CytoscapeInput-nodes-",
      best_module,
      "_0.1.txt"
    )
  ),
  
  weighted = TRUE,
  
  threshold =
    cytoscape_threshold,
  
  nodeNames =
    modProbes,
  
  nodeAttr =
    node_annot
)


cat(
  "Edges exported:",
  nrow(cyt$edgeData),
  "\n"
)


cat(
  "Nodes exported:",
  nrow(cyt$nodeData),
  "\n"
)


# ============================================================
# 38. COMPLETE WGCNA NETWORK FOR CYTOSCAPE
# ============================================================

all_genes <-
  colnames(datExpr)


dimnames(TOM) <- list(
  all_genes,
  all_genes
)


# Mantener las 10 conexiones TOM
# más fuertes de cada gen

top_n <- 10


edge_list <- do.call(
  
  rbind,
  
  lapply(
    
    seq_along(all_genes),
    
    function(i) {
      
      weights <- TOM[
        i,
      ]
      
      
      # No conectar gen consigo mismo
      
      weights[i] <- NA
      
      
      available <-
        sum(
          !is.na(weights)
        )
      
      
      strongest <- order(
        
        weights,
        
        decreasing = TRUE,
        
        na.last = NA
        
      )[
        1:min(
          top_n,
          available
        )
      ]
      
      
      data.frame(
        
        SourceNode =
          all_genes[i],
        
        TargetNode =
          all_genes[strongest],
        
        weight =
          weights[strongest],
        
        stringsAsFactors =
          FALSE
      )
    }
  )
)


# ============================================================
# 39. REMOVE DUPLICATE EDGES
# ============================================================

edge_list$key <- apply(
  
  edge_list[
    ,
    c(
      "SourceNode",
      "TargetNode"
    )
  ],
  
  1,
  
  function(x) {
    
    paste(
      sort(x),
      collapse = "__"
    )
    
  }
)


edge_list <- edge_list[
  !duplicated(edge_list$key),
  ,
  drop = FALSE
]


edge_list$key <- NULL


edge_list$direction <-
  "undirected"


# ============================================================
# 40. SAVE ALL EDGES
# ============================================================

write.table(
  
  edge_list,
  
  file = file.path(
    outpathcount,
    "Cytoscape_ALL_CLUSTERS_edges_0.1.txt"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 41. CREATE ALL NODE TABLE
# ============================================================

node_table <- data.frame(
  
  gene =
    all_genes,
  
  module =
    moduleColors,
  
  stringsAsFactors =
    FALSE
)


write.table(
  
  node_table,
  
  file = file.path(
    outpathcount,
    "Cytoscape_ALL_CLUSTERS_nodes_0.1.txt"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "Genes:",
  nrow(node_table),
  "\n"
)


cat(
  "Connections:",
  nrow(edge_list),
  "\n"
)


cat(
  "\nGenes by module:\n"
)


print(
  table(
    node_table$module
  )
)


# ============================================================
# 42. GENES IN greenyellow MODULE
# ============================================================

if ("greenyellow" %in% node_table$module) {
  
  
  greenyellow_genes <- node_table$gene[
    
    node_table$module ==
      "greenyellow"
    
  ]
  
  
  cat(
    "\nGenes in greenyellow:",
    length(greenyellow_genes),
    "\n"
  )
  
  
  print(
    greenyellow_genes
  )
  
  
  write.table(
    
    greenyellow_genes,
    
    file = file.path(
      outpathcount,
      "genes_module_greenyellow_0.1.txt"
    ),
    
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
} else {
  
  warning(
    "The greenyellow module was not found."
  )
}


# ============================================================
# 43. LOAD SIGNIFICANT DESeq2 GENES
#
# CR95
# Xylella vs saline
# padj < 0.1
# ============================================================

deseq_sig <- read.table(
  
  file.path(
    deseq_path,
    "DESeq2_CR95_xylella_vs_saline_SIGNIFICANT_0.1.txt"
  ),
  
  header = TRUE,
  
  sep = "\t",
  
  check.names = FALSE,
  
  stringsAsFactors = FALSE
)


# ============================================================
# 44. CHECK DESeq2 COLUMNS
# ============================================================

required_columns <- c(
  "Geneid",
  "log2FoldChange",
  "pvalue",
  "padj"
)


missing_columns <- setdiff(
  required_columns,
  colnames(deseq_sig)
)


if (length(missing_columns) > 0) {
  
  stop(
    "Missing columns in DESeq2 file: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ============================================================
# 45. ENSURE padj < 0.1
#
# Even though file already contains significant genes,
# this guarantees consistency.
# ============================================================

deseq_sig <- deseq_sig[
  
  !is.na(deseq_sig$padj) &
    deseq_sig$padj < 0.1,
  
  ,
  
  drop = FALSE
]


cat(
  "\nDESeq2 significant genes loaded:",
  nrow(deseq_sig),
  "\n"
)


cat(
  "UP:",
  sum(
    deseq_sig$log2FoldChange > 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "DOWN:",
  sum(
    deseq_sig$log2FoldChange < 0,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 46. ADD DESeq2 INFORMATION TO ALL WGCNA NODES
# ============================================================

cytoscape_DE <-
  node_table


idx <- match(
  
  cytoscape_DE$gene,
  
  deseq_sig$Geneid
)


# ============================================================
# 47. SIGNIFICANT YES / NO
# ============================================================

cytoscape_DE$DESeq2_significant <- ifelse(
  
  is.na(idx),
  
  "No",
  
  "Yes"
)


# ============================================================
# 48. LOG2 FOLD CHANGE
# ============================================================

cytoscape_DE$log2FoldChange <-
  deseq_sig$log2FoldChange[
    idx
  ]


# ============================================================
# 49. REGULATION
#
# UP   = log2FoldChange > 0
# DOWN = log2FoldChange < 0
#
# All matched genes already have padj < 0.1
# ============================================================

cytoscape_DE$regulation <-
  "Not_significant"


cytoscape_DE$regulation[
  
  !is.na(idx) &
    cytoscape_DE$log2FoldChange > 0
  
] <- "Up"


cytoscape_DE$regulation[
  
  !is.na(idx) &
    cytoscape_DE$log2FoldChange < 0
  
] <- "Down"


# ============================================================
# 50. ADD P-VALUE AND PADJ
# ============================================================

cytoscape_DE$pvalue <-
  deseq_sig$pvalue[
    idx
  ]


cytoscape_DE$padj <-
  deseq_sig$padj[
    idx
  ]


# ============================================================
# 51. REVIEW DESeq2 GENES IN NETWORK
# ============================================================

cat(
  "\nRegulation in complete WGCNA network:\n"
)


print(
  table(
    cytoscape_DE$regulation
  )
)


cat(
  "\nRegulation by WGCNA module:\n"
)


print(
  table(
    cytoscape_DE$module,
    cytoscape_DE$regulation
  )
)


cat(
  "\nSignificant DESeq2 genes entering WGCNA network:",
  
  sum(
    cytoscape_DE$DESeq2_significant ==
      "Yes"
  ),
  
  "\n"
)


cat(
  "UP entering network:",
  
  sum(
    cytoscape_DE$regulation ==
      "Up"
  ),
  
  "\n"
)


cat(
  "DOWN entering network:",
  
  sum(
    cytoscape_DE$regulation ==
      "Down"
  ),
  
  "\n"
)


# ============================================================
# 52. SAVE COMPLETE CYTOSCAPE NODE TABLE
# WITH DESeq2 INFORMATION
# ============================================================

write.table(
  
  cytoscape_DE,
  
  file = file.path(
    outpathcount,
    "Cytoscape_ALL_CLUSTERS_DESeq2_CR95_0.1.txt"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 53. greenyellow MODULE + DESeq2
# ============================================================

if ("greenyellow" %in% cytoscape_DE$module) {
  
  
  greenyellow_table <- cytoscape_DE[
    
    cytoscape_DE$module ==
      "greenyellow",
    
    ,
    
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Save ALL genes in greenyellow
  # with DESeq information
  # ----------------------------------------------------------
  
  write.table(
    
    greenyellow_table,
    
    file = file.path(
      outpathcount,
      "greenyellow_ALL_genes_DESeq2_CR95_0.1.txt"
    ),
    
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Only significant DESeq2 genes
  # padj < 0.1
  # ----------------------------------------------------------
  
  greenyellow_DE <- greenyellow_table[
    
    greenyellow_table$DESeq2_significant ==
      "Yes",
    
    ,
    
    drop = FALSE
  ]
  
  
  cat(
    "\nSignificant DESeq2 genes in greenyellow:",
    nrow(greenyellow_DE),
    "\n"
  )
  
  
  cat(
    "UP:",
    sum(
      greenyellow_DE$regulation ==
        "Up"
    ),
    "\n"
  )
  
  
  cat(
    "DOWN:",
    sum(
      greenyellow_DE$regulation ==
        "Down"
    ),
    "\n"
  )
  
  
  print(
    greenyellow_DE
  )
  
  
  # ----------------------------------------------------------
  # Save significant greenyellow genes
  # ----------------------------------------------------------
  
  write.table(
    
    greenyellow_DE,
    
    file = file.path(
      outpathcount,
      "DE_genes_module_greenyellow_CR95_0.1.txt"
    ),
    
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
}


# ============================================================
# 54. FINAL SUMMARY
# ============================================================

cat(
  "\n=============================================\n"
)


cat(
  "WGCNA ANALYSIS COMPLETE\n"
)


cat(
  "WGCNA network: top 5000 variable genes\n"
)


cat(
  "DESeq2 overlay: CR95 Xylella vs saline\n"
)


cat(
  "DESeq2 significance: padj < 0.1\n"
)


cat(
  "UP   = log2FoldChange > 0\n"
)


cat(
  "DOWN = log2FoldChange < 0\n"
)


cat(
  "All generated output files end in _0.1\n"
)


cat(
  "=============================================\n"
)

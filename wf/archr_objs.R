library(ArchR)
library(ggplot2)
library(harmony)
library(patchwork)
library(Seurat)

# globals ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

project_name <- args[1]
genome <- args[2]
tile_size <- as.integer(args[3])
min_tss <- as.numeric(args[4])
min_frags <- as.integer(args[5])
lsi_iterations <- as.integer(args[6])

for (i in strsplit(args[7], ",")) {
  lsi_resolution <- as.numeric(i)
}
for (i in strsplit(args[8], ",")) {
  lsi_varfeatures <- as.integer(i)
  }
for (i in strsplit(args[9], ",")) {
  clustering_resolution <- as.numeric(i)
}

umap_mindist <- as.numeric(args[10])

runs <- strsplit(args[11:length(args)], ",")
inputs <- c()
for (run in runs) {
  inputs[run[1]] <- run[2]
  }

out_dir <- paste0(project_name, "_ArchRProject")

# functions --------------------------------------------------------------------

build_atlas_seurat_object <- function(
  run_id,
  matrix,
  metadata,
  spatial_path) {
  # Prepare and combine gene matrix, metadata, and image for seurat object
  # for runs within a project.

  image <- Read10X_Image(
    image.dir = spatial_path,
    filter.matrix = TRUE
  )
  metadata <- metadata[metadata$Sample == run_id, ]

  matrix <- matrix[, c(grep(pattern = run_id, colnames(matrix)))]
  matrix@Dimnames[[2]] <- metadata@rownames

  object <- CreateSeuratObject(
    counts = matrix,
    assay  = "Spatial",
    meta.data = as.data.frame(metadata)
  )
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- "Spatial"
  object[["slice1"]] <- image
  return(object)
}

spatial_plot <- function(seurat_object, name) {
  clusters <- sort(unique(seurat_object$Clusters))
  colors <- ArchRPalettes$stallion2[seq_len(length(clusters))]
  names(colors) <- clusters
  SpatialDimPlot(
    seurat_object,
    group.by = "Clusters",
    pt.size.factor = 1,
    cols = colors,
    stroke = 0) +
  ggtitle(name) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 10))
}

feature_plot <- function(seurat_obj, feature, name) {
  SpatialFeaturePlot(
    object = seurat_obj,
    features = feature,
    alpha = c(0.2, 1),
    pt.size.factor = 1
    ) +
  ggtitle(paste0(feature, " : ", name)) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15)
  )
}

# create archr project --------------------------------------------------------

addArchRGenome(genome)
addArchRThreads(threads = 24)

arrow_files <- createArrowFiles(
   inputFiles = inputs,
   sampleNames = names(inputs),
   minTSS = min_tss,
   minFrags = min_frags,
   maxFrags = 1e+07,
   addTileMat = TRUE,
   addGeneScoreMat = TRUE,
   offsetPlus = 0,
   offsetMinus = 0,
   TileMatParams = list(tileSize = tile_size)


proj <- ArchRProject(
  ArrowFiles = arrow_files,
  outputDirectory = out_dir
)

# Add an additional Conditions column
for (run in runs) {
  proj$Condition[proj$Sample == run[1]] <- run[3]
}

# Filter on-tissue
all_ontissue <- c()
for (run in runs) {
  positions <- read.csv(run[4], header = FALSE)
  positions$V1 <- paste(run[1], "#", positions$V1, "-1", sep = "")
  on_tissue <- positions$V1 [which(positions$V2 == 1)]
  all_ontissue <- c(all_ontissue, on_tissue)
}
proj <- proj[proj$cellNames %in% all_ontissue]

# iterate plotting ------------------------------------------------------------

# make dataframe with  Cartesian Product of three parameter lists
parameter_set <- expand.grid(
  lsi_resolution,
  lsi_varfeatures,
  clustering_resolution
)
print(parameter_set)

# init 'dict' to store dimplots, vector for umap plots
umapplots <- c()
dimplots <- list()

for (row in 1:nrow(parameter_set)) {

  set <- parameter_set[c(row),c(1,2,3)]
  lsi_resolution_i <- set[[1]]
  varfeatures_i <- set[[2]]
  clustering_resolution_i <- set[[3]]

  print(c(lsi_resolution_i, varfeatures_i, clustering_resolution_i))

  # work with a copy of the original project
  proj_i <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = lsi_iterations,
    clusterParams = list(
      resolution = c(lsi_resolution_i),
      sampleCells = 10000,
      n.start = 10
    ),
    varFeatures = varfeatures_i,
    dimsToUse = 1:30,
    force = TRUE
  )
  if (length(runs) > 1) {
    proj_i <- addHarmony(
      ArchRProj = proj_i,
      reducedDims = "IterativeLSI",
      name = "Harmony",
      groupBy = "Sample",
      force = TRUE
    )
    name <- "Harmony"
  } else {
    name <- "IterativeLSI"
  }
  proj_i <- addClusters(
    input = proj_i,
    reducedDims = name,
    method = "Seurat",
    name = "Clusters",
    resolution = c(clustering_resolution_i),
    force = TRUE
  )
  proj_i <- addUMAP(
    ArchRProj = proj_i,
    reducedDims = name,
    name = "UMAP",
    nNeighbors = 30,
    minDist = umap_mindist,
    metric = "cosine",
    force = TRUE
  )

  # plot umaps by sample and cluster
  p1 <- plotEmbedding(
    ArchRProj = proj_i,
    colorBy = "cellColData",
    name = "Sample",
    embedding = "UMAP"
  ) +
  ggtitle(
    paste(
      "colored by Sample",
      lsi_resolution_i,
      varfeatures_i,
      clustering_resolution_i
    )
  ) +
  theme(plot.title = element_text(size = 10)) +
  theme(legend.key.size = unit(.5, "cm")) +
  theme(legend.text=element_text(size = 6)) +
  guides(colour = guide_legend(
    override.aes = list(size = 2, alpha = 1),
    nrow = 2)
  )

  p2 <- plotEmbedding(
    ArchRProj = proj_i,
    colorBy = "cellColData",
    name = "Clusters",
    embedding = "UMAP"
  ) +
  ggtitle(
    paste(
      "colored by Cluster",
      lsi_resolution_i,
      varfeatures_i,
      clustering_resolution_i
    )
  ) +
  theme(plot.title = element_text(size = 10)) +
  theme(legend.key.size = unit(.5, "cm")) +
  theme(legend.text=element_text(size = 6)) +
  guides(colour = guide_legend(
    override.aes = list(size = 2, alpha = 1),
    nrow = 2)
  )
  umapplots[[row]] <- p1 + p2

  proj_i <- addImputeWeights(proj_i)

  # create metadata object for Seurat object
  metadata <- getCellColData(ArchRProj = proj_i)
  rownames(metadata) <- str_split_fixed(
  str_split_fixed(
    row.names(metadata),
    "#",
    2)[, 2],
  "-",
  2)[, 1]
  metadata["log10_nFrags"] <- log(metadata$nFrags)

  # create gene matrix for Seurat object
  gene_matrix <- getMatrixFromProject(
    ArchRProj = proj_i,
    useMatrix = "GeneScoreMatrix"
  )
  matrix <- imputeMatrix(
    mat = assay(gene_matrix),
    imputeWeights = getImputeWeights(proj_i)
  )
  gene_row_names <- gene_matrix@elementMetadata$name
  rownames(matrix) <- gene_row_names

  seurat_objs <- c()
  for (run in runs) {

    obj <- build_atlas_seurat_object(
      run_id = run[1],
      matrix = matrix,
      metadata = metadata,
      spatial_path = run[5]
    )
    seurat_objs <- c(seurat_objs, obj)

    p1 <- spatial_plot(
      obj,
      name = paste(
        run[1],
        lsi_resolution_i,
        varfeatures_i,
        clustering_resolution_i
      )
    )
    dimplots[[run[1]]][[row]] <- p1
    }
}

# save umap plots in a single pdf
pdf("umap_plots.pdf")
for (i in seq_along((umapplots))) {
  print(umapplots[[i]])
}
dev.off()

# save spatialdim plots in a single pdf
pdf("spatialdim_plots.pdf")
for (i in seq_along((dimplots))) {
  print(grid.arrange(grobs = dimplots[[i]], ncol = 2))
}
dev.off()

# save qc plots in a single pdf
pdf("qc_plots.pdf")
for (obj in seurat_objs) {
  name <- unique(obj@meta.data[["Sample"]])
  nfrags_plot <- feature_plot(obj, "log10_nFrags", name)
  tss_plot <- feature_plot(obj, "TSSEnrichment", name)

    print(nfrags_plot)
    print(tss_plot)
    par(newpage = TRUE)

}
dev.off()

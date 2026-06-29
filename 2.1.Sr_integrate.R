
suppressMessages(library(Seurat))
# suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(tidyverse))
# suppressMessages(library(ConfigParser))
suppressMessages(library(reticulate))
suppressMessages(library(KernSmooth))
suppressMessages(library(DoubletFinder))
# suppressMessages(library(DoubletDecon))
suppressMessages(library(fields))
suppressMessages(library(psych))
suppressMessages(library(pheatmap))
#suppressMessages(library(scater))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(scales))
suppressMessages(library(cowplot))

opt <- list(
  id = "PF",
  outdir = "../result_review/PD_integrated_rds/", 
  path = "./lib_PD.csv", 
  genome = "GRCh38",
  min_cells = 3,
  high_nGene = 30000,
  integrate = "Seurat",
  resolution = 0.6,
  Nfeatures = 2000,
  y_high_dispersion_cutoff = 1.0,
  y_low_dispersion_cutoff = 0.1,
  x_high_dispersion_cutoff = 5,
  x_low_dispersion_cutoff = 0.125,
  high_percent.mito = 10.0,
  high_percent.HB = 5.0,
  pca = 20,
  sct = "SCT",
  low_nGene = 200,
  qc = FALSE,
  cycle = FALSE,
  gene_positive = NULL,
  gene_dispersion = NULL,
  gene_dispersion = NULL
)

###configure
ANNO <- '/disk/user/liaoshuilin/Pipeline-bai/Cr/lib/Anno/'

CC = c("#DC143C", "#0000FF", "#20B2AA", "#FFA500", "#9370DB", "#98FB98", "#F08080", "#1E90FF", "#7CFC00", "#FFFF00", "#808000", "#FF00FF", "#FA8072", "#7B68EE", "#9400D3", "#800080", "#A0522D", "#D2B48C", "#D2691E", "#87CEEB", "#40E0D0", "#5F9EA0", "#FF1493", "#0000CD", "#008B8B", "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A", "#4682B4", "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE", "#FF6347", "#6A5ACD", "#9932CC", "#8B008B", "#8B4513", "#DEB887")

## === function ===
## ======

#mkdir
dirMK <- function(Dir) {
  if (!dir.exists(Dir)) {
    dir.create(Dir)
  }
}

#convert pdf2png
f2g <- function(mypdf, mypng) {
  mypng <- sub('.png', '', mypng)
  system(paste0('/usr/bin/pdftoppm -png ', mypdf, ' ', mypng, ' -singlefile'))
}

#print log information
prdate <- function(LOG) {
  print (date())
  print (paste0(LOG, ' start...'))
}

##Integrate
#Create Object
CreateObj <- function() {
  libs <- read.table(file = opt$path, header = F, sep = ',')
  colnames(libs) <- c('Sample', 'Path')
  ob.list <- list()
  numsample <<- 0
  Allsamples <<- c()
  for (eachS in libs$Sample) {
    numsample <<- numsample + 1
    Allsamples <<- c(Allsamples, eachS)
    #		ob = paste('ob', eachS, sep = '_')
    ExprPath = libs[libs$Sample == eachS, 'Path']
    if (endsWith(ExprPath, '.csv')) {
      pbmc <- read.csv(ExprPath, header = T, row.names = 1, stringsAsFactors = FALSE, check.names = F)
    } else if (endsWith(ExprPath, '_matrix')) {
      pbmc <- Read10X(data.dir = ExprPath)
    } else {
      stop('Original File type is not correct!')
    }
    colnames(pbmc) <- str_replace_all(colnames(pbmc), '-1', paste0('-', numsample))
    ob <- CreateSeuratObject(counts = pbmc, project = eachS, min.cells = opt$min_cells)
    ob$stim <- eachS
    if (! is.null(opt$gene_positive)) {
      positive.genes_total <- readLines(gene_positive)
      positive.genes <- intersect(positive.genes_total, rownames(ob@assays$RNA@counts))
      ob[["percent.positive"]] <- PercentageFeatureSet(ob, features = positive.genes)
      ob <- subset(ob, subset = percent.positive == 0)
    }
    ob <- subset(ob, subset = nFeature_RNA > as.numeric(opt$low_nGene))
    ob <- NormalizeData(ob)
    ob <- FindVariableFeatures(ob, selection.method = "vst", nfeatures = opt$Nfeatures)
    ob.list[[eachS]] <- ob
  }
  return (ob.list)
}

#Anchors
AnchorParse <- function(object) {
  if (! is.null(opt$gene_anchors)) {
    return (readLines(opt$gene_anchors))
  } else {
    return (SelectIntegrationFeatures(object.list = object, nfeatures = as.numeric(opt$Nfeatures)))
  }
}

## === integrate ===
## ======

prdate('Seurat Integrate')
#1. Analydir
if(!exists(opt$outdir)){
  dirMK(opt$outdir)
}

for (eachDir in c('Anchors', 'Diff', 'Marker', 'CellsRatio')) {
  dirMK(paste0(opt$outdir, '/', eachDir))
}

#2. Gene Annotation
gene.names <- read.table(file = paste0(ANNO, '/', opt$genome, '_genes.tsv'), sep = '\t', header = F, stringsAsFactors = FALSE)
gene.names$V2 <- make.unique(gene.names$V2)
colnames(gene.names) <- c('Gene', 'Gene_name')

#3. Create Seurat object
prdate('3')
ob.list <- CreateObj()

#3.1 Integrate
prdate('3.1')
if (opt$sct == 'SCT') {
  combined <- lapply(X = ob.list, FUN = SCTransform)
  anchors <- AnchorParse(combined)
  combined <- PrepSCTIntegration(object.list = combined, anchor.features = anchors)
  combined <- FindIntegrationAnchors(object.list = combined, normalization.method = 'SCT', anchor.features = anchors, reference = 1)
  combined <- IntegrateData(anchorset = combined, normalization.method = 'SCT', dims = 1:as.numeric(opt$pca))
  if (! is.null(opt$gene_dispersion)) {
    high_varience_gene <- readLines(opt$gene_dispersion)
  } else {
    high_varience_gene <- VariableFeatures(combined)
  }
} else {
  anchors <- AnchorParse(ob.list)
  if (opt$sct == 'CCA') {
    combined <- FindIntegrationAnchors(object.list = ob.list, dims = 1:as.numeric(opt$pca), reference = 1)
  } else {
    combined <- FindIntegrationAnchors(object.list = ob.list, anchor.features = anchors, reference = 1)
  }
  combined <- IntegrateData(anchorset = combined, dims = 1:as.numeric(opt$pca))
  if (! is.null(opt$gene_dispersion)) {
    high_varience_gene <- readLines(opt$gene_dispersion)
  } else {
    high_varience_gene <- VariableFeatures(combined)
  }
  combined <- ScaleData(combined, features = high_varience_gene, verbose = FALSE)
}
anchors <- data.frame(anchors)
colnames(anchors) <- c('Gene_name')
write.table(anchors, file = paste0(opt$outdir, '/Anchors/', opt$id, '_anchors_gene.csv'), quote = F, row.names = F)
DefaultAssay(combined) <- 'integrated'
combined <- RunPCA(combined, npcs = as.numeric(opt$pca), verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:as.numeric(opt$pca))
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunTSNE(object = combined, dims.use = 1:as.numeric(opt$pca), do.fast = TRUE, check_duplicates = FALSE)
combined <- RunUMAP(object = combined, dims = 1:as.numeric(opt$pca), umap.method = 'umap-learn', metric = 'correlation', check_duplicates = FALSE)
saveRDS(combined, file = paste0(opt$outdir, '/', opt$id, '_combined.rds'))

#3.3 High_varience_gene
prdate('3.3')
high_varience_gene <- data.frame(high_varience_gene)
colnames(high_varience_gene)[1] <- 'Gene_name'
#	high_varience_gene <- left_join(high_varience_gene, gene.names, by = 'Gene_name')
write.table(high_varience_gene, file = paste0(opt$outdir, '/Anchors/', opt$id, '_high_varience_gene.csv'), sep = ',', quote = F, row.names = F)

count_bloor = "N"
if(count_bloor == "Y"){

  #4.1 Count
  prdate('4.1')
  DataCluster <- data.frame(Idents(combined))
  colnames(DataCluster) <- c('Cluster')
  write.table(DataCluster, file = paste0(opt$outdir, '/Anchors/', opt$id, '_cluster.csv'), row.names = T, col.names = T, quote = F, sep = ',')
  ##4.1.1 Samples cluster percent
  prdate('4.1.1')
  sample_clus <- data.frame(Idents(combined))
  sample_clus <- cbind(rownames(sample_clus), sample_clus)
  colnames(sample_clus) <- c('Barcode', 'Cluster')
  num <- numsample + 1
  b <- mutate(sample_clus, Sample = as.numeric(str_split(sample_clus$Barcode, '-', simplify = TRUE)[,2]))
  c <- b %>% 
    group_by(Cluster) %>% 
    dplyr::count(Sample) %>% 
    spread(key = Sample, value = n)
  print (head(c))
  colnames(c)[2:num] <- Allsamples
  bf <- c[,-1]
  rowsum <- rowSums(bf, na.rm = T)
  bf_r <- bf/rowsum
  bf_r <- cbind(as.numeric(rownames(bf_r)) - 1, bf_r)
  colnames(bf_r)[1] <- 'Cluster'
  colnames(bf_r)[2:num] <- Allsamples
  write.table(c, file = paste0(opt$outdir, '/CellsRatio/', opt$id, '_cluster_count.csv'), sep = ',', row.names = F, quote = F)
  write.table(bf_r, file = paste0(opt$outdir, '/CellsRatio/', opt$id, '_cluster_percent.csv'), sep = ',', row.names = F, quote = F)
  colour <- if (numsample > 40) hue_pal()(numsample) else CC[1:numsample]
  td <- gather(bf_r, key = 'Cluster Name', value = 'Cells Ratio', -Cluster)
  td[,1] <- factor(td[,1], levels = sort(as.numeric(bf_r$Cluster)))
  td[,2] <- factor(td[,2], levels = Allsamples)
  pdf(paste0(opt$outdir, '/CellsRatio/', opt$id, '_cluster_percent.pdf'))
  print (ggplot(td, aes(x = td[,1], y = td[,3], fill = td[,2])) + geom_bar(position = 'stack', stat = 'identity') + labs(x = 'Cluster', y = 'Cells Ratio') + theme(panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill ='transparent', color = 'transparent'), axis.text = element_text(color = 'black')) + scale_y_continuous(expand = c(0.001, 0.001)) + scale_fill_manual(values = colour) + guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol = 1, title = 'Sample')))
  dev.off()
  f2g(paste0(opt$outdir, '/CellsRatio/', opt$id, '_cluster_percent.pdf'), paste0(opt$outdir, '/CellsRatio/', opt$id, '_cluster_percent.png'))
  
  #4.1.2 Cluster Sample percent
  prdate('4.1.2')
  bf_c <- t(t(bf)/rowSums(t(bf), na.rm = T))
  row.names(bf_c) <- c(1:dim(bf_c)[1])
  bf_c <- as.data.frame(cbind((as.numeric(rownames(bf_c)) - 1), bf_c))
  colnames(bf_c)[1] <- "Cluster"
  write.table(bf_c, file = paste0(opt$outdir, '/CellsRatio/', opt$id, '_sample_percent.csv'), sep = ',', row.names = F, quote = F)
  cluster_number <- length(row.names(bf_c))
  colour <- if(cluster_number > 40) hue_pal()(cluster_number) else CC[1:cluster_number]
  td_c <- gather(bf_c, key = 'Sample Name', value = 'Cells Ratio', -Cluster)
  td_c[,1] <- factor(td_c[,1], levels = sort(as.numeric(bf_c$Cluster)))
  td_c[,2] <- factor(td_c[,2], levels = Allsamples)
  pdf(paste0(opt$outdir, '/CellsRatio/', opt$id, '_sample_percent.pdf'))
  print (ggplot(td_c, aes(x = td_c[,2], y = td_c[,3], fill = td_c[,1])) + geom_bar(position = 'stack', stat = 'identity') + labs(x = 'Sample', y = 'Cells Ratio') + theme(panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent', color = 'transparent'), axis.text = element_text(color = 'black')) + scale_y_continuous(expand = c(0.001, 0.001)) + scale_fill_manual(values = colour) + guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol = 1, title = 'Cluster')))
  dev.off()
  f2g(paste0(opt$outdir, '/CellsRatio/', opt$id, '_sample_percent.pdf'), paste0(opt$outdir, '/CellsRatio/', opt$id, '_sample_percent.png'))
}

#4.2 Tsne
prdate('4.2')
tsne.data <- combined@reductions$tsne@cell.embeddings
tsne.data <- cbind(rownames(tsne.data), tsne.data)
colnames(tsne.data)[1] <- 'Barcode'
write.table(tsne.data, file = paste0(opt$outdir, '/Anchors/', opt$id, '_tSNE.csv'), row.names = F, col.names = T, sep = ',', quote = F)
pdf(paste0(opt$outdir, '/Anchors/', opt$id, '_sample_tSNE.pdf'))
print (DimPlot(combined, reduction = 'tsne', group.by = 'stim'))
dev.off()
f2g(paste0(opt$outdir, '/Anchors/', opt$id, '_sample_tSNE.pdf'), paste0(opt$outdir, '/Anchors/', opt$id, '_sample_tSNE.png'))
pdf(paste0(opt$outdir, '/Anchors/', opt$id, '_cluster_tSNE.pdf'))
print (DimPlot(combined, reduction = "tsne", label = TRUE))
dev.off()
f2g(paste0(opt$outdir, '/Anchors/', opt$id, '_cluster_tSNE.pdf'), paste0(opt$outdir, '/Anchors/', opt$id, '_cluster_tSNE.png'))
  
#4.3 Umap
prdate('4.3')
umap.data <- combined@reductions$umap@cell.embeddings
umap.data <- cbind(rownames(umap.data), umap.data)
colnames(umap.data)[1] <- 'Barcode'
write.table(umap.data, file = paste0(opt$outdir, '/Anchors/', opt$id, '_UMAP.csv'), row.names = F, col.names = T, sep = ',', quote = F)
pdf(paste0(opt$outdir, '/Anchors/', opt$id, '_sample_UMAP.pdf'))
print (DimPlot(combined, reduction = 'umap', group.by = 'stim'))
dev.off()
f2g(paste0(opt$outdir, '/Anchors/', opt$id, '_sample_UMAP.pdf'), paste0(opt$outdir, '/Anchors/', opt$id, '_sample_UMAP.png'))
pdf(paste0(opt$outdir, '/Anchors/', opt$id, '_cluster_UMAP.pdf'))
print (DimPlot(combined, reduction = 'umap', label = TRUE))
dev.off()
f2g(paste0(opt$outdir, '/Anchors/', opt$id, '_cluster_UMAP.pdf'), paste0(opt$outdir, '/Anchors/', opt$id, '_cluster_UMAP.png'))
  
#4.4 Marker
prdate('4.4')
DefaultAssay(combined) <- 'RNA'
combined <- JoinLayers(combined)
combined@misc$markers <- FindAllMarkers(object = combined, only.pos = F, min.pct = 0.25, layer = "data")
combined.markers <- combined@misc$markers
colnames(combined.markers)[7] <- 'Gene_name'
combined.markers <- left_join(combined.markers, gene.names, by = 'Gene_name')
combined.markers <- combined.markers[,c(8,7,6,2,1,5,3,4)]
write.table(combined.markers,file = paste0(opt$outdir, '/Diff/', 'Cluster_diff.csv'), quote = F, row.names = F, col.names = T, sep = ',')
cluster <- unique(combined.markers$cluster)
for (i in cluster) {
  data <- filter(combined.markers, cluster == i)
  write.table(data, paste0(opt$outdir, '/Diff/', 'Cluster_', i, '_diff.csv'), quote = F, row.names = F, col.names = T, sep = ',')
  dataSig <- filter(combined.markers, cluster == i, p_val_adj < 0.05, avg_log2FC > 0)
  write.table(dataSig, paste0(opt$outdir, '/Diff/', 'Cluster_', i, '_diff_significant.csv'), quote = F, row.names = F, col.names = T, sep = '\t')
}
  
#4.4.1 Vio and Tsne Figure
prdate('4.4.1')
top_markers <- combined.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
# colnames(top_markers)[2] <- 'gene'
clus <- unique(top_markers$cluster)
for (each in clus) {
  top_clus <- top_markers[top_markers$cluster == each,]
  genes <- top_clus$Gene_name
  genes = genes[!stringr::str_starts(genes, '-')]
  pdf(paste0(opt$outdir, '/Marker/', opt$id, '_Cluster_', each, '_violin.pdf'), height = 5, width = 6)
  print (VlnPlot(combined, genes, ncol = 2, pt.size = 0) + xlab('Cluster') + ylab('log(UMI)'))
  dev.off()
  f2g(paste0(opt$outdir, '/Marker/', opt$id, '_Cluster_', each, '_violin.pdf'), paste0(opt$outdir, '/Marker/', opt$id, '_Cluster_', each, '_violin.png'))
  pdf(paste0(opt$outdir, '/Marker/', opt$id, '_Cluster_', each, '_tsne.pdf'), height = 6, width = 7)
  print (FeaturePlot(object = combined, genes, cols = c('grey', 'blue'), reduction = 'tsne'))
  dev.off()
  f2g(paste0(opt$outdir, '/Marker/', opt$id, '_Cluster_', each, '_tsne.pdf'), paste0(opt$outdir, '/Marker/', opt$id, '_Cluster_', each, '_tsne.png'))
}
  
#4.4.2 Heatmap
prdate('4.4.2')
top_genes <- top_markers$Gene_name
combined <- ScaleData(object = combined, features = top_genes)
pdf(paste0(opt$outdir, '/Marker/', opt$id, '_Heatmap.pdf'), height = 10, width = 15)
print (DoHeatmap(object = combined,  features = top_genes))
dev.off()
f2g(paste0(opt$outdir, '/Marker/', opt$id, '_Heatmap.pdf'), paste0(opt$outdir, '/Marker/', opt$id, '_Heatmap.png'))

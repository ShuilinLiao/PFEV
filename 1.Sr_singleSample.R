suppressMessages(library(Seurat))
suppressMessages(library(optparse))
suppressMessages(library(ConfigParser))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(tidyverse))
suppressMessages(library(reticulate))
suppressMessages(library(KernSmooth))
suppressMessages(library(DoubletFinder))
suppressMessages(library(fields))
suppressMessages(library(psych))
suppressMessages(library(pheatmap))
# suppressMessages(library(scater))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(scales))
suppressMessages(library(cowplot))
suppressMessages(library(Matrix))
suppressMessages(library(SoupX))
suppressMessages(library(DropletUtils))

rm(list = ls())
opt <- list(outdir = "../result_review/single_rds",
            path = "../data/GSE130888/LTPD4", # alterable
            id = "LTPD4", # alterable
            genome = "GRCh38",
            min_cells = 3,
            high_nGene = 30000,
            integrate = NULL,
            resolution = 0.4,
            Nfeatures = 2000,
            y_high_dispersion_cutoff = 1.0,
            y_low_dispersion_cutoff = 0.1,
            x_high_dispersion_cutoff = 5,
            x_low_dispersion_cutoff = 0.125,
            high_percent.mito = 30.0,
            high_percent.HB = 5.0,
            pca = 20,
            sct = "SCT",
            low_nGene = 200,
            doublet = "DoubletFinder",
            qc = FALSE,
            cycle = FALSE)

###configure
ANNO <- '/disk/user/liaoshuilin/Pipeline-bai/Cr/lib/Anno/'


###Local configure
CC = c("#DC143C", "#0000FF", "#20B2AA", "#FFA500", "#9370DB", "#98FB98", "#F08080", "#1E90FF", "#7CFC00", "#FFFF00", "#808000", "#FF00FF", "#FA8072", "#7B68EE", "#9400D3", "#800080", "#A0522D", "#D2B48C", "#D2691E", "#87CEEB", "#40E0D0", "#5F9EA0", "#FF1493", "#0000CD", "#008B8B", "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A", "#4682B4", "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE", "#FF6347", "#6A5ACD", "#9932CC", "#8B008B", "#8B4513", "#DEB887")


###function
#mkdir
dirMK <- function(Dir) {
  if (!dir.exists(Dir)) {
    dir.create(Dir)
  }
}

#prepare dataset
rds_prepare <- function(infile) {
  # infile = opt$path
  pbmc.data <- Read10X(data.dir = infile)
  colnames(pbmc.data) <- paste(colnames(pbmc.data), '-', 1, sep = '')
  pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = opt$min_cells, project = opt$id)
  rdata <- list(pbmc, nrow(pbmc.data), ncol(pbmc.data))
  return (rdata)
}

#99% genes
perc99 <- function(Data) {
  # Data = pbmc
  nGene_99 <- as.integer(quantile(Data@meta.data$nFeature_RNA, probs = 0.99))
  if (nGene_99 > 25000) {
    print('you need to set nGene manually, nGene_99>15000!')
  } else {
    for (j in seq(500, 25000, 500)) {
      nj <- j + 500
      if (nGene_99 < nj){
        nGene_99 <- nj
        break
      }
    }
  }
  if (as.numeric(opt$high_nGene) == 30000) {
    high_nGene <- nGene_99
  } else {
    high_nGene <- as.numeric(opt$high_nGene)
  }
  return (high_nGene)
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

###################################Single Analysis########################################

gene.names <- read.table(file = paste0(paste0(opt$path), '/', 'features.tsv.gz'), sep ='\t', header = F, stringsAsFactors = FALSE)
gene.names <- gene.names[,1:2]
# gene.names$V2 <<- make.unique(gene.names$V2)
colnames(gene.names) <- c('Gene','Gene_name')

######Step1 Prepare rds dataset
dirMK(opt$outdir)
dirMK(paste0(opt$outdir, "/", opt$id))
for ( eachDir in c('QC', 'PCA', 'DIMENSION', 'DIFF', 'MARKER') ){
  dirMK(paste0(opt$outdir, "/", opt$id, '/', eachDir))
}

#1.1 Prepare
prdate('1.1')
RDSF <- rds_prepare(opt$path)
pbmc <- RDSF[[1]]
gene_number <- RDSF[[2]]
cell_number <- RDSF[[3]]
nGene_median <- as.integer(median(pbmc@meta.data$nFeature_RNA))
pbmc_raw <- pbmc

#1.2 compute 99 persent of nGene(integer) and percent(two decimal)
prdate('1.2')
nGene_99 <- as.integer(quantile(pbmc@meta.data$nFeature_RNA,probs = 0.99))

#1.3 choose nGene_99 as n*500 
prdate('1.3')
high_nGene <- perc99(pbmc)


######Step2 QC
#2.1 makedir
prdate('2.1')
dirMK(paste0(opt$outdir, '/QC'))


#2.2 summarize the base information: gene and umi per cell
prdate('2.2')
meta_data <- data.frame(pbmc@meta.data)
meta_data <- cbind(rownames(meta_data), meta_data)
colnames(meta_data)[1] <- 'Barcode'
write.table(meta_data, file = paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI.csv'), sep = ',', quote = F, row.names = F)

#2.3 mitochondrion and erythrocyte percentage and filtering
prdate('2.3')
if (opt$genome %in% c('GRCh38', 'mm10', 'Rattus_norvegicus')){
  #configure
  config <- ConfigParser$new(Sys.getenv(), optionxform = identity)
  config$read(paste0(ANNO, '/Seurat_Mito_Hb.INI'))
  
  #mitochondrion and erythrocyte genes
  mito.genes_total <- as.vector(unlist(str_split(config$get(opt$genome, NA, 'MITO'), ',')))
  HB.genes_total <- as.vector(unlist(str_split(config$get(opt$genome, NA, 'HB'), ',')))
  #consensus mitochondrion genes
  mito_m  <- match(mito.genes_total, rownames(pbmc@assays$RNA))
  mito.genes <- rownames(pbmc@assays$RNA)[mito_m]
  mito.genes <- mito.genes[!is.na(mito.genes)]
  pbmc[['percent.mito']] <- PercentageFeatureSet(pbmc, features = mito.genes)
  #consensus erythrocyte genes
  HB_m <- match(HB.genes_total, rownames(pbmc@assays$RNA))
  HB.genes <- rownames(pbmc@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]
  pbmc[['percent.HB']] <- PercentageFeatureSet(pbmc, features = HB.genes)
  
  #Figure for mitochondrion and erythrocyte
  pdf(paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI_mito_HB.pdf'), width = 10)
  print(VlnPlot(pbmc, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mito','percent.HB'), ncol = 4))
  dev.off()
  f2g(paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI_mito_HB.pdf'), paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI_mito_HB.png'))
  
  #Figure for relation among UMI and mitochondrion and erythrocyte
  p1 <- FeatureScatter(object = pbmc, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
  p2 <- FeatureScatter(object = pbmc, feature1 = 'nCount_RNA', feature2 = 'percent.mito')
  p3 <- FeatureScatter(object = pbmc, feature1 = 'nCount_RNA', feature2 = 'percent.HB')
  pdf(paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI_mito_HB_relation.pdf'), width = 10)
  print(CombinePlots(plots = list(p1, p2, p3), ncol = 3, legend = 'none'))
  dev.off()
  f2g(paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI_mito_HB_relation.pdf'), paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI_mito_HB_relation.png'))
  
  #99% mitochondrion
  percent.mito_99 <- round(quantile(pbmc@meta.data$percent.mito, probs = 0.99), 2)
  #filter percentage mitochondrion
  if (percent.mito_99 > 30) {
    print('percent.mito_99>30!')
    percent.mito_99 <- 30
  } else {
    for (m in seq(0, 30, 5)) {
      nm <- m + 5
      if (percent.mito_99 < nm) {
        percent.mito_99 <- nm
        break
      }
    }
  }
  if (as.numeric(opt$high_percent.mito) == 30.0) {
    high_percent.mito <- percent.mito_99
  } else {
    high_percent.mito <- as.numeric(opt$high_percent.mito)
  }
  
  #Filter Cells
  pbmc <- subset(pbmc, subset = ((nFeature_RNA > as.numeric(opt$low_nGene)) & (nFeature_RNA < as.numeric(opt$high_nGene)) & (percent.mito < high_percent.mito) & (percent.HB < as.numeric(opt$high_percent.HB))))
} else {
  #Figure: base information
  #1
  pdf(paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI.pdf'), width = 10)
  print(VlnPlot(pbmc, features = c('nCount_RNA', 'nFeature_RNA'), ncol = 4))
  dev.off()
  f2g(paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI.pdf'), paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI.png'))
  #2
  pdf(paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI_relation.pdf'), width = 10)
  print(FeatureScatter(object = pbmc, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA'))
  dev.off()
  f2g(paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI_relation.pdf'), paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI_relation.png'))
  
  #filter genes
  pbmc <- subset(pbmc, subset = ((nFeature_RNA > as.numeric(opt$low_nGene)) & (nFeature_RNA < as.numeric(opt$high_nGene))))
}


#2.4 Doublet
prdate('2.4')

if (opt$doublet == "DoubletFinder") {
  #2.4.1 DoubletFinder
  pbmc <- NormalizeData(object = pbmc, normalization.method = 'LogNormalize', scale.factor = 10000)
  pbmc <- ScaleData(pbmc)
  pbmc <- FindVariableFeatures(object = pbmc, selection.method ='vst', mean.function = ExpMean, 
                               dispersion.function = LogVMR, 
                               mean.cutoff = c(as.numeric(opt$x_low_dispersion_cutoff), as.numeric(opt$x_high_dispersion_cutoff)), 
                               dispersion.cutoff = c(as.numeric(opt$y_high_dispersion_cutoff), Inf), nfeatures = as.numeric(opt$Nfeatures))
  pbmc <- RunPCA(pbmc)
  pbmc <- RunUMAP(pbmc, dims = 1:20)
  
  pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:20)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  
  bef_num = dim(pbmc@meta.data)[1]
  
  ## pK Identification (no ground-truth) ------------------------------------------------------------------------------------------
  sweep.res.list_kidney <- paramSweep(pbmc, PCs = 1:10, sct = FALSE) 
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn <- find.pK(sweep.stats_kidney)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- pbmc@meta.data$SCT_snn_res.0.5 
  homotypic.prop <- modelHomotypic(annotations)   
  nExp_poi <- round(0.075*nrow(pbmc@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  pbmc <- doubletFinder(pbmc, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE) 
  df_col <- grep("^DF.classifications", colnames(pbmc@meta.data), value = TRUE)
  pbmc <- subset(pbmc, cells = rownames(pbmc@meta.data)[pbmc@meta.data[[df_col]] == "Singlet"])
  
  aft_num = dim(pbmc@meta.data)[1]
  print(paste0("Before doublet remove, the cell number is: ", bef_num))
  print(paste0("After doublet remove, the cell number is: ", aft_num))
}

#2.5 Normalize Data
prdate('2.5')
pbmc <- NormalizeData(object = pbmc, normalization.method = 'LogNormalize', scale.factor = 10000)

#2.6 Variable Genes
prdate('2.6')
pbmc <- FindVariableFeatures(object = pbmc, selection.method ='vst', mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(as.numeric(opt$x_low_dispersion_cutoff), as.numeric(opt$x_high_dispersion_cutoff)), dispersion.cutoff = c(as.numeric(opt$y_high_dispersion_cutoff), Inf), nfeatures = as.numeric(opt$Nfeatures))
varFet <- VariableFeatures(object = pbmc)
mycell_number <- ncol(pbmc)
pbmc <- ScaleData(object = pbmc, features = varFet, vars.to.regress = c('nCount_RNA', 'percent.HB', 'percent.mito'), display.progress = FALSE)

tod<- GetAssayData(pbmc_raw, assay = "RNA", layer = "counts") 
toc<-GetAssayData(pbmc, assay = "RNA", layer = "counts") 
tod <- tod[rownames(toc),]

all <- toc
all <- CreateSeuratObject(all)
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)
all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all, resolution = 0.5)
all <- RunUMAP(all, dims = 1:20)
matx <- all@meta.data

# SoupX
scNoDrops = SoupChannel(tod, toc, calcSoupProfile = FALSE)
soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
scNoDrops = setSoupProfile(scNoDrops, soupProf)
sc<-scNoDrops
sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
sc = setContaminationFraction(sc, 0.2)
sc = autoEstCont(sc)
pbmc = sc

#2.8 Export QC information
prdate('2.7')
gene_number_filtered <- nrow(pbmc@assays$RNA)
cell_number_filtered <- ncol(pbmc@assays$RNA)
nGene_median_filtered <- as.integer(median(pbmc@meta.data$nFeature_RNA))

#2.8.1 Export filtered genes
gene_filtered <- data.frame(rownames(pbmc@assays$RNA))
colnames(gene_filtered) <- 'Gene_name'
gene_filtered <- left_join(gene_filtered, gene.names, by = 'Gene_name')
gene_filtered <- gene_filtered[, c(2,1)]
write.table(gene_filtered, file = paste0(opt$outdir, '/QC/', opt$id, '_gene_filtered.csv'), sep = ',', quote = F, row.names = F)

#2.8.2 Export setting information
if (opt$genome %in% c('GRCh38', 'mm10', 'Rattus_norvegicus')) {	
  filtered_information <- c(opt$min_cells, opt$low_nGene, opt$high_nGene, high_percent.mito, opt$high_percent.HB, gene_number, cell_number, nGene_median, gene_number_filtered, cell_number_filtered, nGene_median_filtered)
  filtered_information <- t(data.frame(filtered_information))
  colnames(filtered_information) <- c('min.cells', 'low.thresholds_nGene', 'high.thresholds_nGene', 'high.thresholds_percent.mito', 'high.thresholds_percent.HB', 'gene_number', 'cell_number', 'nGene_median', 'gene_number_filtered', 'cell_number_filtered', 'nGene_median_filtered')
} else {
  filtered_information <- c(opt$min_cells, opt$low_nGene, opt$high_nGene, gene_number, cell_number, nGene_median, gene_number_filtered, cell_number_filtered, nGene_median_filtered)
  filtered_information <- t(data.frame(filtered_information))
  colnames(filtered_information) <- c('min.cells', 'low.thresholds_nGene', 'high.thresholds_nGene', 'gene_number', 'cell_number', 'nGene_median', 'gene_number_filtered', 'cell_number_filtered', 'nGene_median_filtered')
}
write.table(filtered_information, file = paste0(opt$outdir, '/QC/', opt$id, '_filtered_information.csv'), sep = ',', quote = F, row.names = F)

#2.8.3 Export RDS
meta_data <- data.frame(pbmc@meta.data)
meta_data <- cbind(rownames(meta_data),meta_data)
colnames(meta_data)[1] <- 'Barcode'
write.table(meta_data[,1], file = paste0(opt$outdir, '/QC/', opt$id, '_barcode_filtered.csv'), sep = ',', row.names = F, quote = F, col.names = F)
write.table(meta_data, file = paste0(opt$outdir, '/QC/', opt$id, '_nGene_nUMI_filtered.csv'), sep = ',', quote = F, row.names = F)

#2.8.4 Figure for filtered information
pdf(paste0(opt$outdir, '/QC/', opt$id, '_after_filtered_nGene_nUMI_mito_HB.pdf'), width = 10)
print(VlnPlot(pbmc, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'percent.HB'), ncol = 4))
dev.off()
f2g(paste0(opt$outdir, '/QC/', opt$id, '_after_filtered_nGene_nUMI_mito_HB.pdf'), paste0(opt$outdir, '/QC/', opt$id, '_after_filtered__nGene_nUMI_mito_HB.png'))

######Step4 Save final RDS
prdate('4')
saveRDS(pbmc, file = paste0(opt$outdir, '/', opt$id, '.rds'))

######Step5 Seurat object to 3 files
Obj2files <- function(sce, dir){
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  ct <- GetAssayData(object = sce, assay = "RNA", layer = "counts")
  writeMM(obj = ct, file = paste0(dir, "/matrix.mtx"))
  write.table(data.frame(rownames(ct), rownames(ct)),
              file = paste0(dir, "/genes.tsv"),
              quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

  write.table(colnames(ct),
              file = paste0(dir, "/barcodes.tsv"),
              quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  message("Files saved to: ", dir)
}

prdate('5')
out = adjustCounts(sc)
DropletUtils:::write10xCounts(paste0(opt$outdir, "/", opt$id, "/FILES_matrix/"), out, version="3")

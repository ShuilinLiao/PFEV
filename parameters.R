sample_order_HC <- c("HC1", "HC2", "HC3")
sample_order_PD <- c("LTPD1", "LTPD2", "LTPD4", "STPD1", "STPD2", "STPD3", "STPD4", "STPD5", "STPD6")

group_order_HC <- c("HC", "HC", "HC")
group_order_PD <- c("LTPD", "LTPD", "LTPD", "STPD", "STPD", "STPD", "STPD", "STPD", "STPD")

cell_order_HC <- c("Mesoth", "Fibro", "Myofib", "Endo", "Mono", "Bcell", "Tcell")
cell_order_PD <- c("Perito", "Mono", "Tcell", "Bcell", "Unclassified")

cols_celltypes_HC <- c("#a586b7", "#fec582", "#38b9e6", "#cb8c93", "#f4efa6", "#89bbaa", "#bd2520")
cols_celltypes_PD <- c("#98dbce", "#fffcd9", "#eb8c3c", "#82cddd", "#9d9d9d")            

## == marker

def_markers_PD <- list(
  Perito =  c("MSLN", "TJP1", "PDPN", "ICAM1",  "UPK3B", "CDH1", "KRT19", "LRRN4", "CALB2", "CD200", 
              "WT1", "KRT5", "KRT7", "KRT8", "GPM6A", "PDGFRB"), 
  Mono = c("CD14", "FCN1", "C1QA", "C1QB", "C1QC", "MNDA", "CD68"), 
  Tcell = c("NKG7", "CCL5", "PRF1", "CD3D", "CD8A", "PDCD1", "GZMA", "PTPRC", "IGHG3", "MKI67", "TOP2A", "IL7R"),  
  Bcell = c("CD40", "CD19", "MS4A1", "IGHD")
)

def_markers_PD_8 <- list(
  PeritoS =  c("MSLN", "TJP1", "PDPN", "ICAM1",  "UPK3B", "CDH1", "KRT19", "LRRN4", "CALB2", "CD200", 
              "WT1", "KRT5", "KRT7", "KRT8", "GPM6A", "PDGFRB"), 
  PeritoL =  c("CDH2", "COL1A1", "COL3A1", "COL4A1"), # "VIM", 
  Mono = c("S100A4", "FCGR3A", "S100A8", "S100A9", "S100A12", "APOE",
           "CD14", "FCN1", "C1QA", "C1QB", "C1QC", "MNDA", "CD68"), 
  MonoDC = c("MRC1", "CD1A", "ITGAX", "CD1C"),
  Bcell = c("CD40", "CD19", "MS4A1", "IGHD"),
  NKT = c("NKG7", "CCL5", "PRF1", "CD3D", "CD8A", "PDCD1", "GZMA", "PTPRC"),
  Th = c("MKI67", "TOP2A"), # "IGHG3", 
  Nt = c("IL7R")
)


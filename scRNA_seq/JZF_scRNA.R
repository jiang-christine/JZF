library(sce.allrat)
library(combined_sce.allrat)
library(ggplot2)
library(clustree)
library(cowplot)
library(qs)
library(data.table)
library(dplyr)
library(FanCodeV1)
library(Seurat)
setwd("/home/jpc/JZF/scRNA_seq")


data_dir='raw_files' 
#fs=list.files('raw_files/','^GSM')
fs=list.files('./raw_files/')
fs=list.files('./raw_files/')
samples=fs

library(Matrix)

# 2. 创建空列表存储单个样本的sce.allrat对象
sce.allrat_list <- list()

# 3. 循环读取每个样本并创建sce.allrat对象
for (sample in samples) {
  # 构建样本路径
  sample_path <- file.path(data_dir, sample)
  
  # 读取10x Genomics数据
  sample_data <- Read10X(data.dir = sample_path, gene.column = 1) #
  
  # 创建sce.allrat对象，并添加样本标识
  sce.allrat_obj <- CreateSeuratObject(
    counts = sample_data,
    project = sample,  # 添加项目标识
    min.cells = 3,     # 过滤表达基因数少于3个细胞的基因
    min.features = 200 # 过滤检测到少于200个基因的细胞
  )
  
  # 添加样本信息到元数据
  sce.allrat_obj$sample <- sample
  
  # 添加到列表
  sce.allrat_list[[sample]] <- sce.allrat_obj
}

# 4. 合并所有样本
# 方法1：直接合并（适合样本间差异小）
combined_sce.allrat <- merge(
  x = sce.allrat_list[[1]],
  y = unlist(sce.allrat_list[-1]),
  add.cell.ids = samples  # 添加样本ID作为细胞前缀
)

sce.all <-basic_qc(combined_sce.allrat)
library(future)
library(harmony)
library(clustree)
sce.all <-run_harmony(sce.all)
plan(multisession, workers = 1)
DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.2",label = T, raster = FALSE)
DimPlot(sce.all , reduction = "umap", group.by = "sample",label = T, raster = FALSE)

ElbowPlot(sce.all, ndims = 20) 
qsave(sce.all,file="sce_all.qs",nthreads = 16)

sce.all <- subset(sce.all , subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent_mito < 10)

# library(DoubletFinder)
# sce.all <- MarkDoublets(sce.all, PCs=1:10, split.by=NULL) #因内存不够，没有运行成功

library(scDblFinder)
#sce.all <- scDblFinder(sce.all) 
sce <- as.SingleCellExperiment(sce.all)
sce <- scDblFinder(sce)
sce.all$scDblFinder.score <- sce$scDblFinder.score
sce.all$scDblFinder.class <- sce$scDblFinder.class

input_sce<- subset(sce.all,subset= scDblFinder.class=="singlet")


FeaturePlot(sce.all,features =c("scDblFinder.score"),reduction = "umap",raster = FALSE,pt.size = 0.01)


library(FanCodeV1)

library(clustree)
clustree(sce.all)

basic_markers <- list(
 "Immu" = c("PTPRC"),  # CD45

 # 淋巴细胞谱系
 "T" = c("CD3D", "CD3E", "CD3G"),
 "CD4_T" = c("CD4"),
 "CD8_T" = c("CD8A", "CD8B"),
 "B" = c("CD79A", "CD79B", "CD19", "MS4A1"),
 "Plasma" = c("JCHAIN", "IGKC"),
 "NK" = c("NKG7", "KLRG1", "KLRC1"),
 "ILC" = c("KIT", "KLRF1", "IL7R"),

 # 髓系细胞
 "Myeloid" = c("CD14", "LYZ", "S100A8"), #, "FCGR3A"
 "Mono" = c( "FCGR3A",  "MARCO"), #"CD68","CD14",
 "DC" = c("CD1C", "LAMP3"),
 "Macro" = c("CD68",  "MSR1", "CD163", "MRC1"),#"LYZ",
 "Neutro" = c("FCGR3B", "CXCR2", "SLC25A37", "G0S2", "CXCR1", "ADGRG3"),
 "Mast" = c("TPSAB1", "CPA3", "MS4A2"),
 # "cDC1" = c("CLEC9A", "XCR1", "BATF3", "IRF8", "CADM1", "THBD"),
 # "cDC2" = c("CD1C", "FCER1A", "CLEC10A", "IRF4", "ITGAX"),  # ITGAX = CD11c
 # 
 # # 浆细胞样树突状细胞
 # "pDC" = c("LILRA4", "GZMB", "IRF7", "TCF4", "CLEC4C", "IL3RA", "SPIB"),
 # 
 # # 成熟/活化树突状细胞
 # "mDC" = c("LAMP3", "CCR7", "CD40", "CD83", "CCL19", "CCL22", "IDO1"),
 # 
 # 基质细胞
 "Fibro" = c("PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "DCN", "FAP"),
 "Myofibro" = c("ACTA2", "TAGLN", "MYH11"),
 "Endo" = c("PECAM1", "VWF", "CDH5", "ENG", "TIE1", "CLDN5"),
 "Peri" = c("RGS5",  "CSPG4", "MCAM"), #"PDGFRB",

"Epi" = c("EPCAM", "KRT18", "KRT19"))

# 注释掉的标记
# "Cell Cycle" = c("TOP2A", "MKI67"),
# "Platelets" = c("PF4", "PPBP")


DotPlot(input_sce, features = basic_markers,cols = c("grey","red"), cluster.idents = F) + 
  RotatedAxis() + 
  theme( panel.border = element_rect(color="black"), panel.spacing = unit(1,"mm"), strip.text = element_text(margin=margin(b=3, unit="mm")),strip.placement ='outlet',axis.line = element_blank(), ) + 
  labs(x="", y="")



celltype=data.frame(ClusterID=0:28,
                    celltype= 0:28) 

celltype[celltype$ClusterID %in% c(0,3,7,17,28),2]= "T cells"
celltype[celltype$ClusterID %in% c(5),2]= "B cells"
celltype[celltype$ClusterID %in% c(25),2]='Mast cell'
celltype[celltype$ClusterID %in% c(21),2]='Plasma'
celltype[celltype$ClusterID %in% c(1,2,8,20),2]='NK cell'
celltype[celltype$ClusterID %in% c(4,14,19),2]='Myeloid_unknown'
celltype[celltype$ClusterID %in% c(11),2]= "DC"
celltype[celltype$ClusterID %in% c(12),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(10,13,14,27),2]='Mono/Macro'
celltype[celltype$ClusterID %in% c(22),2]= "Epi"
celltype[celltype$ClusterID %in% c(6,15,16,24,26,9),2]= "Endo"
celltype[celltype$ClusterID %in% c(23),2]='Myofibroblast'
celltype[celltype$ClusterID %in% c(18),2]='Unknown'

input_sce@meta.data$Celltype = "NA"
for(i in 1:nrow(celltype)){
  input_sce@meta.data[which(input_sce @meta.data$RNA_snn_res.1 == celltype$ClusterID[i]),'Celltype'] <- celltype$celltype[i]}
table(input_sce@meta.data$Celltype)

DimPlot(input_sce, reduction = "umap", group.by = "Celltype",label = T, raster = FALSE)

DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.1",label = T, raster = FALSE)

#############----------------from python
cellinfo <- my_load_tibble("umap_df.csv")


cellinfo$cell <- paste(cellinfo$sample,cellinfo$X,sep="_")
cellinfo$cellID <- sub("-.*$", "", cellinfo$cell)
length(intersect(cellinfo$cellID,rownames(metadat_11)))

setdiff(cellinfo$cellID,intersect(cellinfo$cellID,rownames(metadat_11)))
rownames(metadat_11)[duplicated(rownames(metadat_11))][1:10]

metadat_11 <- sce.all@meta.data
#metadat_11$cellid <- gsub("_", "", rownames(metadat_11))
metadat_11$cellid <-rownames(metadat_11)

length(intersect(cellinfo$cell,metadat_11$cellid))
co_cell <- intersect(cellinfo$cell,rownames(metadat_11))

sce11 <- subset(sce.all,cells =co_cell)
metadat_12 <- sce11@meta.data
metadat_12<-as_tibble(metadat_12,rownames = "cell")
cellinfo_12 <- cellinfo %>% filter(cell %in% co_cell)
length(intersect(cellinfo$cell,metadat_12$cell))
metadat_13<- inner_join(metadat_12,cellinfo_12[,c(23,36:39)],by="cell")

sce11 <-AddMetaData(sce11,metadata = metadat_13)


# 2. 提取UMAP坐标矩阵（移除cellID列）
umap_matrix <- as.matrix(metadat_13[, c("UMAP1", "UMAP2")])
rownames(umap_matrix) <- metadat_13$cell

# 3. 创建DimReduc对象
umap_reduction <- CreateDimReducObject(
  embeddings = umap_matrix,
  key = "UMAP_",  # 降维键名前缀
  assay = "RNA"   # 指定关联的assay（根据实际修改）
)

sce11@reductions$umap_python <- umap_reduction
DimPlot(sce11, group.by= "cell_type_lvl1",reduction = "umap_python", label = T, raster = FALSE)+
  DimPlot(sce11, group.by= "cell_type_lvl1",reduction = "umap", label = T, raster = FALSE)+ NoLegend()








sweep.res <- paramSweep_v3(sce.all, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)  # 选择最优pK值

# 4. 计算双细胞比例 (每1000细胞约7.5%)
nExp_poi <- round(0.075 * nrow(sce.all@meta.data))  # 调整比例

# 5. 运行DoubletFinder
sce.all <- doubletFinder_v3(
  seu = sce.all,
  PCs = 1:10,
  pN = 0.25,                   # 默认值
  pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])),  # 自动选择
  nExp = nExp_poi,
  reuse.pANN = FALSE
)

# 6. 查看结果 (列名可能不同)
doublet_col <- grep("DF.classifications", colnames(sce.all@meta.data), value = TRUE)
table(sce.all@meta.data[[doublet_col]])

# 7. 移除双细胞
singlets <- subset(sce.all, cells = colnames(sce.all)[sce.all@meta.data[[doublet_col]] == "Singlet"])



basic.markers <- list(
  # 淋巴细胞谱系
  "T" = c("CD3D", "CD3E", "CD8A", "CD4", "IL7R", "TRAC"),
  "B" = c("CD79A", "MS4A1", "CD19", "CD22"),
  "NK" = c("NCAM1", "NKG7", "GNLY"),
  
  # 髓系细胞
  "Myeloid" = c("CD14", "LYZ", "S100A8", "FCGR3A"),
  
  # 上皮细胞
  "Epithelial" = c("EPCAM", "KRT8", "KRT18", "CDH1"),
  
  # 基质细胞亚群
  "Fibroblasts" = c("COL1A1", "DCN", "LUM", "PDGFRA"),
  "Endothelial" = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  "Pericytes" = c("ACTA2", "MYH11", "DES", "TAGLN")
)

sce.all@meta.data$seurat_clusters <- sce.all@meta.data$RNA_snn_res.0.2
Idents(sce.all) <- "seurat_clusters"

DotPlot(sce.all, features = basic.markers,cols = c("grey","red"), cluster.idents = F) + 
  RotatedAxis() + 
  theme( panel.border = element_rect(color="black"), panel.spacing = unit(1,"mm"), strip.text = element_text(margin=margin(b=3, unit="mm")),strip.placement ='outlet',axis.line = element_blank(), ) + 
  labs(x="", y="")

library(qs)
scobj <- qread("scobj.qs",nthreads = 16)
 


sweep.res <- paramSweep(sce.all, PCs = 1:10, sct = FALSE,num.cores = 1)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# 查看最佳 pK 值（用于后面）
best.pK <- bcmvn[which.max(bcmvn$BCmetric), "pK"]
# 比例通常 5%-10%，可以自己估计或多试几个值
homotypic.prop <- modelHomotypic(sce.all$sce.allrat_clusters)  # 计算同型双细胞比例
nExp_poi <- round(0.075 * ncol(sce.all))  # 设定期望双细胞数量为细胞总数的 7.5%
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))  # 过滤掉同型双细胞
sce.all <- doubletFinder_v3(
  sce.all,
  PCs = 1:10,
  pN = 0.25,
  pK = best.pK,
  nExp = nExp_poi.adj,
  reuse.pANN = FALSE,
  sct = FALSE
)

sce.all

qsave(input_sce,file="data/input_sce.qs",nthreads = 16) #用R分群的数据，去双细胞



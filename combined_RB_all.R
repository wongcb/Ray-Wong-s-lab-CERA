########Spatial transcriptomic profiling of human retinoblastoma############
########created by Luozixian Wang and Raymond Wong##########################
########CERA, UNIMELB, 05/02/2024###########################################
##load the packages
library(AUCell)
library(BiocManager)
library(bmixture)
library(CARD)
library(CellChat)
library(clusterProfiler)
library(cowplot)
library(DDRTree)
library(DESeq2)
library(devtools)
library(doParallel)
library(dplyr)
library(EnhancedVolcano)
library(enrichplot)
library(foreach)
library(future)
library(garnett)
library(ggalluvial)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(harmony)
library(infercnv)
library(kableExtra)
library(KernSmooth)
library(Matrix)
library(MuSiC)
library(NMF)
library(NNLM)
library(patchwork)
library(pathview)
library(promises)
library(RColorBrewer)
library(rdist)
library(remotes)
library(reticulate)
library(scCustomize)
library(SCDC)
library(SCENIC)
library(SCopeLoomR)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratObject)
library(SeuratWrappers)
library(SingleCellExperiment)
library(spdep)
library(STutility)
library(swne)
library(tibble)
library(tidyverse)
library(zellkonverter)
options(stringsAsFactors = FALSE)

#####QC and processing#####
#####RB####################
## filtered spots (under tissue)
infoTable <- read.csv('infoTable.csv')
se <- InputFromTable(infoTable)

##mito + ribo content
# Collect all genes coded on the mitochondrial genome
mt.genes <- grep(pattern = "^MT-", x = rownames(se), value = TRUE)
se$percent.mito <- (Matrix::colSums(se@assays$RNA@counts[mt.genes, ])/Matrix::colSums(se@assays$RNA@counts))*100

# Collect all genes coding for ribosomal proteins
rp.genes <- grep(pattern = "^RPL|^RPS", x = rownames(se), value = TRUE)
se$percent.ribo <- (Matrix::colSums(se@assays$RNA@counts[rp.genes, ])/Matrix::colSums(se@assays$RNA@counts))*100

# QC for spots with unique gene and mitochondria percentage
se.QC <- SubsetSTData(se, expression = nFeature_RNA > 20 & percent.mito < 30)
cat("Spots removed: ", ncol(se) - ncol(se.QC), "\n")

UMIplot.se.QC <- ST.FeaturePlot(se.QC, features = "nFeature_RNA", cols = c("lightgray", "mistyrose", "red", "dark red", "black"), pt.size = 2, ncol = 2)
UMIplot.se <- ST.FeaturePlot(se, features = "nFeature_RNA", cols = c("lightgray", "mistyrose", "red", "dark red", "black"), pt.size = 2, ncol = 2)
custom_theme <- theme(legend.position = c(0.45, 0.8), # Move color legend to top
                      legend.direction = "horizontal", # Flip legend
                      legend.text = element_text(angle = 30, hjust = 1), # rotate legend axis text
                      strip.text = element_blank(), # remove strip text
                      plot.title = element_blank(), # remove plot title
                      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) # remove plot margins
UMIplot.se.QC & theme_bw()

## loading images
se.QC <- LoadImages(se.QC, time.resolve = FALSE, verbose = TRUE)

## manual annotation 
se.QC <- ManualAnnotation(se.QC)

## QC check counts metrics 
QC.annotated <- VlnPlot(se.QC, features = c("nFeature_RNA", "nCount_RNA"), group.by = "labels")
saveRDS(se.QC, 'se.QC.RDS')

##subset out the unnannotated
Idents(se.QC) <- 'labels'
se.QC.anno<- SubsetSTData(se.QC, idents = c('RB 1', 'RB 2', 'RB 3', 'RB 4'))

## SCTransform
# Add a section column to your meta.data
se.QC.anno$section <- paste0("section_", GetStaffli(se.QC.anno)[[, "sample", drop = T]])
table(se.QC.anno$section)
se.QC.anno <- SCTransform(se.QC.anno, vars.to.regress = "section")
FeatureOverlay(se.QC.anno, features = "labels", ncols = 2)

# simple workflow for Dimensionality reduction (PCA), UMAP embedding, Clustering
se.QC.anno <- se.QC.anno %>% 
  RunPCA() %>%
  RunUMAP(reduction = "pca", dims = 1:30)
se.QC.anno <- se.QC.anno %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(reduction = "pca", dims = 1:30)
se.QC.anno$seurat_clusters_pca <- se.QC.anno$seurat_clusters

#Clustering  
se.QC.anno.UMAPp1 <- DimPlot(se.QC.anno, group.by = "labels", reduction = "umap")
se.QC.anno.UMAPp2 <- DimPlot(se.QC.anno, group.by = "seurat_clusters_pca", label = TRUE, label.size = 8, reduction = "umap")
se.QC.anno.UMAPp1 - se.QC.anno.UMAPp2
se.QC.anno.clusterSTp1 <- ST.FeaturePlot(se.QC.anno, features = "seurat_clusters_pca", indices = 1, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
se.QC.anno.clusterSTp2 <- ST.FeaturePlot(se.QC.anno, features = "seurat_clusters_pca", indices = 2, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
se.QC.anno.clusterSTp3 <- ST.FeaturePlot(se.QC.anno, features = "seurat_clusters_pca", indices = 3, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
se.QC.anno.clusterSTp4 <- ST.FeaturePlot(se.QC.anno, features = "seurat_clusters_pca", indices = 4, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
cowplot::plot_grid(se.QC.anno.clusterSTp1, se.QC.anno.clusterSTp2, se.QC.anno.clusterSTp3, se.QC.anno.clusterSTp4, ncol = 4)

# use harmony for integration across sections
se.QC.anno <- RunHarmony(se.QC.anno, group.by.vars = "labels", reduction = "pca", dims.use = 1:30, assay.use = "SCT", verbose = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters()
se.QC.anno$seurat_clusters_harmony <- se.QC.anno$seurat_clusters

## cell cycle analysis 
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
se.QC.anno <- CellCycleScoring(se.QC.anno, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
## stored in 'Phase' Idents

#####Retina####################
## filtered spots (under tissue)
infoTable <- read.csv('infoTable.csv')
se <- InputFromTable(infoTable)

##mito + ribo content
# Collect all genes coded on the mitochondrial genome
mt.genes <- grep(pattern = "^MT-", x = rownames(se), value = TRUE)
se$percent.mito <- (Matrix::colSums(se@assays$RNA@counts[mt.genes, ])/Matrix::colSums(se@assays$RNA@counts))*100

# Collect all genes coding for ribosomal proteins
rp.genes <- grep(pattern = "^RPL|^RPS", x = rownames(se), value = TRUE)
se$percent.ribo <- (Matrix::colSums(se@assays$RNA@counts[rp.genes, ])/Matrix::colSums(se@assays$RNA@counts))*100
saveRDS(se, 'pooled.retina.se.RDS')

# QC for spots with unique gene and mitochondria percentage
se.QC <- SubsetSTData(se, expression = nFeature_RNA > 200 & percent.mito < 30)
cat("Spots removed: ", ncol(se) - ncol(se.QC), "\n")

UMIplot.se.QC <- ST.FeaturePlot(se.QC, features = "nFeature_RNA", cols = c("lightgray", "mistyrose", "red", "dark red", "black"), pt.size = 2, ncol = 2)
UMIplot.se <- ST.FeaturePlot(se, features = "nFeature_RNA", cols = c("lightgray", "mistyrose", "red", "dark red", "black"), pt.size = 2, ncol = 2)
custom_theme <- theme(legend.position = c(0.45, 0.8), # Move color legend to top
                      legend.direction = "horizontal", # Flip legend
                      legend.text = element_text(angle = 30, hjust = 1), # rotate legend axis text
                      strip.text = element_blank(), # remove strip text
                      plot.title = element_blank(), # remove plot title
                      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) # remove plot margins
UMIplot.se + theme_bw() +labs(title='nFeature_RNA in se') | UMIplot.se.QC + theme_bw() +labs(title='nFeature_RNA in se.QC (>200 genes + <30 mito%)')
dev.off()

## loading images
se.QC <- LoadImages(se.QC, time.resolve = FALSE, verbose = TRUE)

## manual annotation 'off tissue' = spots with no tissue/crap
se.QC <- ManualAnnotation(se.QC)
Idents(se.QC) <- 'labels'
saveRDS(se.QC, 'pooled.retina.se.QC.RDS')

##subset out the 'off tissue' spots
se.QC.anno<- SubsetSTData(se.QC, idents = 'Default')

## manual label eye compartments
# Others #folded, unknown, discarded
# Retina
# Choroid
# Sclera
# Optic nerve
se.QC.anno <- ManualAnnotation(se.QC.anno)
Idents(se.QC.anno) <- 'labels'
saveRDS(se.QC.anno, 'pooled.retina.se.QC.anno.all.RDS')

## subset out the 'Others' spots => folded, unknown, discarded
se.QC.anno<- SubsetSTData(se.QC.anno, idents = c('Choroid', 'Retina', 'Sclera', 'Optic nerve'))

## SCTransform
# Add a section column to your meta.data
se.QC.anno$section <- paste0("section_", GetStaffli(se.QC.anno)[[, "sample", drop = T]])
table(se.QC.anno$section)
se.QC.anno <- SCTransform(se.QC.anno, vars.to.regress = "section")

###Overlay of annotation in H+E
FeatureOverlay(se.QC.anno, features = "labels", ncols = 2, sample = 1:4)

# simple workflow for Dimensionality reduction (PCA), UMAP embedding, Clustering
se.QC.anno <- se.QC.anno %>% 
  RunPCA() %>%
  RunUMAP(reduction = "pca", dims = 1:30)
se.QC.anno <- se.QC.anno %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(reduction = "pca", dims = 1:30)
se.QC.anno$seurat_clusters_pca <- se.QC.anno$seurat_clusters

#Clustering  
se.QC.anno.UMAPp1 <- DimPlot(se.QC.anno, group.by = "labels", reduction = "umap")
se.QC.anno.UMAPp2 <- DimPlot(se.QC.anno, group.by = "seurat_clusters_pca", label = TRUE, label.size = 8, reduction = "umap")
se.QC.anno.clusterSTp1 <- ST.FeaturePlot(se.QC.anno, features = "seurat_clusters_pca", indices = 1, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
se.QC.anno.clusterSTp2 <- ST.FeaturePlot(se.QC.anno, features = "seurat_clusters_pca", indices = 2, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
se.QC.anno.clusterSTp3 <- ST.FeaturePlot(se.QC.anno, features = "seurat_clusters_pca", indices = 3, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
se.QC.anno.clusterSTp4 <- ST.FeaturePlot(se.QC.anno, features = "seurat_clusters_pca", indices = 4, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())

## cell cycle analysis 
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
se.QC.anno <- CellCycleScoring(se.QC.anno, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
## stored in 'Phase' Idents

######SP <- se.QC.anno from RB##########
######STeye <- se.QC.anno from Retina###
##figure 1B
STeye <- MaskImages(object = STeye)
ImagePlot(STeye, method = "raster", type = "masked", ncols = 2)
FeatureOverlay(STeye, features = "labels", ncols = 2)
SP <- MaskImages(object = SP)
ImagePlot(SP, method = "raster", type = "masked", ncols = 2)
ST.FeaturePlot(object = SP, features = "seurat_clusters_harmony", pt.size = 3) + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 8), legend.position = "bottom")

##figure 1C
DimPlot(SP, reduction = "umap.harmony", label = TRUE, pt.size = 2.0, label.size = 12)

##figure 1E-H
##3D stack
#By runnig Create3DStack, we can create a z-stack of “2D point patterns” which we’ll use to interpolate expression values over and visualzie expression in 2D space.
SP <- Create3DStack(SP)
# 3d Plot RB genes
FeaturePlot3D(SP, features = "UBE2C", pt.size = 0.6, max.cutoff = 5)
FeaturePlot3D(SP, features = "RB1", pt.size = 0.6, max.cutoff = 5)
FeaturePlot3D(SP, features = "MDM2", pt.size = 0.6, max.cutoff = 5)
FeaturePlot3D(SP, features = "MYCN", pt.size = 0.6, max.cutoff = 5)

##figure 2A, 3B-E, module score
#modulescore
#Retinal progenitor cell type
RPCscore <- Seurat::AddModuleScore(object = SP, features = c("SOX2", "HES1", "MKI67", "HES5", "FZD5", "PAX6"), name = "RPC_enriched")
colnames(RPCscore@meta.data)[21] <- "Retinal progenitor cell Score"
pRPC <- VlnPlot(RPCscore,features = "Retinal progenitor cell Score", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
pRPC
#Pigment epithelium cell type
PECscore <- Seurat::AddModuleScore(object = SP, features = c("SERPINF1", "MITF", "BEST1", "TTR"), name = "PEC_enriched")
colnames(PECscore@meta.data)[19] <- "Pigment epithelium cell Score"
pPEC <- VlnPlot(PECscore,features = "Pigment epithelium cell Score", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
pPEC
#Rod cell type
RODscore <- Seurat::AddModuleScore(object = SP, features = c("RHO", "PDE6A", "CNGA1", "NRL", "GNAT1", "GNB1", "SAG", "ELOVL4", "PDE6B", "GNGT1"), name = "ROD_enriched")
colnames(RODscore@meta.data)[25] <- "Rod Cell Score"
pROD <- VlnPlot(RODscore,features = "Rod Cell Score", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
pROD
#Cone cell type
CONEscore <- Seurat::AddModuleScore(object = SP, features = c("ARR3", "GNGT2", "PDE6H", "GUCA1C", "GNAT2", "RXRG", "THRB", "PDC", "GNB3", "CRX"), name = "CONE_enriched")
colnames(CONEscore@meta.data)[25] <- "Cone Cell Score"
pCONE <- VlnPlot(CONEscore,features = "Cone Cell Score", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
pCONE
#Bipolar cell type
BCscore <- Seurat::AddModuleScore(object = SP, features = c("GRM6", "VSX2"), name = "BC_enriched")
colnames(BCscore@meta.data)[16] <- "Bipolar cell Score"
pBC <- VlnPlot(BCscore,features = "Bipolar cell Score", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#Horizontal cell type
HCscore <- Seurat::AddModuleScore(object = SP, features = c("ONECUT1", "ONECUT2", "TFAP2B"), name = "HC_enriched")
colnames(HCscore@meta.data)[18] <- "Horizontal Cell Score"
pHC <- VlnPlot(HCscore,features = "Horizontal Cell Score", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#Amacrine cell type
ACscore <- Seurat::AddModuleScore(object = SP, features = c("GAD1", "CALB1", "NRXN2", "TFAP2A", "PROX1", "GAD2"), name = "AC_enriched")
colnames(ACscore@meta.data)[21] <- "Amacrine cell Score"
pAC <- VlnPlot(ACscore,features = "Amacrine cell Score", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#Retinal ganglion cell type
RGCscore <- Seurat::AddModuleScore(object = SP, features = c("POU4F2", "GAP43", "NEFL", "SNCG", "ATOH7", "EBF3", "THY1", "NRN1"), name = "RGC_enriched")
colnames(RGCscore@meta.data)[23] <- "Retinal ganglion cell Score"
pRGC <- VlnPlot(RGCscore,features = "Retinal ganglion cell Score", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#Cone precursor cell type
CPscore <- Seurat::AddModuleScore(object = SP, features = c("CRX", "RXRG", "THRB"), name = "CP_enriched")
colnames(CPscore@meta.data)[18] <- "Cone Precursor Cell Score"
pCP <- VlnPlot(CPscore,features = "Cone Precursor Cell Score", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#glial cells
Gscore <- Seurat::AddModuleScore(object = SP, features = c("CD68", "HLA-DPA1", "HLA-DPB1", "CLU"), name = "G_enriched")
colnames(Gscore@meta.data)[19] <- "Glial Cell Score"
pG <- VlnPlot(Gscore,features = "Glial Cell Score", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#cancer-associated fibroblasts (CAF)
CAFscore <- Seurat::AddModuleScore(object = SP, features = c("ACTA2", "VIM", "FGF9"), name = "CAF_enriched")
colnames(CAFscore@meta.data)[18] <- "Cancer-Associated Fibroblast Score"
pCAF <- VlnPlot(CAFscore,features = "Cancer-Associated Fibroblast Score", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#proliferative cells
proscore <- Seurat::AddModuleScore(object = SP, features = c("MKI67", "TOP2A", "KIF14"), name = "pro_enriched")
colnames(proscore@meta.data)[18] <- "Highly-Proliferated CP Score"
ppro <- VlnPlot(proscore,features = "Highly-Proliferated CP Score", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))

plot_grid(pRPC, pPEC, pROD, pCONE, pHC, pAC, pRGC, pBC, labels = "AUTO", nrow = 2)
pCP
pG
pCAF
ppro

##figure 3A
RidgePlot(SP, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

##figure 3C
pbar <- SP@meta.data %>%
  group_by(seurat_clusters_harmony,Phase) %>%
  dplyr::count() %>%
  group_by(seurat_clusters_harmony) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters_harmony,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")
pbar

##figure 4A
de.markers <- FindAllMarkers(SP, only.pos = TRUE)
top10 <- de.markers %>%
  dplyr::filter(p_val_adj < 0.01) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(wt = -p_val_adj, n = 10)
SP.DEheatmap<- DoHeatmap(SP, features = top10$gene)
SP.DEheatmap

#supplementary figure 2
RB_check <- c("RB1", "MDM2", "E2F1", "MYCN", "UBE2C", "INK4A", "ARF", "MDM4")
DotPlot(object = SP, features = RB_check, cols = c("cornsilk", "red2"), dot.scale = 14, 
        group.by = "seurat_clusters_harmony", cluster.idents = F, scale.by = "size") + 
  theme(legend.position = "right", 
        axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18), legend.title = element_text(size = 18))

#supplementary figure 3A
####intergrate the RB and Eye dataset
hRetina <- subset(STeye, idents = "Retina")
Idents(STeye)
hChoroid <- subset(STeye, idents = "Choroid")
hSclera <- subset(STeye, idents = "Sclera")
hOptic_nerve <- subset(STeye, idents = "Optic nerve")
hRetina$type <- "Healthy_Retina"
hChoroid$type <- "Choroid"
hSclera$type <- "Sclera"
hOptic_nerve$type <- "Optic_nerve"
SP$type <- "RB"
#integration between the splitted datasets
split_retina.anchors <- FindIntegrationAnchors(object.list = list(hRetina, hChoroid, hSclera, hOptic_nerve, SP), dims = 1:20)
split_retina.combined <- IntegrateData(anchorset = split_retina.anchors, dims = 1:20)
DefaultAssay(split_retina.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
split_retina.combined <- ScaleData(split_retina.combined, verbose = FALSE)
split_retina.combined <- RunPCA(split_retina.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
split_retina.combined <- RunUMAP(split_retina.combined, reduction = "pca", dims = 1:20)
split_retina.combined <- FindNeighbors(split_retina.combined, reduction = "pca", dims = 1:20)
split_retina.combined <- FindClusters(split_retina.combined, resolution = 0.5)
# Visualization
p1 <- DimPlot(split_retina.combined, reduction = "umap", group.by = "type", pt.size = 1, label.size = 3.5)
p1

##supplementary figure 3B, C
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

split_retina.combined_cycle <- CellCycleScoring(split_retina.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="SCT")
DimPlot(split_retina.combined_cycle, reduction = "umap", label = TRUE, pt.size = 1, label.size = 3.5)
head(split_retina.combined_cycle[[]])
# Visualize the distribution of cell cycle markers across
RidgePlot(split_retina.combined_cycle, features = c("PCNA", "MCM6", "TOP2A", "MKI67"), ncol = 2, group.by = "type")
#barchart
pbar_cycle_split <- split_retina.combined_cycle@meta.data %>%
  group_by(type,Phase) %>%
  dplyr::count() %>%
  group_by(type) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=type,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per sample type")
pbar_cycle_split


##supplementary figure 3D
#DE analysis
DefaultAssay(split_retina.combined) <- "RNA"
head(split_retina.combined[[]])
levels(split_retina.combined)
Idents(split_retina.combined) <- split_retina.combined@meta.data$type
#FindAllMarkers
Eye_RB_DEG <- FindAllMarkers(split_retina.combined, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- Eye_RB_DEG %>%
  dplyr::filter(p_val_adj < 0.01) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(wt = -p_val_adj, n = 10)
Eye_RB.DEheatmap <- DoHeatmap(split_retina.combined, slot = "counts", features = top10$gene)
Eye_RB.DEheatmap


####inferCNV######
counts_matrix = GetAssayData(SP, slot="counts")
#switch active identity to cell type
Idents(STeye) <- "labels2"
DefaultAssay(STeye) <- "RNA"
STeye@meta.data$tissue_type <- STeye@active.ident
DimPlot(STeye, reduction = "umap", label = TRUE, pt.size = 1.3, label.size = 4.5)
#subset the retina part
STretina <- subset(STeye, idents = "Retina", invert = FALSE)
## merge
RB_STretina.combined <- merge(SP, y = STretina, add.cell.ids = c("SP", "STeye"), project = "Merge")
RB_STretina.combined@meta.data[2982:2987,]
RB_STretina.combined@meta.data$tissuetype <-c(RB_STretina.combined@meta.data[1:2982,15],RB_STretina.combined@meta.data[2983:3426,18])
RB_STretina.combined@meta.data[2977:2987,]
saveRDS(RB_STretina.combined,"RB_STretina.combined.rds")
#prepare the input data for infercnv
dfcount = GetAssayData(RB_STretina.combined, slot="counts")
dfids <- rownames(RB_STretina.combined@meta.data)
tissue_types <- RB_STretina.combined@meta.data$tissuetype
tissue_info <- data.frame(CellID = dfids, CellType = tissue_types, stringsAsFactors = FALSE)
write.table(tissue_info, file = "RB_STretina_tissue_info.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#estimate CNV
#create the infercnv object for cell type compare
infercnv_obj_RB_retina = CreateInfercnvObject(raw_counts_matrix=dfcount,
                                              annotations_file="RB_STretina_tissue_info.txt",
                                              delim="\t",
                                              gene_order_file="hg38_gencode_v27.txt",
                                              ref_group_names=c("Retina"))
# perform infercnv operations to reveal cnv signal
infercnv_obj_RB_retina = infercnv::run(infercnv_obj_RB_retina,
                                       cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                       out_dir="RB_STretina_output_dir",  # dir is auto-created for storing outputs
                                       cluster_by_groups=T,   # cluster
                                       denoise=T,
                                       HMM=T, 
                                       num_threads = 16,
                                       write_expr_matrix = T)

#figure 1D in the "RB_STretina_output_dir"

#CNV scores
#The CNV score was defined as the mean squares of CNV values across the genome
cnv_score_table_v3 = data.table::fread("RB_STretina_output_dir/infercnv.observations.txt", 
                                       data.table = F) %>% 
  column_to_rownames(var = 'V1')
library(scales)
cnvScore <- function(data){
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    rescale(to=c(-1, 1)) %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}
cnv_score_v3 <- cnvScore(cnv_score_table_v3)
cellAnnota <- subset(RB_STretina.combined@meta.data, select = c('tissuetype'))
cnv_score_v3 <- cbind(cnv_score_v3, cellAnnota[row.names(cnv_score_v3),])
names(cnv_score_v3) <- c('cnv_score_v3', 'celltype')
color <- ggsci::pal_aaas()(10)

##supplementary figure 4
ggplot(cnv_score_v3, aes(celltype, cnv_score_v3, color = celltype)) +
  geom_boxplot() +
  scale_color_manual(values = color) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "NA") +
  labs(x = '', y = 'CNV Scores', title = '') +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1)) +
  stat_compare_means()

########Conditional autoregressive-based deconvolution (CARD)#####
#spatial
stRB <- SP
#single cell: Wu, C. et al. Single-cell characterization of malignant phenotypes and microenvironment alteration in retinoblastoma. Cell Death Dis. 13, 438 (2022).
scRB <- readRDS("RB_data.RDS")
#generate sp_count and sp_location
sp_count <- stRB@assays$RNA@counts
RB4spatial <- SubsetSTData(stRB, expression = labels %in% "RB 4")
ImagePlot(stRB, method = "raster")
ImagePlot(RB4spatial, method = "raster")
st.RB4 <- GetStaffli(RB4spatial)
st.RB4
RB4_coordinate <- st.RB4@meta.data
RB4sp_location <- data.frame(RB4_coordinate)
RB4sp_location$adj_x <- NULL
RB4sp_location$adj_y <- NULL
RB4sp_location$pixel_x <- NULL
RB4sp_location$pixel_y <- NULL
RB4sp_location$original_x <- NULL
RB4sp_location$original_y <- NULL
RB4sp_location$warped_x <- NULL
RB4sp_location$warped_y <- NULL
RB4sp_location$sample <- NULL
colnames(RB4sp_location) <- c("x", "y")
plot(RB4sp_location, method = "raster")
RB4sp_count <- sp_count[, colnames(sp_count) %in% rownames(RB4sp_location)]
RB3spatial <- SubsetSTData(stRB, expression = labels %in% "RB 3")
ImagePlot(stRB, method = "raster")
ImagePlot(RB3spatial, method = "raster")
st.RB3 <- GetStaffli(RB3spatial)
st.RB3
RB3_coordinate <- st.RB3@meta.data
RB3sp_location <- data.frame(RB3_coordinate)
RB3sp_location$adj_x <- NULL
RB3sp_location$adj_y <- NULL
RB3sp_location$pixel_x <- NULL
RB3sp_location$pixel_y <- NULL
RB3sp_location$original_x <- NULL
RB3sp_location$original_y <- NULL
RB3sp_location$warped_x <- NULL
RB3sp_location$warped_y <- NULL
RB3sp_location$sample <- NULL
colnames(RB3sp_location) <- c("x", "y")
plot(RB3sp_location, method = "raster")
RB3sp_count <- sp_count[, colnames(sp_count) %in% rownames(RB3sp_location)]
RB2spatial <- SubsetSTData(stRB, expression = labels %in% "RB 2")
ImagePlot(stRB, method = "raster")
ImagePlot(RB2spatial, method = "raster")
st.RB2 <- GetStaffli(RB2spatial)
st.RB2
RB2_coordinate <- st.RB2@meta.data
RB2sp_location <- data.frame(RB2_coordinate)
RB2sp_location$adj_x <- NULL
RB2sp_location$adj_y <- NULL
RB2sp_location$pixel_x <- NULL
RB2sp_location$pixel_y <- NULL
RB2sp_location$original_x <- NULL
RB2sp_location$original_y <- NULL
RB2sp_location$warped_x <- NULL
RB2sp_location$warped_y <- NULL
RB2sp_location$sample <- NULL
colnames(RB2sp_location) <- c("x", "y")
plot(RB2sp_location, method = "raster")
RB2sp_count <- sp_count[, colnames(sp_count) %in% rownames(RB2sp_location)]
RB1spatial <- SubsetSTData(stRB, expression = labels %in% "RB 1")
ImagePlot(stRB, method = "raster")
ImagePlot(RB1spatial, method = "raster")
st.RB1 <- GetStaffli(RB1spatial)
st.RB1
RB1_coordinate <- st.RB1@meta.data
RB1sp_location <- data.frame(RB1_coordinate)
RB1sp_location$adj_x <- NULL
RB1sp_location$adj_y <- NULL
RB1sp_location$pixel_x <- NULL
RB1sp_location$pixel_y <- NULL
RB1sp_location$original_x <- NULL
RB1sp_location$original_y <- NULL
RB1sp_location$warped_x <- NULL
RB1sp_location$warped_y <- NULL
RB1sp_location$sample <- NULL
colnames(RB1sp_location) <- c("x", "y")
plot(RB1sp_location, method = "raster")
RB1sp_count <- sp_count[, colnames(sp_count) %in% rownames(RB1sp_location)]
#generate sc_count and sc_meta
scRB
scRB@meta.data$celltype
scRB@assays$RNA
scRB@active.ident
scRB@meta.data$orig.ident
count.data <- GetAssayData(object = scRB[["RNA"]], slot = "counts")
count.data
sc_count <- count.data
cellType <- scRB@meta.data$celltype
cellType
cellID <- scRB@active.ident
cellID
sampleInfo <- scRB@meta.data$orig.ident
sampleInfo
meta.data <- data.frame(cellID, cellType, sampleInfo)
meta.data$cellID <- NULL
meta.data
meta.data <- tibble::rownames_to_column(meta.data, "cellID")
rownames(meta.data) <- meta.data$cellID
sc_meta <- meta.data
print(dim(sc_count))
print(dim(sc_meta))
table(sc_meta$cellType)
##get rid of Neural Cell and Others in the sc dataset to re-do the deconvonlution
scRB <- readRDS("RB_data.RDS")
scRB_updated = UpdateSeuratObject(object = scRB)
scRB = UpdateSeuratObject(object = scRB)
scRB$celltype
celltype <- c("CP","HP-CP","Other","Neural cell","CAF","Rod-like","Glial","Cone-like")
names(celltype) <- levels(scRB_updated)
Idents(scRB_updated) <- scRB_updated$celltype
scRB_updated <- RenameIdents(scRB_updated, celltype)
DimPlot(object = scRB_updated, reduction = "umap", label = T)
scRB_updated <- subset(x = scRB_updated, idents = c("Neural cell", "Other"), invert = TRUE)
scRB_updated$celltype
scRB
scRB_updated
DimPlot(object = scRB_updated, reduction = "umap", label = T)
remove(scRB)
scRB_updated@meta.data$celltype
scRB_updated@assays$RNA
scRB_updated@active.ident
scRB_updated@meta.data$orig.ident
count.data <- GetAssayData(object = scRB_updated[["RNA"]], slot = "counts")
count.data
sc_count <- count.data
cellType <- scRB_updated@meta.data$celltype
cellType
cellID <- scRB_updated@active.ident
cellID
sampleInfo <- scRB_updated@meta.data$orig.ident
sampleInfo
meta.data <- data.frame(cellID, cellType, sampleInfo)
meta.data$cellID <- NULL
meta.data
meta.data <- tibble::rownames_to_column(meta.data, "cellID")
rownames(meta.data) <- meta.data$cellID
sc_meta <- meta.data
print(dim(sc_count))
print(dim(sc_meta))
table(sc_meta$cellType)
###RB4
###create CARD objects
CARD_RB4obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = RB4sp_count,
  spatial_location = RB4sp_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5)
#Run deconvolution
CARD_RB4object = CARD_deconvolution(CARD_RB4obj)
##save results
saveRDS(CARD_RB4obj, "CARD_RB4obj_modified.RDS")
saveRDS(CARD_RB4object, "CARD_RB4object_modified.RDS")
##visualization
print(CARD_RB4object@Proportion_CARD[1:2,])
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p4 <- CARD.visualize.pie(proportion = CARD_RB4object@Proportion_CARD,spatial_location = CARD_RB4object@spatial_location, colors = colors)
print(p4)
###RB3
#create CARD objects
CARD_RB3obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = RB3sp_count,
  spatial_location = RB3sp_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 
#Run deconvolution
CARD_RB3object = CARD_deconvolution(CARD_RB3obj)
##save results
saveRDS(CARD_RB3obj, "CARD_RB3obj_modified.RDS")
saveRDS(CARD_RB3object, "CARD_RB3object_modified.RDS")
##visualization
print(CARD_RB3object@Proportion_CARD[1:2,])
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p3 <- CARD.visualize.pie(proportion = CARD_RB3object@Proportion_CARD,spatial_location = CARD_RB3object@spatial_location, colors = colors)
print(p3)
##RB2
#create CARD objects
CARD_RB2obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = RB2sp_count,
  spatial_location = RB2sp_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 
#Run deconvolution
CARD_RB2object = CARD_deconvolution(CARD_RB2obj)
##save results
saveRDS(CARD_RB2obj, "CARD_RB2obj_modified.RDS")
saveRDS(CARD_RB2object, "CARD_RB2object_modified.RDS")
##visualization
print(CARD_RB2object@Proportion_CARD[1:2,])
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p2 <- CARD.visualize.pie(proportion = CARD_RB2object@Proportion_CARD,spatial_location = CARD_RB2object@spatial_location, colors = colors)
print(p2)
##RB1
#create CARD objects
CARD_RB1obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = RB1sp_count,
  spatial_location = RB1sp_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 
#Run deconvolution
CARD_RB1object = CARD_deconvolution(CARD_RB1obj)
##save results
saveRDS(CARD_RB1obj, "CARD_RB1obj_modified.RDS")
saveRDS(CARD_RB1object, "CARD_RB1object_modified.RDS")
##visualization
print(CARD_RB1object@Proportion_CARD[1:2,])
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p1 <- CARD.visualize.pie(proportion = CARD_RB1object@Proportion_CARD,spatial_location = CARD_RB1object@spatial_location, colors = colors)
print(p1)
#combine the pie charts
grid.arrange(p1, p2, p3, p4, nrow = 2)
#generate the results
RB1_CARD <- data.frame(CARD_RB1object@Proportion_CARD)
RB2_CARD <- data.frame(CARD_RB2object@Proportion_CARD)
RB3_CARD <- data.frame(CARD_RB3object@Proportion_CARD)
RB4_CARD <- data.frame(CARD_RB4object@Proportion_CARD)
write.csv(RB1_CARD, "RB1_CARD.CSV")
write.csv(RB2_CARD, "RB2_CARD.CSV")
write.csv(RB3_CARD, "RB3_CARD.CSV")
write.csv(RB4_CARD, "RB4_CARD.CSV")
spatial_cluster<-stRB$new.cluster.ids
write.csv(spatial_cluster, "spatial_cluster.CSV")

##figure 2B-E
## select the cell type that we are interested
ct.visualize = c("CP", "MKI67+ CP", "Glial","CAF")
## visualize the spatial distribution of the cell type proportion
p3_2 <- CARD.visualize.prop(
  proportion = RB3@Proportion_CARD,        
  spatial_location = RB3@spatial_location, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 2)                                 ### number of columns in the figure panel
print(p3_2)
p4_2 <- CARD.visualize.prop(
  proportion = RB4@Proportion_CARD,        
  spatial_location = RB4@spatial_location, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 2)                                 ### number of columns in the figure panel
print(p4_2)

###RNA velocity
#prepare for scVelo
#save metadata table
# save metadata table:
SP$barcode <- colnames(SP)
SP$UMAP_1 <- SP@reductions$umap@cell.embeddings[,1]
SP$UMAP_2 <- SP@reductions$umap@cell.embeddings[,2]
write.csv(SP@meta.data, file='metadata.csv', quote=F, row.names=F)
# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(SP, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0(out_data_dir, 'counts.mtx'))
# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(SP@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)
##go to RNA_velocity_pipeline.txt for codes: figure 3B, D, F, H
##cellrank go to CellRank_pipeline.txt for codes: figure 3E,G

##neighbourhood analysis
#Define neighbor region and following analysis
#C2
SP <- SetIdent(SP, value = "seurat_clusters_harmony")
SP <- RegionNeighbours(SP, id = "2", keep.within.id = T, verbose = TRUE)
#DE analysis of the cluster border
SP <- SetIdent(SP, value = "nbs_2")
nbs_2.markers <- FindMarkers(SP, ident.1 = "2", ident.2 = "nbs_2")
nbs_2.markers$gene <- rownames(nbs_2.markers)
RB.subset <- SubsetSTData(SP, expression = nbs_2 %in% c("2", "nbs_2"))
C2_sorted.marks <- nbs_2.markers %>% arrange(-avg_log2FC) %>% top_n(n = 40, wt = abs(avg_log2FC))
#C0
FeatureOverlay(SP, features = "seurat_clusters_harmony", sampleids = c(1,2,3,4), ncols = 2, pt.size = 3)
SP <- SetIdent(SP, value = "seurat_clusters_harmony")
SP <- RegionNeighbours(SP, id = 0, keep.within.id = T, verbose = TRUE)
#DE analysis of the cluster border
SP <- SetIdent(SP, value = "nbs_0")
nbs_0.markers <- FindMarkers(SP, ident.1 = "0", ident.2 = "nbs_0")
nbs_0.markers$gene <- rownames(nbs_0.markers)
RB.subset <- SubsetSTData(SP, expression = nbs_0 %in% c("0", "nbs_0"))
C0_sorted.marks <- nbs_0.markers %>% arrange(-avg_log2FC) %>% top_n(n = 40, wt = abs(avg_log2FC))

##figure 5A-B
FeatureOverlay(SP, features = "nbs_2", ncols = 2, 
               sampleids = c(1,2,3,4), cols = c("lightgray", "red"), pt.size = 3)
FeatureOverlay(SP, features = "nbs_0", ncols = 2, 
               sampleids = c(1,2,3,4), cols = c("lightgray", "red"), pt.size = 3)

##figure 5C-D
nbs0_list
nbs_0.markers
nbs_2.markers
EnhancedVolcano(nbs_0.markers,
                lab = rownames(nbs_0.markers), 
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Significant genes C0 neighbor vs C0',
                pCutoff = 0.1,
                FCcutoff = 0.58,
                pointSize = 4,
                labSize = 7,
                colAlpha = 1,
                legendLabels = c('Not Sig','Log2FC','padj','padj&Log2FC'),
                legendPosition = "bottom",
                legendLabSize = 16,
                legendIconSize = 5.0,
                drawConnectors = T,
                widthConnectors = 0.75)
EnhancedVolcano(nbs_2.markers,
                lab = rownames(nbs_2.markers), 
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Significant genes C2 neighbor vs C2',
                pCutoff = 0.1,
                FCcutoff = 0.58,
                pointSize = 4,
                labSize = 7,
                colAlpha = 1,
                legendLabels = c('Not Sig','Log2FC','padj','padj&Log2FC'),
                legendPosition = "bottom",
                legendLabSize = 16,
                legendIconSize = 5.0,
                drawConnectors = T,
                widthConnectors = 0.75)


#DE analysis
####gene set enrichment analysis with clusterprofiler
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
C0M <- FindMarkers(SP, ident.1 = "C0")
C1M <- FindMarkers(SP, ident.1 = "C1")
C2M <- FindMarkers(SP, ident.1 = "C2")
C3M <- FindMarkers(SP, ident.1 = "C3")
C4M <- FindMarkers(SP, ident.1 = "C4")
C5M <- FindMarkers(SP, ident.1 = "C5")
C6M <- FindMarkers(SP, ident.1 = "C6")
C7M <- FindMarkers(SP, ident.1 = "C7")
C8M <- FindMarkers(SP, ident.1 = "C8")
C9M <- FindMarkers(SP, ident.1 = "C9")
#nbs_2.markers <- FindMarkers(RBanno, ident.1 = "2", ident.2 = "nbs_2")
original_nbs2_list <- nbs_2.markers$avg_log2FC
names(original_nbs2_list) <- rownames(nbs_2.markers)
nbs2_list <- na.omit(original_nbs2_list)
nbs2_list = sort(nbs2_list, decreasing = TRUE)
#nbs_0.markers <- FindMarkers(RBanno, ident.1 = "0", ident.2 = "nbs_0")
original_nbs0_list <- nbs_0.markers$avg_log2FC
names(original_nbs0_list) <- rownames(nbs_0.markers)
nbs0_list <- na.omit(original_nbs0_list)
nbs0_list = sort(nbs0_list, decreasing = TRUE)
#umap.harmony clusters
original_C0_list <- C0M$avg_log2FC
names(original_C0_list) <- rownames(C0M)
C0_list <- na.omit(original_C0_list)
C0_list = sort(C0_list, decreasing = TRUE)
C0_list
original_C1_list <- C1M$avg_log2FC
names(original_C1_list) <- rownames(C1M)
C1_list <- na.omit(original_C1_list)
C1_list = sort(C1_list, decreasing = TRUE)
C1_list
original_C2_list <- C2M$avg_log2FC
names(original_C2_list) <- rownames(C2M)
C2_list <- na.omit(original_C2_list)
C2_list = sort(C2_list, decreasing = TRUE)
C2_list
original_C3_list <- C3M$avg_log2FC
names(original_C3_list) <- rownames(C3M)
C3_list <- na.omit(original_C3_list)
C3_list = sort(C3_list, decreasing = TRUE)
C3_list
original_C4_list <- C4M$avg_log2FC
names(original_C4_list) <- rownames(C4M)
C4_list <- na.omit(original_C4_list)
C4_list = sort(C4_list, decreasing = TRUE)
C4_list
original_C5_list <- C5M$avg_log2FC
names(original_C5_list) <- rownames(C5M)
C5_list <- na.omit(original_C5_list)
C5_list = sort(C5_list, decreasing = TRUE)
C5_list
original_C6_list <- C6M$avg_log2FC
names(original_C6_list) <- rownames(C6M)
C6_list <- na.omit(original_C6_list)
C6_list = sort(C6_list, decreasing = TRUE)
C6_list
original_C7_list <- C7M$avg_log2FC
names(original_C7_list) <- rownames(C7M)
C7_list <- na.omit(original_C7_list)
C7_list = sort(C7_list, decreasing = TRUE)
C7_list
original_C8_list <- C8M$avg_log2FC
names(original_C8_list) <- rownames(C8M)
C8_list <- na.omit(original_C8_list)
C8_list = sort(C8_list, decreasing = TRUE)
C8_list
original_C9_list <- C9M$avg_log2FC
names(original_C9_list) <- rownames(C9M)
C9_list <- na.omit(original_C9_list)
C9_list = sort(C9_list, decreasing = TRUE)
C9_list
#BP for biological process, MF for molecular function, CC for cell compartment, ALL for all
keytypes(org.Hs.eg.db)
C0gse <- gseGO(geneList=C0_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C1gse <- gseGO(geneList=C1_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C2gse <- gseGO(geneList=C2_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C3gse <- gseGO(geneList=C3_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C4gse <- gseGO(geneList=C4_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C5gse <- gseGO(geneList=C5_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C6gse <- gseGO(geneList=C6_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C7gse <- gseGO(geneList=C7_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C8gse <- gseGO(geneList=C8_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C9gse <- gseGO(geneList=C9_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
nbs0gse <- gseGO(geneList=nbs0_list, 
                 ont ="BP", 
                 keyType = "SYMBOL", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "none")
nbs2gse <- gseGO(geneList=nbs2_list, 
                 ont ="BP", 
                 keyType = "SYMBOL", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "none")
require(DOSE)
pCP0D <- dotplot(C0gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP0C <- cnetplot(C0gse, categorySize="pvalue", foldChange=C0_list, showCategory = 3)

pCP1D <- dotplot(C1gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP1C <- cnetplot(C1gse, categorySize="pvalue", foldChange=C1_list, showCategory = 3)

pCP2D <- dotplot(C2gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP2C <- cnetplot(C2gse, categorySize="pvalue", foldChange=C2_list, showCategory = 3)

pCP3D <- dotplot(C3gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP3C <- cnetplot(C3gse, categorySize="pvalue", foldChange=C3_list, showCategory = 3)

pCP4D <- dotplot(C4gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP4C <- cnetplot(C4gse, categorySize="pvalue", foldChange=C4_list, showCategory = 3)

pCP5D <- dotplot(C5gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP5C <- cnetplot(C5gse, categorySize="pvalue", foldChange=C5_list, showCategory = 3)

pCP6D <- dotplot(C6gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP6C <- cnetplot(C6gse, categorySize="pvalue", foldChange=C6_list, showCategory = 3)

pCP7D <- dotplot(C7gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP7C <- cnetplot(C7gse, categorySize="pvalue", foldChange=C7_list, showCategory = 3)

pCP8D <- dotplot(C8gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP8C <- cnetplot(C8gse, categorySize="pvalue", foldChange=C8_list, showCategory = 3)

pCP9D <- dotplot(C9gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP9C <- cnetplot(C9gse, categorySize="pvalue", foldChange=C9_list, showCategory = 3)

pCPnbs0D <- dotplot(nbs0gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCPnbs0C <- cnetplot(nbs0gse, categorySize="pvalue", foldChange=nbs0_list, showCategory = 3)

pCPnbs2D <- dotplot(nbs2gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCPnbs2C <- cnetplot(nbs2gse, categorySize="pvalue", foldChange=nbs2_list, showCategory = 3)

##figure 4B-E
pCP0D
pCP0C
pCP2D
pCP2C

##figure 5E-F
pCPnbs0C
pCPnbs2C

##supplementary figure 5
pCP0D
pCP0C
pCP1D
pCP1C
pCP2D
pCP2C
pCP3D
pCP3C
pCP4D
pCP4C
pCP5D
pCP5C
pCP6D
pCP6C
pCP7D
pCP7C
pCP8D
pCP8C
pCP9D
pCP9C

###SCENIC (Single-Cell rEgulatory Network Inference and Clustering)
#input from seruat objects
#read the datasets
#read the spatial dataset
RB <- SP
harmony_cluster <- c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9")
names(harmony_cluster) <- levels(RB)
RB <- RenameIdents(RB, harmony_cluster)
DimPlot(RB, reduction = "umap.harmony", label = TRUE, pt.size = 1.8, label.size = 10)
## Get data from seurat object:
exprMat <- as.matrix(RB@assays$RNA@data)
dim(exprMat)
exprMat[1:4,1:4]
cellInfo <- RB@meta.data[,c(15,2,3)]
cellInfo
colnames(cellInfo) = c('Harmony_Cluster', 'nGene', 'uUMI')
head(cellInfo)
table(cellInfo$Harmony_Cluster)
#initialize settings
library(SCENIC)
#reference database
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather")
# mc9nr: Motif collection version 9: 24k motifs
# dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
options(timeout=1000)
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}
#input the reference model
org <- "hgnc"
dbDir <- "cisTarget_databases"
dbDir <- path.expand(dbDir)
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
dbs
lapply(dbs,function(x) file.exists(file.path(dbDir, x)))
myDatasetTitle <- "SCENIC example on Retinoblastoma"
motifAnnotations_hgnc <- motifAnnotations
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=1)
#alternatively use the new dataset
dbs[dbs == "hg19-500bp-upstream-7species.mc9nr.feather"] <- "hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather"
dbs[dbs == "hg19-tss-centered-10kb-7species.mc9nr.feather"] <- "hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather"
scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases", dbs=dbs, nCores=1)
saveRDS(scenicOptions, file = "int/scenicOptions.RDS")
### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)
### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
library(foreach)
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 50
scenicOptions@settings$defaultTsne$perpl <- 9
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
savedSelections <- shiny::runApp(aucellApp)
# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
scenicOptions@settings$devType="png"
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings
nPcs <- c(5,15,50)
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=T, cex=.5)
par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5)
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 32
scenicOptions@settings$defaultTsne$perpl <- 09
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
#save the result model
exprMat <- get_dgem(open_loom(loomPath))
dgem <- exprMat
head(colnames(dgem))
scenicOptions@fileNames$output["loomFile",] <- "output/mouseBrain_SCENIC.loom"
export2loom(scenicOptions, exprMat)
# Export:
scenicOptions@fileNames$output["loomFile",] <- "output/mouseBrain_SCENIC.loom"
export2loom(scenicOptions, exprMat)
devtools::install_github("aertslab/SCopeLoomR")
library(SCopeLoomR)
dgem <- exprMat
head(colnames(dgem))
scenicOptions@fileNames$output["loomFile",] <- "output/RB_SCENIC.loom"
export2loom(scenicOptions, exprMat, addAllTsnes = T, hierarchy=c("SCENIC", "RB"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
#use the save data to explore the results
library(Seurat) 
library(SCENIC)
library(doParallel)
scenicOptions=readRDS(file="int/scenicOptions.Rds")
### Exploring output 
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org
#result if TF enrichment
# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes") 
as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T))
#number of motifs of every gene
as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T))
#visualize the motif characterisation of specific gene
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="ETS1"]
viewMotifs(tableSubset)
#motif characterisation after adding regulons
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="ETS1" & highConfAnnot==TRUE]
#use the saved figures
library(Seurat) 
library(SCENIC)
library(AUCell)
library(doParallel)
library(SCopeLoomR)
scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicLoomPath <- getOutName(scenicOptions, "RB_SCENIC.loom")
loom <- open_loom("output/RB_SCENIC.loom")
# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
viewMotifs(tableSubset)
#density in tSNE
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Dlx5", "Sox10", "Sox9","Irf1", "Stat6")],], plots="Expression")
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)
#show multiple regulons
par(mfrow=c(1,2))
regulonNames <- c( "ETS1","SOX4")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)

regulonNames <- list(red=c("Sox10", "Sox8"),
                     green=c("Irf1"),
                     blue=c( "Tef"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")
#GRN
regulons <- loadInt(scenicOptions, "regulons")
regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
#step2 alternative
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
viewMotifs(motifEnrichment_selfMotifs_wGenes)
tableSubset1 <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="CEBPD"]
viewMotifs(tableSubset1)
tableSubset2 <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="ETS1"]
viewMotifs(tableSubset2) 
tableSubset3 <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="JUNB"]
viewMotifs(tableSubset3) 
tableSubset4 <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="JUND"]
viewMotifs(tableSubset4) 
tableSubset5 <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="SOX4"]
viewMotifs(tableSubset5) 
tableSubset6 <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="SREBF2"]
viewMotifs(tableSubset6) 
tableSubset7 <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="TCF4"]
viewMotifs(tableSubset7) 

##figure 4F
#average regulon activity in clusters
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Harmony_Cluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity", cluster_columns = F)


##Cellchat
#read the datasets
# Prepare input data for CellChat analysis
new.cluster.ids <- c("C0", "C1", "C2", "C3", "C4", 
                     "C5",
                     "C6", "C7", "C8", "C9")
names(new.cluster.ids) <- levels(SP)
SP <- RenameIdents(SP, new.cluster.ids)
DimPlot(SP, reduction = "umap.harmony", label = TRUE, pt.size = 1.3, label.size = 4.5)
data.input = GetAssayData(SP, slot = "data", assay = "SCT")
meta = data.frame(labels = Idents(SP), row.names = names(Idents(SP)))
###use the whole transcriptomic profile for cell chat
cellchat.all <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat.all <- addMeta(cellchat.all, meta = meta)
cellchat.all <- setIdent(cellchat.all, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat.all@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat.all@idents)) # number of cells in each cell group
##set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat.all@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat.all <- subsetData(cellchat.all)
future::plan("multisession", workers = 14)
#pre-processing of the expression data for cell-cell communication
cellchat.all <- identifyOverExpressedGenes(cellchat.all)
cellchat.all <- identifyOverExpressedInteractions(cellchat.all)
#The communication probability and infer cellular communication network
cellchat.all <- computeCommunProb(cellchat.all)
cellchat.all <- filterCommunication(cellchat.all, min.cells = 10)
cellchat.all <- computeCommunProbPathway(cellchat.all)
cellchat.all <- aggregateNet(cellchat.all)
saveRDS(cellchat.all, file = "cellchat.rds")

##figure 6A
groupSize <- as.numeric(table(cellchat.all@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.all@net$count, vertex.weight = rowSums(cellchat.all@net$count), 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions", edge.width.max = 14)
netVisual_circle(cellchat.all@net$weight, vertex.weight = rowSums(cellchat.all@net$weight), 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", edge.width.max = 14)

##figure 6B
#Manifold and classification learning analysis of signaling networks
#Identify signaling groups based on their functional similarity
cellchat.all <- computeNetSimilarity(cellchat.all, type = "functional")
cellchat.all <- netEmbedding(cellchat.all, slot.name = "netP", type = "functional")
#Manifold learning of the signaling networks for a single dataset
cellchat.all <- netClustering(cellchat.all, type = "functional")
#Classification learning of the signaling networks for a single dataset
#Visualization in 2D-space
pfun <- netVisual_embedding(cellchat.all, type = "functional", label.size = 3.5)
pfunzoom <- netVisual_embeddingZoomIn(cellchat.all, type = "functional", nCol = 2)
pfun
pfunzoom

##figure 6C
pathways.show <- c("NOTCH")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.all, signaling = pathways.show, layout = "chord")
pathways.show <- c("VEGF")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.all, signaling = pathways.show, layout = "chord")
pathways.show <- c("APP")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.all, signaling = pathways.show, layout = "chord")
pathways.show <- c("FN1")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.all, signaling = pathways.show, layout = "chord")
# Compute the network centrality scores
cellchat.all <- netAnalysis_computeCentrality(cellchat.all, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.all, signaling = "NOTCH", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.all, signaling = "VEGF", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.all, signaling = "APP", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.all, signaling = "FN1", width = 8, height = 2.5, font.size = 10)

#figure 6D
#Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.all, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.all, pattern = "incoming")
ht1 + ht2

##figure 6E
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat.all, sources.use = c("C0", "C1", "C3", "C4", "C5", "C6", "C7", "C8", "C9"), targets.use = "C2", remove.isolate = FALSE)
netVisual_bubble(cellchat.all, sources.use = "C0", targets.use = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"), remove.isolate = FALSE)

##supplementary figure 6
netVisual_bubble(cellchat.all, sources.use = "C2", targets.use = c("C0", "C1", "C3", "C4", "C5", "C6", "C7", "C8", "C9"), remove.isolate = FALSE)
netVisual_bubble(cellchat.all, sources.use = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"), targets.use = "C0", remove.isolate = FALSE)

######save session info
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
####end of the session####
####author: Raymond Wong, Luozixian Wang######
####05/02/2024################################

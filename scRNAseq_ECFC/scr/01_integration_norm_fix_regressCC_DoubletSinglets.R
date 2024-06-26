# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)
library(ggrepel)

# read in the data --------------------------------------------------------
# load the LUT
LUT <- read_tsv("data/LUT_sample.tsv")

# in this case wa are going to use the fix threshold filtered data
df_file <- dir("../out_large/Belfast_scRNAseq_ECFC_senescence/object/") %>%
  str_subset(pattern = c("datasc_fix_filter_norm_doublet_")) %>%
  str_subset(pattern = c("Singlet")) %>%
  data.frame(sample_file = .) %>%
  mutate(sample_id = str_remove_all(sample_file,pattern = "datasc_fix_filter_norm_doublet_Singlet_|.rds")) %>%
  left_join(.,LUT,by = "sample_id")

# set them up in a list
data.list <- lapply(df_file$sample_file, function(x){
  readRDS(paste0("../out_large/Belfast_scRNAseq_ECFC_senescence/object/",x))
}) %>%
  setNames(df_file$sample)

# here the steps recommended by the seurat SCTransform integration
# save the list of features for the integration
features <- SelectIntegrationFeatures(object.list = data.list)
combined.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)

# this command creates an 'integrated' data assay
data.combined <- IntegrateData(anchorset = combined.anchors)

# fix the meta
# data.combined$treat <- factor(data.combined$treat,levels = c("CTRL","CSF"))
# data.combined$infection <- factor(data.combined$infection,levels = c("WT","SOX10"))

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
# # Run the standard workflow for visualization and clustering
# data.combined <- ScaleData(data.combined, verbose = FALSE)
# add the cell cycle analysis
DefaultAssay(data.combined) <- "RNA"

# convert them to mouse gene
# library(homologene)
# s.genes <- homologene(cc.genes$s.genes, inTax = 9606, outTax = 10090) %>% pull("10090")
# g2m.genes <- homologene(cc.genes$g2m.genes, inTax = 9606, outTax = 10090) %>% pull("10090")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m

data.combined <- CellCycleScoring(data.combined, s.features = s.genes, g2m.features = g2m.genes)
# head(data.combined@meta.data)
DefaultAssay(data.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
# data.combined <- ScaleData(data.combined,vars.to.regress = c("percent.mt","nCount_RNA"))
data.combined <- ScaleData(data.combined,vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"))
# data.combined <- ScaleData(data.combined,vars.to.regress = c("percent.mt","nCount_RNA"))
# data.combined <- ScaleData(data.combined)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 0.3)

# show all the varibale available in the metadata
head(data.combined@meta.data)

# show the data after integration
DimPlot(data.combined, reduction = "umap", group.by = "DF")
ggsave("out/image/UMAP_integration_cluster_norm_fix_CellCycle_label_regressCC_DoubletSinglet.pdf",width = 5,height = 4)

df_meta <- data.combined@meta.data %>%
  data.frame() %>%
  rownames_to_column("barcode")
df_umap <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")

data2 <- left_join(df_umap,df_meta,"barcode")
data2_avg <- data2 %>% group_by(seurat_clusters) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=seurat_clusters),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data2_avg,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~orig.ident,nrow = 1)
ggsave("out/image/UMAP_integration_cluster_norm_fix_regressCC_DoubletSinglet_label_split.pdf",width = 12,height = 3)

ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = data2 %>%
               arrange(desc(DF)),aes(x = UMAP_1,y = UMAP_2,col=DF),size=0.3,alpha=0.5) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = data2_avg,aes(x = UMAP_1,y = UMAP_2,label = DF),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~orig.ident,nrow = 1)
ggsave("out/image/UMAP_integration_cluster_norm_fix_regressCC_DoubletSinglet_split.pdf",width = 12,height = 3)

# DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters",label = T)

# DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters",split.by = "seurat_clusters")

# save the table of normalized expression for the integrated data
write_tsv(data.combined@meta.data %>%
            rownames_to_column("cell_id"),"out/table/data.combined_fix_normalized_meta_regressCC_DoubletSinglet.tsv")
#
write_tsv(data.combined@assays$RNA@data %>%
            as.matrix() %>%
            data.frame() %>%
            rownames_to_column(var = "gene"),"../out_large/Belfast_scRNAseq_ECFC_senescence/table/data.combined_fix_normalized_data_regressCC_DoubletSinglet.tsv")

# at this point is possible to move to SCINA for the classification of the cell clusters
# data.combined@assays$RNA@data

# Identify conserved cell type markers ------------------------------------
# For performing differential expression after integration, we switch back to the original
# data
# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(data.combined) <- "RNA"
# notice that in this case the data object of the RNA slot is already filled with the normalzied data, therefore in this case (following the Normalize workfloe for the integration) there is no need to run the NormalizeData on the RNA slot of the integrated object
data.combined@assays$RNA@data[20:50,1:10]

# scale the data see the note on evernote on why this can be done also at this point. this is needed because the scale.data is empty
data.combined@assays$RNA@scale.data
all.genes <- rownames(data.combined)
data.combined <- ScaleData(data.combined,vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"))
# confirm now the slot is loaded
data.combined@assays$RNA@scale.data[1:10,1:10]

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
# data.combined.markers <- FindAllMarkers(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# # save the table of all markers
# data.combined.markers %>%
#   write_tsv("out/fixed_threshold_doublet//FindAllMarkers_data.combined_fix_regressCC_DoubletSinglet.tsv")

# # try the other approach using the
# # the issue is that some of the clusters have very few cells in the old condition that is throwing an error
# 
# # run the identification of the markers across all the clusters
# cluster_id <- sort(unique(data.combined@meta.data$seurat_clusters))
# 
# data.combined@meta.data %>%
#   group_by(seurat_clusters,orig.ident) %>%
#   summarise(n = n())
# 
# # in this case there is one sample per condition, as the treat and infection is crossed. the best grouping varaible in this case would be the orig.ident
# list_conserved_markers <- lapply(cluster_id, function(x){
#   FindConservedMarkers(data.combined, ident.1 = x, grouping.var = "orig.ident", verbose = T,min.cells.group = 0)
# }) %>%
#   setNames(cluster_id) %>%
#   bind_rows(.id = "cluster_number")
# 
# # save the table of all markers
# list_conserved_markers %>%
#   rownames_to_column(var = "x") %>%
#   separate(x,into = c("gene","id"),sep = "\\.") %>%
#   write_tsv("out/FindConservedMarkers_data.combined_fix_regress.tsv")
# both table can be used for the annotation of the clusters
# save the partial combined object before finalizing the final annotation
saveRDS(object = data.combined,file = "../out_large/Belfast_scRNAseq_ECFC_senescence/object/data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.rds")
data.combined <- readRDS(file = "../out_large/Belfast_scRNAseq_ECFC_senescence/object/data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.rds")

# -------------------------------------------------------------------------
# main umap
DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters",label = T)
ggsave("out/image/UMAP_integration_cluster_norm_fix_subcluster_label_regressCC_DoubletSinglet.pdf",width = 5,height = 4)

# main umap cell cycle
DimPlot(data.combined, reduction = "umap", group.by = "Phase")
ggsave("out/image/UMAP_integration_cluster_norm_fix_subcluster_phase_regressCC_DoubletSinglet.pdf",width = 5,height = 4)

# split by sample
df_meta <- data.combined@meta.data %>%
  rownames_to_column("barcode")
df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")

data2 <- left_join(df_UMAP,df_meta,"barcode")
data2_avg <- data2 %>% group_by(seurat_clusters) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# single plot
data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("out/image/UMAP_integration_cluster_norm_fix_subcluster_label_regressCC_DoubletSinglet_ggplot.pdf",width = 5,height = 3)

# split by sample
data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  facet_grid(~orig.ident)+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("out/image/UMAP_integration_cluster_norm_fix_subcluster_label_split_regressCC_DoubletSinglet.pdf",width = 12,height = 3)

# proportion of cell per cluster
df_summary <- df_meta %>%
  group_by(orig.ident,seurat_clusters) %>%
  summarise(n = n()) %>%
  group_by(orig.ident) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)

write_tsv(df_summary,"out/table/summary_data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.tsv")

# df_summary %>%
#   mutate(orig.ident = factor(orig.ident)) %>%
#   # mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
#   ggplot(aes(x=seurat_clusters,y=prop,fill=orig.ident))+geom_col()+theme_bw()
# ggsave("out/image/summary_stacked_data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.pdf",width = 5,height = 4)

df_summary %>%
  # mutate(orig.ident = factor(orig.ident,levels = c("Young_1","Adult_1","Old_1"))) %>%
  # mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  ggplot(aes(x=seurat_clusters,y=prop,fill=orig.ident))+geom_col()+theme_bw()+
  scale_fill_manual(values = c("#36AAF9FF","#4662D7FF","#FABA39FF","#F66B19FF"))
ggsave("out/image/summary_stacked_data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.pdf",width = 5,height = 4)

# df_summary %>%
#   mutate(orig.ident = factor(orig.ident)) %>%
#   ggplot(aes(x=seurat_clusters,y=prop,fill=orig.ident))+geom_col(position = "dodge")+theme_bw()
# ggsave("out/image/summary_dodge_data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.pdf",width = 8,height = 4)

df_summary %>%
  # mutate(orig.ident = factor(orig.ident,levels = c("Young_1","Adult_1","Old_1"))) %>%
  ggplot(aes(x=seurat_clusters,y=prop,fill=orig.ident))+geom_col(position = "dodge")+theme_bw()+
  scale_fill_manual(values = c("#36AAF9FF","#4662D7FF","#FABA39FF","#F66B19FF"))
ggsave("out/image/summary_dodge_data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.pdf",width = 8,height = 4)

shortlist_features_list_short <- list(
  IMMUNE = c("IGHG1","CD38","CD8A","CD2","SKAP1","LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1"),
  ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
  PANEL_1 = c("CD34","PROCR","BST1","NRP1")
)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_short <- DotPlot(data.combined, features = shortlist_features_list_short, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
test_short
ggsave("out/image/Dotplot_clusters_NOT_annotated_short_panle_integration_norm_fix_long_regressCC_DoubletSinglet.pdf",width = 12,height = 4)

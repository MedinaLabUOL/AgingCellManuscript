# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)
library(ggrepel)
library(cowplot)

# read in the data --------------------------------------------------------
# load the LUT
LUT <- read_csv("data/LUT_sample.csv")

# in this case wa are going to use the fix threshold filtered data
df_file <- dir("../out_large/Belfast_scRNAseq_ECRetinaMouse/object/") %>%
  str_subset(pattern = c("datasc_fix_filter_norm_doublet_")) %>%
  str_subset(pattern = c("Singlet")) %>%
  data.frame(sample_file = .) %>%
  mutate(sample_id = str_remove_all(sample_file,pattern = "datasc_fix_filter_norm_doublet_Singlet_|.rds")) %>%
  left_join(.,LUT,by = "sample_id")

# set them up in a list
data.list <- lapply(df_file$sample_file, function(x){
  readRDS(paste0("../out_large/Belfast_scRNAseq_ECRetinaMouse/object/",x))
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
library(homologene)
s.genes <- homologene(cc.genes$s.genes, inTax = 9606, outTax = 10090) %>%
  pull("10090")
g2m.genes <- homologene(cc.genes$g2m.genes, inTax = 9606, outTax = 10090) %>%
  pull("10090")

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
data.combined <- FindClusters(data.combined, resolution = 0.2)

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
  ggrepel::geom_text_repel(data = data2_avg,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+
  theme_bw()+
  facet_wrap(~orig.ident,nrow = 1)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2))
ggsave("out/image/UMAP_integration_cluster_norm_fix_regressCC_DoubletSinglet_label_split.pdf",width = 10,height = 3)

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
  facet_wrap(~orig.ident,nrow = 1)+
theme(strip.background = element_blank(), 
      panel.border = element_rect(colour = "black", fill = NA))
ggsave("out/image/UMAP_integration_cluster_norm_fix_regressCC_DoubletSinglet_split.pdf",width = 10,height = 3)

# DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters",label = T)

# DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters",split.by = "seurat_clusters")

# save the table of normalized expression for the integrated data
write_tsv(data.combined@meta.data %>%
            rownames_to_column("cell_id"),"out/table/data.combined_fix_normalized_meta_regressCC_DoubletSinglet.tsv")
#
write_tsv(data.combined@assays$RNA@data %>%
            as.matrix() %>%
            data.frame() %>%
            rownames_to_column(var = "gene"),"../out_large/Belfast_scRNAseq_ECRetinaMouse/table/data.combined_fix_normalized_data_regressCC_DoubletSinglet.tsv")

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
saveRDS(object = data.combined,file = "../out_large/Belfast_scRNAseq_ECRetinaMouse/object/data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.rds")
data.combined <- readRDS(file = "../out_large/Belfast_scRNAseq_ECRetinaMouse/object/data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.rds")

# plots for elisa ---------------------------------------------------------
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

data2 %>%
  mutate(orig.ident = factor(orig.ident,levels = c("Young_1","Adult_1","Old_1"))) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  facet_grid(~orig.ident)+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2))

ggsave("out/image/UMAP_integration_cluster_norm_fix_subcluster_label_split_regressCC_DoubletSinglet.pdf",width = 11,height = 3)

data2 %>%
  mutate(orig.ident = factor(orig.ident,levels = c("Young_1","Adult_1","Old_1"))) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2))
ggsave("out/image/UMAP_integration_cluster_norm_fix_subcluster_label_regressCC_DoubletSinglet2.pdf",width = 5,height = 3)

data2 %>%
  mutate(orig.ident = factor(orig.ident,levels = c("Young_1","Adult_1","Old_1"))) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=Phase),alpha=0.5,size=0.1)+
  # geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1))
ggsave("out/image/UMAP_integration_cluster_norm_fix_subcluster_phase_regressCC_DoubletSinglet2.pdf",width = 4,height = 3)

# proportion of cell per cluster
df_summary <- df_meta %>%
  group_by(orig.ident,seurat_clusters) %>%
  summarise(n = n()) %>%
  group_by(orig.ident) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)

write_tsv(df_summary,"out/table/summary_data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.tsv")

df_summary %>%
  mutate(orig.ident = factor(orig.ident,levels = c("Young_1","Adult_1","Old_1"))) %>%
  # mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  ggplot(aes(x=seurat_clusters,y=prop,fill=orig.ident))+geom_col()+theme_bw()
ggsave("out/image/summary_stacked_data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.pdf",width = 5,height = 4)

df_summary %>%
  mutate(orig.ident = factor(orig.ident,levels = c("Young_1","Adult_1","Old_1"))) %>%
  ggplot(aes(x=seurat_clusters,y=prop,fill=orig.ident))+geom_col(position = "dodge")+theme_bw()
ggsave("out/image/summary_dodge_data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.pdf",width = 8,height = 4)

shortlist_features_list_long <- list("ROD" = c("Rho","Nrl","Nr2e3","Cngb1","Cplx4","Cnga1","Nxnl1"), 
                                     "CONE" = c("Opn1sw","Arr3","Gnat2","Thrb"),
                                     "BP" = c("Cabp5","Grm6","Isl1","Vsx2","Lhx3"),
                                     "AMA" = c("Pax6","Gad2","Nrxn2","Lrrn3"),
                                     "RGC" = c("Nrn1"),
                                     "MULLER" = c("S100b","Gpr37","Dkk3"),
                                     "MG" = c("Mertk","Siglech","P2ry13","Cx3cr1","Tyrobp","Trem2","Gpr34"),
                                     "MAC_inf" = c("H2-Aa","Cd74","H2-Ab1","H2-Eb1","Il1b"),
                                     "MONO" = c("Lyz2","Lgals3"),
                                     "LYMPHO" = c("Cd3g","Cd3d","Cd3e","Nkg7","Klrd1","Klrc2"),
                                     "ENDO" = c("Vwf","Cdh5","Tek","Pecam1","Flt1","Kdr","Nos3","Mcam","Mmrn1","Cldn5","Bmx","Angpt2","Gja4","Tie1","Robo4","Ecscr"),
                                     "PERiCYTE" = c("Pdgfrb","Des","Cspg4","Acta2","Rgs5","Abcc9","Kcnj8"),
                                     "FIB" = c("Col14a1","Col1a1","Itga8","Nr2f2","Cygb","Dcn","Prrx1","Aebp1","Meg3"),
                                     "RPE" = c("Ttr","Ptgds","Rlbp1","Lrat","Rpe65","Rdh5","Krt18","Rdh10","Bmp4","Rbp1"))

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_long <- DotPlot(data.combined, features = shortlist_features_list_long, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
test_long
ggsave("out/image/Dotplot_clusters_NOT_annotated_big_panle_integration_norm_fix_long_regressCC_DoubletSinglet.pdf",width = 25,height = 4)

# attempt annotation ------------------------------------------------------
# add new annotation based on the dotplot
id_cluster <- data.combined@meta.data$seurat_clusters
data.combined$cell_id <- case_when(id_cluster == 4 ~ "Peri",
                                   # id_cluster == 5 ~ "Fib",
                                   id_cluster == 11 ~ "RPE",
                                   id_cluster %in% c(0,1,12) ~ "Endo",
                                   id_cluster == 3 ~ "Muller",
                                   id_cluster == 8 ~ "MAC_inf",
                                   id_cluster == 13 ~ "IMM_13",
                                   id_cluster == 14 ~ "IMM_14",
                                   id_cluster %in% c(6) ~ "BP_6",
                                   id_cluster %in% c(7) ~ "BP_7",
                                   id_cluster %in% c(9) ~ "BP_9",
                                   id_cluster %in% c(10) ~ "BP_10",
                                   id_cluster == 2 ~ "Rod",
                                   id_cluster == 5 ~ "Cone")

saveRDS(object = data.combined,file = "../out_large/Belfast_scRNAseq_ECRetinaMouse/object/data.combined_annotated_norm_fix_regressCC_DoubletSinglet.rds")
data.combined <- readRDS(file = "../out_large/Belfast_scRNAseq_ECRetinaMouse/object/data.combined_annotated_norm_fix_regressCC_DoubletSinglet.rds")

# plot after annotation ---------------------------------------------------
# main umap
DimPlot(data.combined, reduction = "umap", group.by = "cell_id",label = T)
ggsave("out/image/UMAP_integration_cellID_norm_fix_subcluster_label_regressCC_DoubletSinglet.pdf",width = 5,height = 4)

# split by sample
df_meta <- data.combined@meta.data %>%
  rownames_to_column("barcode")
df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")

data2 <- left_join(df_UMAP,df_meta,"barcode")
data2_avg <- data2 %>% group_by(cell_id) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

data2 %>%
  mutate(orig.ident = factor(orig.ident,levels = c("Young_1","Adult_1","Old_1"))) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=cell_id),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = cell_id)) +
  facet_grid(~orig.ident)+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2))

ggsave("out/image/UMAP_integration_cellID_norm_fix_subcluster_label_split_regressCC_DoubletSinglet.pdf",width = 11,height = 3)

data2 %>%
  mutate(orig.ident = factor(orig.ident,levels = c("Young_1","Adult_1","Old_1"))) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=cell_id),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = cell_id)) +
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2))
ggsave("out/image/UMAP_integration_cellID_norm_fix_subcluster_label_regressCC_DoubletSinglet2.pdf",width = 5,height = 3)

# proportion of cell per cluster
df_summary <- df_meta %>%
  group_by(orig.ident,cell_id) %>%
  summarise(n = n()) %>%
  group_by(orig.ident) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)

write_tsv(df_summary,"out/table/summary_data.combined_annotated_norm_fix_regressCC_DoubletSinglet.tsv")

# library(scales) 
# show_col(viridis::turbo(10))
# viridis::turbo(10)

df_summary %>%
  mutate(orig.ident = factor(orig.ident,levels = c("Young_1","Adult_1","Old_1"))) %>%
  group_by(cell_id) %>%
  mutate(sum_prop = sum(prop)) %>%
  ungroup()%>%
  mutate(cell_id = fct_reorder(cell_id,sum_prop,.desc = T)) %>%
  ggplot(aes(x=cell_id,y=prop,fill=orig.ident))+geom_col()+theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  scale_fill_manual(values = c("#4662D7FF","#FABA39FF","#7A0403FF"))
ggsave("out/image/summary_stacked_data.combined_annotated_norm_fix_regressCC_DoubletSinglet.pdf",width = 5,height = 4)

df_summary %>%
  mutate(orig.ident = factor(orig.ident,levels = c("Young_1","Adult_1","Old_1"))) %>%
  group_by(cell_id) %>%
  mutate(sum_prop = sum(prop)) %>%
  ungroup()%>%
  mutate(cell_id = fct_reorder(cell_id,sum_prop,.desc = T)) %>%
  ggplot(aes(x=cell_id,y=prop,fill=orig.ident))+geom_col(position = "dodge")+theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  scale_fill_manual(values = c("#4662D7FF","#FABA39FF","#7A0403FF"))
ggsave("out/image/summary_dodge_data.combined_annotated_norm_fix_regressCC_DoubletSinglet.pdf",width = 8,height = 4)

cell_id_order <- df_summary %>%
  mutate(orig.ident = factor(orig.ident,levels = c("Young_1","Adult_1","Old_1"))) %>%
  group_by(cell_id) %>%
  summarise(sum_prop = sum(prop)) %>%
  arrange(desc(sum_prop)) %>%
  pull(cell_id)

shortlist_features_list_long <- list("ROD" = c("Rho","Nrl","Nr2e3","Cngb1","Cplx4","Cnga1","Nxnl1"), 
                                     "CONE" = c("Opn1sw","Arr3","Gnat2","Thrb"),
                                     "BP" = c("Cabp5","Grm6","Isl1","Vsx2","Lhx3"),
                                     "AMA" = c("Pax6","Gad2","Nrxn2","Lrrn3"),
                                     "RGC" = c("Nrn1"),
                                     "MULLER" = c("S100b","Gpr37","Dkk3"),
                                     "MG" = c("Mertk","Siglech","P2ry13","Cx3cr1","Tyrobp","Trem2","Gpr34"),
                                     "MAC_inf" = c("H2-Aa","Cd74","H2-Ab1","H2-Eb1","Il1b"),
                                     "MONO" = c("Lyz2","Lgals3"),
                                     "LYMPHO" = c("Cd3g","Cd3d","Cd3e","Nkg7","Klrd1","Klrc2"),
                                     "ENDO" = c("Vwf","Cdh5","Tek","Pecam1","Flt1","Kdr","Nos3","Mcam","Mmrn1","Cldn5","Bmx","Angpt2","Gja4","Tie1","Robo4","Ecscr"),
                                     "PERiCYTE" = c("Pdgfrb","Des","Cspg4","Acta2","Rgs5","Abcc9","Kcnj8"),
                                     "FIB" = c("Col14a1","Col1a1","Itga8","Nr2f2","Cygb","Dcn","Prrx1","Aebp1","Meg3"),
                                     "RPE" = c("Ttr","Ptgds","Rlbp1","Lrat","Rpe65","Rdh5","Krt18","Rdh10","Bmp4","Rbp1"))

# notice that this is done only on the subset of the young (control) cells
Idents(data.combined) <- "cell_id"
test_long <- DotPlot(data.combined, features = shortlist_features_list_long, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
test_long
ggsave("out/image/Dotplot_clusters_annotated_big_panle_integration_norm_fix_long_regressCC_DoubletSinglet.pdf",width = 25,height = 4)

# notice that this is done only on the subset of the young (control) cells
test <- data.combined$cell_id
test_update <- case_when(str_detect(test,pattern = "^BP")~"BP",
          T~test)
data.combined$cell_id2 <- unname(test_update)

Idents(data.combined) <- "cell_id2"
test_long2 <- DotPlot(data.combined, features = shortlist_features_list_long[1:7], dot.scale = 8,cluster.idents = T) +
  RotatedAxis()

test_long3 <- DotPlot(data.combined, features = shortlist_features_list_long[8:14], dot.scale = 8,cluster.idents = T) +
  RotatedAxis()

test_long2/test_long3

ggsave("out/image/Dotplot_clusters_annotated_big_panle_integration_norm_fix_long_regressCC_DoubletSinglet2.pdf",width = 16,height = 8)

# single plot
test <- DotPlot(data.combined, features = shortlist_features_list_long,cluster.idents = T)

df_test <- lapply(shortlist_features_list_long,function(x){
  test$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

df_test %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_size(range = c(0, 6)) +
  # facet_grid(~cell_type,scales = "free",space="free")+
  facet_wrap(~cell_type,scales = "free_x",nrow=2)+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")
ggsave("out/image/Dotplot_clusters_annotated_big_panle_integration_norm_fix_long_regressCC_DoubletSinglet3.pdf",width = 16,height = 8)


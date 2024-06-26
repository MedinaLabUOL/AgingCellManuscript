# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(circlize)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS(file = "../out_large/Belfast_scRNAseq_ECRetinaMouse/object/data.combined_annotated_norm_fix_regressCC_DoubletSinglet.rds")
# load the LUT
LUT <- read_csv("data/LUT_sample.csv")

# wrangling ---------------------------------------------------------------
# add the treatmend varibale to the metadata
meta <- data.combined@meta.data %>%
  rownames_to_column("barcodes")

# update differential expression before annotation ------------------------
# perform the DE using the cluster identity
head(data.combined)

# confirm the correct ident
Idents(data.combined) <- "seurat_clusters"

# make sure the dafault assay is RNA and make sure the data have been scaled also on the RNA slot
DefaultAssay(data.combined)
data.combined@assays$RNA@scale.data[1:10,1:10]

#
data.combined.markers <- FindAllMarkers(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
data.combined.markers %>%
  write_tsv("out/table/FindAllMarkers_data.combined_fix_regressCC_DoubletSinglet_seuratClusters.tsv")

# confirm the correct ident
Idents(data.combined) <- "cell_id"
#
data.combined.markers_cellID <- FindAllMarkers(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# save the table of all markers
data.combined.markers_cellID %>%
  write_tsv("out/table/FindAllMarkers_data.combined_fix_regressCC_DoubletSinglet_cellID.tsv")

# find conserved markers --------------------------------------------------
# try the other approach using the
# the issue is that some of the clusters have very few cells in the old condition that is throwing an error

# run the identification of the markers across all the clusters
cluster_id <- sort(unique(data.combined@meta.data$seurat_clusters))

data.combined@meta.data %>%
  group_by(seurat_clusters,orig.ident) %>%
  summarise(n = n()) %>%
  print(n=50)

Idents(data.combined) <- "seurat_clusters"
# in this case there are three conditions
# remove cluster 11 because there too few cells
cluster_id2 <- cluster_id %>%
  str_subset(pattern = "11",negate = T)

list_conserved_markers_treat <- lapply(cluster_id2, function(x){
  FindConservedMarkers(data.combined, ident.1 = x, grouping.var = "treat", verbose = T,min.cells.group = 0)
}) %>%
  setNames(cluster_id2) %>%
  bind_rows(.id = "cluster_number")

list_conserved_markers_treat %>%
  rownames_to_column("gene") %>%
  separate(gene,into = c("gene","id"),sep = "\\.\\.") %>%
  write_tsv("out/table/FindConservedMarkers_treat_data.combined_fix_regressCC_DoubletSinglet_seuratClusters.tsv")

list_conserved_markers_treat %>%
  rownames_to_column("gene") %>%
  separate(gene,into = c("gene","id"),sep = "\\.\\.") %>%
  group_by(cluster_number,gene) %>%
  mutate(avg_log2FC = mean(c(Old_avg_log2FC,Adult_avg_log2FC,Young_avg_log2FC)),
         avg_dPct = mean(c(Old_pct.1 - Old_pct.2 ,Adult_pct.1-Adult_pct.2, Young_pct.1-Young_pct.2))) %>%
  filter(avg_log2FC>0) %>%
  arrange(cluster_number,max_pval,desc(avg_log2FC)) %>%
  write_tsv("out/table/FindConservedMarkers_treat_data.combined_fix_regressCC_DoubletSinglet_seuratClusters_reordered.tsv")

# alos use the annotation from the cell_id
# run the identification of the markers across all the cell_id
cell_id <- sort(unique(data.combined@meta.data$cell_id))

data.combined@meta.data %>%
  group_by(cell_id,orig.ident) %>%
  summarise(n = n()) %>%
  print(n=50)

Idents(data.combined) <- "cell_id"
# remove cluster 11 because there too few cells
cell_id2 <- cell_id %>%
  str_subset(pattern = "RPE",negate = T)

list_conserved_markers_treat_cellID <- lapply(cell_id2, function(x){
  FindConservedMarkers(data.combined, ident.1 = x, grouping.var = "treat", verbose = T,min.cells.group = 0)
}) %>%
  setNames(cell_id2) %>%
  bind_rows(.id = "cluster_number")

list_conserved_markers_treat_cellID %>%
  rownames_to_column("gene") %>%
  separate(gene,into = c("gene","id"),sep = "\\.\\.") %>%
  write_tsv("out/table/FindConservedMarkers_treat_data.combined_fix_regressCC_DoubletSinglet_cellID.tsv")

list_conserved_markers_treat_cellID %>%
  rownames_to_column("gene") %>%
  separate(gene,into = c("gene","id"),sep = "\\.\\.") %>%
  group_by(cluster_number,gene) %>%
  mutate(avg_log2FC = mean(c(Old_avg_log2FC,Adult_avg_log2FC,Young_avg_log2FC)),
         avg_dPct = mean(c(Old_pct.1 - Old_pct.2 ,Adult_pct.1-Adult_pct.2, Young_pct.1-Young_pct.2))) %>%
  filter(avg_log2FC>0) %>%
  arrange(cluster_number,max_pval,desc(avg_log2FC)) %>%
  write_tsv("out/table/FindConservedMarkers_treat_data.combined_fix_regressCC_DoubletSinglet_cellID_reordered.tsv")

# plot the markers as an heatmap ------------------------------------------
# slect the gene to plot focus on the seuratclusters first
GOI <- read_tsv("out/table/FindConservedMarkers_treat_data.combined_fix_regressCC_DoubletSinglet_seuratClusters_reordered.tsv") %>%
  filter(str_detect(gene,pattern = "MT-|-AS1|\\.",negate = T)) %>%
  select(gene,cluster_number,max_pval,avg_log2FC) %>%
  group_by(cluster_number) %>%
  arrange(cluster_number,max_pval,avg_log2FC) %>%
  mutate(rank = rank(max_pval)) %>%
  top_n(30,wt = -rank) %>%
  pull(gene) %>%
  unique()

#produce the heatmap
test1 <- DoHeatmap(data.combined,group.by = "seurat_clusters", features = GOI) + NoLegend()
# # beware! the Identity column is not correct
data_heatmap_long <- test1$data %>%
  drop_na() %>%
  select(-Identity)

# add the original meta
data_heatmap_long_fix <- left_join(data_heatmap_long,meta,by = c("Cell"="barcodes"))

mat <- data_heatmap_long_fix %>%
  arrange(seurat_clusters,treat) %>%
  dplyr::select(Feature,Cell,Expression) %>%
  pivot_wider(names_from = Cell,values_from = Expression) %>%
  column_to_rownames("Feature")

ref_meta <- meta %>%
  slice(match(colnames(mat),meta$barcodes))

id_clusters <- names(table(meta$seurat_clusters))
id_colors <- hue_pal()(length(id_clusters))
#
names(id_colors) <- id_clusters 
column_ha <- HeatmapAnnotation(cell = ref_meta$seurat_clusters,
                               treat = ref_meta$treat,
                               col = list(treat = c("Old"="red","Adult"="yellow","Young"="Blue"),
                                          cell = id_colors))

ht <- Heatmap(mat,
              name = "test",
              use_raster = F,
              show_row_names = F,
              show_column_names = F,
              cluster_columns = F,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
              top_annotation = column_ha)

pdf("out/image/heatmap_panel_expression_FindConservedMarkers_seurat_clusters.pdf",width = 8,height = 6)
ht
dev.off()

# do the same for the cellID
# slect the gene to plot focus on the seuratclusters first
GOI_cellID <- read_tsv("out/table/FindConservedMarkers_treat_data.combined_fix_regressCC_DoubletSinglet_cellID_reordered.tsv") %>%
  filter(str_detect(gene,pattern = "MT-|-AS1|\\.",negate = T)) %>%
  select(gene,cluster_number,max_pval,avg_log2FC) %>%
  group_by(cluster_number) %>%
  arrange(cluster_number,max_pval,avg_log2FC) %>%
  mutate(rank = rank(max_pval)) %>%
  top_n(30,wt = -rank) %>%
  pull(gene) %>%
  unique()

#produce the heatmap
test1_cellID <- DoHeatmap(data.combined,group.by = "cell_id", features = GOI) + NoLegend()
# # beware! the Identity column is not correct
data_heatmap_long_cellID <- test1_cellID$data %>%
  drop_na() %>%
  select(-Identity)

# add the original meta
data_heatmap_long_fix_cellID<- left_join(data_heatmap_long_cellID,meta,by = c("Cell"="barcodes"))

mat_cellID <- data_heatmap_long_fix_cellID %>%
  arrange(cell_id,treat) %>%
  dplyr::select(Feature,Cell,Expression) %>%
  pivot_wider(names_from = Cell,values_from = Expression) %>%
  column_to_rownames("Feature")

ref_meta_cellID <- meta %>%
  slice(match(colnames(mat_cellID),meta$barcodes))

id_clusters_cellID <- names(table(meta$cell_id))
id_colors_cellID <- hue_pal()(length(id_clusters_cellID))
#
names(id_colors_cellID) <- id_clusters_cellID
column_ha_cellID <- HeatmapAnnotation(cell = ref_meta_cellID$cell_id,
                                      treat = ref_meta_cellID$treat,
                                      col = list(treat = c("Old"="red","Adult"="yellow","Young"="Blue"),
                                                 cell = id_colors_cellID))

ht_cellID <- Heatmap(mat_cellID,
                     name = "test",
                     use_raster = F,
                     show_row_names = F,
                     show_column_names = F,
                     cluster_columns = F,
                     col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                     top_annotation = column_ha_cellID)

pdf("out/image/heatmap_panel_expression_FindConservedMarkers_seurat_cellID.pdf",width = 8,height = 6)
ht_cellID
dev.off()

# run the comparison of the OLD vs PUP and ADULT vs PUP -------------------
# make suer the correct defalult dataset il loaded should be RNA
DefaultAssay(data.combined)

# check that in the RNA slot the data object is indeed loaded with normalized values
data.combined@assays$RNA@data[1:10,1:10]

# define the grouping variables for the comparison of CSF vs CTRL
head(data.combined@meta.data)

# define the new grouping
data.combined$treat.cellID <- paste(data.combined$treat,data.combined$cell_id, sep = "_")
head(data.combined@meta.data)

# update the idents of the object
Idents(data.combined) <- "treat.cellID"

# avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
# run the differential expression over the WT
clusters_id <- as.character(sort(unique(data.combined$cell_id)))

# the option below is neede to avoid errors for max size
options(future.globals.maxSize = 8000 * 1024^2)

table(data.combined@meta.data$cell_id,data.combined@meta.data$treat)

# comparison 1
# skip RPE and IMM_14 form the compariosn as there are too few cells
clusters_id2 <- clusters_id %>%
  str_subset(pattern = "IMM_14|RPE",negate = T)

list_Old_vs_Young <- lapply(clusters_id2,function(x){
  id_1 <- paste0("Old_",x)
  id_2 <- paste0("Young_",x)
  response <- FindMarkers(data.combined,
                          ident.1 = id_1,
                          ident.2 = id_2,
                          verbose = T,
                          logfc.threshold = 0,
                          min.pct = 0.1)
  response %>%
    rownames_to_column("gene") %>%
    mutate(id_1 = id_1,
           id_2 = id_2) %>%
    mutate(cluster = x)
})

list_Old_vs_Young %>%
  setNames(clusters_id2) %>%
  bind_rows(.id = "cell_id") %>%
  write_tsv("out/table/response_Old_vs_Young_data.combined_fix_regressCC_DoubletSinglet.tsv")

# comparison 2
list_Adult_vs_Young <- lapply(clusters_id2,function(x){
  id_1 <- paste0("Adult_",x)
  id_2 <- paste0("Young_",x)
  response <- FindMarkers(data.combined,
                          ident.1 = id_1,
                          ident.2 = id_2,
                          verbose = T,
                          logfc.threshold = 0,
                          min.pct = 0.1)
  response %>%
    rownames_to_column("gene") %>%
    mutate(id_1 = id_1,
           id_2 = id_2) %>%
    mutate(cluster = x)
})

list_Adult_vs_Young %>%
  setNames(clusters_id2) %>%
  bind_rows(.id = "cell_id") %>%
  write_tsv("out/table/response_Adult_vs_Young_data.combined_fix_regressCC_DoubletSinglet.tsv")

# check the position of genes francesca submitted
# GOI_df <- read_csv("data/genes_endothelial.csv")
read_tsv("out/table/response_Adult_vs_Young_data.combined_fix_regressCC_DoubletSinglet.tsv") %>%
  split(f = .$cluster) %>%
  lapply(.,function(x){
    x %>%
      mutate(dPct = pct.1 - pct.2) %>%
      filter(avg_log2FC > 0.5,
             p_val_adj < 0.01)
  })

read_tsv("out/table/response_Old_vs_Young_data.combined_fix_regressCC_DoubletSinglet.tsv") %>%
  split(f = .$cluster) %>%
  lapply(.,function(x){
    x %>%
      mutate(dPct = pct.1 - pct.2) %>%
      filter(avg_log2FC > 0.5,
             p_val_adj < 0.01)
  })

read_tsv("out/table/response_Old_vs_Young_data.combined_fix_regressCC_DoubletSinglet.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no")) %>%
  filter(DE_cat != "no") %>%
  mutate(cluster=factor(cluster)) %>%
  group_by(cluster,DE_cat) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(cluster = fct_reorder(cluster,tot,.desc = T)) %>%
  ggplot(aes(x=cluster,y=n,fill=DE_cat)) + geom_col() +
  theme_bw()
ggsave("out/image/DE_barplot_Old_Young.pdf",width = 8,height = 4)

read_tsv("out/table/response_Adult_vs_Young_data.combined_fix_regressCC_DoubletSinglet.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no")) %>%
  filter(DE_cat != "no") %>%
  mutate(cluster=factor(cluster)) %>%
  group_by(cluster,DE_cat) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(cluster = fct_reorder(cluster,tot,.desc = T)) %>%
  ggplot(aes(x=cluster,y=n,fill=DE_cat)) + geom_col() +
  theme_bw()
ggsave("out/image/DE_barplot_Adult_Young.pdf",width = 8,height = 4)


read_tsv("out/table/response_Adult_vs_Young_data.combined_fix_regressCC_DoubletSinglet.tsv") %>%
  filter(gene %in% c("Irf7","Ddx58"))

gene_id <- read_tsv("out/table/response_Old_vs_Young_data.combined_fix_regressCC_DoubletSinglet.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no")) %>%
  filter(DE_cat != "no",
         cell_id == "Endo")

AverageExpression(object = data.combined,features = unique(gene_id$gene),group.by = c("treat","cell_id"))$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  select(gene,contains("Endo")) %>%
  pivot_longer(names_to = "sample",values_to = "avg_exp",-gene) %>%
  separate(sample,into = c("treat","cell_id"),sep = "_",remove = F) %>%
  mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  left_join(gene_id,by = "gene") %>%
  ggplot(aes(x=treat,y=avg_exp,group=gene,col=DE_cat))+geom_line()+facet_wrap(~gene,scales = "free")+theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)) +
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/Line_plot__Old_Adult_Young.pdf",width = 10,height = 10)

data.combined2 <- data.combined
data.combined2$treat <- factor(data.combined2$treat,levels = c("Young","Adult","Old"))

pdf("out/image/Violin_Old_Adult_Young.pdf",width = 12,height = 18)
VlnPlot(data.combined2, features = gene_id$gene,group.by = "treat")
dev.off()

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ComplexHeatmap)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../out_large/Belfast_scRNAseq_ECFC_senescence/object/data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.rds")

# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
str_subset(rownames(data.combined),pattern = c("IRF7|DDX58"))
# "CXCL9" is missing
GOI <- c("IRF7","DDX58")

# generate the table for the plots ----------------------------------------
# get the metadata from the other object
meta <- data.combined@meta.data %>%
  rownames_to_column(var = "barcodes") %>%
  mutate(orig.ident = factor(orig.ident,levels = c("P08","P09","P19","P26")))

# get the normalized expression
df_exp <- data.combined@assays$RNA@data %>%
  data.frame()%>%
  mutate(gene = rownames(.))%>%
  filter(gene%in%GOI)%>%
  gather(key = barcodes,value = count, -gene) %>%
  mutate(barcodes = str_replace(barcodes,pattern = "\\.",replacement = "-")) %>%
  # try to min/max normalize the count varaible per gene in order to rescale the difference in terms of expression
  group_by(gene) %>%
  # threshold of positiveness is based on the distriubtion of the expression of the signal in tihs case
  mutate(norm_min_max = ((count - min(count))/(max(count)-min(count))),
         exp_cat = case_when(count > 0~"pos",
                             T~"neg")) %>%
  ungroup() %>%
  mutate(count_fix = count + rnorm(nrow(.))/100000)

# get the coordinates
UMAP1_df <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column(var = "barcodes")

# generate the dataset for mapping the data in the umamp
dim(UMAP1_df)
dim(df_exp)
dim(meta)

df_tot <- reduce(list(meta,UMAP1_df,df_exp),left_join, by="barcodes")
df_tot_avg <- df_tot %>% group_by(seurat_clusters) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

dim(df_tot)

# plot the average expression per sample use the variable cell tyep per donor as grouping
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_id)
data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$seurat_clusters)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_GOI <- AverageExpression(data.combined,features = GOI)

# plot basic violins ------------------------------------------------------
# fix the order of the metadata
data.combined$orig.ident <- factor(data.combined$orig.ident,levels = c("P08","P09","P19","P26"))

Idents(data.combined) <- "seurat_clusters"

# basic seurat graphic
VlnPlot(data.combined, features = GOI,group.by = "orig.ident")

# use the ggplot graphic
list_plot <- lapply(GOI, function(x){
  test <- VlnPlot(data.combined, features = x, group.by = "orig.ident",split.by = "seurat_clusters", combine = FALSE,slot = "data")
})

# extract the data from the plot 
df_violin <- lapply(list_plot,function(x){
  df <- x[[1]]$data
  
  # extract the name of the gene
  gene <- colnames(df)[1]
  
  df %>%
    mutate(gene = gene) %>%
    setNames(c("exp","orig.ident","seurat_clusters","gene"))
}) %>%
  bind_rows() %>%
  tibble() %>%
  mutate(seurat_clusters = as.character(seurat_clusters))


df_violin %>%
  ggplot(aes(x = orig.ident, y = exp)) +
  geom_violin(scale = "width")+
  #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) +
  geom_jitter(alpha = 0.05) +
  facet_grid(seurat_clusters~gene,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/violin_data.combined_annotated_norm_fix_regressCC_DoubletSinglet_GOI3_split.pdf",width = 8,height = 6)

# general UMAP ------------------------------------------------------------
# build the plot using both info
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=seurat_clusters),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()
# facet_wrap(~infection)
# ggsave("images/UMAP_data.combined_annotated_norm_fix_subset_AC_CCA_celltype.pdf",width = 6,height = 4)

# split the sample
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=seurat_clusters),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw() +
  facet_wrap(~orig.ident)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
  )

# expression distribution -------------------------------------------------
# crop the 0 expressing cells
df_exp %>%
  ggplot(aes(x=count))+geom_histogram()+facet_grid(~gene)+theme_bw()+scale_x_log10()+geom_vline(xintercept = 1,col="red",linetype="dotted")

# keep the 0 expressing cells
df_exp %>%
  ggplot(aes(x=count))+geom_histogram()+facet_wrap(~gene)+theme_bw()+
  # scale_x_log10()+
  geom_vline(xintercept = 2.5,col="red",linetype="dotted")

# library(scales)
# show_col(c("#4662D7FF","#FABA39FF","#7A0403FF"))

# dotplot -----------------------------------------------------------------
average_GOI$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group",values_to = "avg",-gene) %>%
  separate(group,into = c("orig.ident","seurat_clusters"),sep = "\\.",remove = F) %>%
  # mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  # mutate(treat = case_when(treat == "ctrl"~"control",
  #                          T~"chronic active")) %>%
  # mutate(treat = factor(treat,levels = c("control","chronic active"))) %>%
  ggplot(aes(x = seurat_clusters,y = avg,col = orig.ident)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_dodge(width=0.9)) +
  theme_bw() +
  facet_wrap(~gene,scales = "free_y",ncol = 1)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
  ) +
  scale_color_manual(values = c("#2C7BB6","#C6E5DB","#FDC980","#D7191C"))

# show_col(colorRampPalette(rev(brewer.pal(n = 5, name ="RdYlBu")))(4))
# colorRampPalette(rev(brewer.pal(n = 5, name ="RdYlBu")))(4)
# show_col(c("#2C7BB6","#C6E5DB","#FDC980","#D7191C"))

ggsave("out/image/dotplot_data.combined_annotated_norm_fix_regressCC_DoubletSinglet_GOI3_split.pdf",width = 5,height = 6)

# expression of the GOI on UMAP -------------------------------------------
# by counts
df_tot %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = count)) + geom_point(alpha = 0.5,size = 0.2) +
  # facet_grid(gene~NMDA_time) + 
  facet_grid(orig.ident~gene) + 
  theme_bw() + scale_color_gradient(low = "gray",high = "blue")+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
  )
# ggsave("images/UMAP_data.combined_annotated_norm_fix_subset_AC_CCA_SPP1_CD44_C1QB_count.pdf",width = 10,height = 6)

# by min max normalized counts
df_tot %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  # arrange(exp_cat) %>%
  arrange(norm_min_max) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = norm_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
  # facet_grid(gene~NMDA_time) + 
  facet_grid(gene~orig.ident) + 
  theme_bw() + 
  scale_color_gradient(low = "gray",high = "blue") +
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
# ggsave("images/UMAP_data.combined_annotated_norm_fix_subset_AC_CCA_SPP1_CD44_count_scaled.pdf",width = 8,height = 6)

df_tot %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  # arrange(exp_cat) %>%
  arrange(norm_min_max) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = norm_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
  # facet_grid(gene~NMDA_time) + 
  facet_grid(gene~orig.ident) + 
  theme_bw() + 
  viridis::scale_color_viridis(option = "H")+
  # scale_color_gradient(low = "gray",high = "blue") +
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("out/image/UMAP_data.combined_annotated_norm_fix_regressCC_DoubletSinglet_GOI3_split_scaled.pdf",width = 8,height = 6)

# show only the comparison Young vs Old

# percentage of positiveness for the expression ---------------------------
# plot the category. being 0 or non zero per cell
df_tot %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.2) +
  # facet_grid(gene~NMDA_time) +
  facet_grid(gene~orig.ident) +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw() + 
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_color_manual(values = c("gray","blue"))
ggsave("out/image/UMAP_data.combined_annotated_norm_fix_regressCC_DoubletSinglet_GOI3_split_cat.pdf",width = 8,height = 6)

df_tot %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5) +
  # facet_grid(gene~NMDA_time) +
  facet_grid(gene~orig.ident) +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw() + 
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_color_manual(values = c("gray","blue")) +
  theme(panel.border = element_rect(linewidth = 3, fill = NA),
      axis.title.y = element_text(size = 30),
      axis.ticks=element_line(size=4, linetype="solid", colour="black"),
      axis.title.x = element_text(size = 30),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      strip.text = element_text(size = 15))
ggsave("out/image/UMAP_data.combined_annotated_norm_fix_regressCC_DoubletSinglet_GOI3_split_cat_reinhold.pdf",width = 16,height = 8)

# barplot
df_tot %>%
  # group_by(NMDA_time,cell_type) %>%
  group_by(gene,orig.ident,seurat_clusters) %>%
  summarise(cells = n(),
            pos = sum(exp_cat=="pos")) %>%
  group_by(gene,orig.ident) %>%
  mutate(tot_cells = sum(cells)) %>%
  mutate(prop_pos = pos/cells) %>%
  ggplot(aes(x = seurat_clusters,y = prop_pos,col=orig.ident)) +
  # geom_col(position = "dodge") +
  geom_point(position=position_dodge(width=0.9)) +
  facet_wrap(~gene,ncol = 1,scales = "free_y") +
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
  ) +
  scale_color_manual(values = c("#2C7BB6","#C6E5DB","#FDC980","#D7191C"))
ggsave("out/image/dotplot_data.combined_annotated_norm_fix_regressCC_DoubletSinglet_GOI2_split_prop_cluster.pdf",width = 5,height = 6)

# proportion on total sample split by condition
df_tot %>%
  group_by(gene,orig.ident,seurat_clusters) %>%
  summarise(cells = n(),
            pos = sum(exp_cat=="pos")) %>%
  ungroup() %>%
  group_by(gene,orig.ident) %>%
  mutate(tot_cells = sum(cells)) %>%
  mutate(prop_pos = pos/tot_cells) %>%
  # ggsave("images/barplot_data.combined_annotated_norm_fix_subset_AC_CCA_SPP1_CD44_C1QB_expcat_sample.pdf",width = 10,height = 3)
  ggplot(aes(x = seurat_clusters,y = prop_pos,col=orig.ident)) +
  # geom_col(position = "dodge") +
  geom_point(position=position_dodge(width=0.9)) +
  facet_wrap(~gene,ncol = 1,scales = "free_y") +
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
  ) +
  scale_color_manual(values = c("#2C7BB6","#C6E5DB","#FDC980","#D7191C"))
ggsave("out/image/dotplot_data.combined_annotated_norm_fix_regressCC_DoubletSinglet_GOI2_split_prop_tot.pdf",width = 5,height = 6)

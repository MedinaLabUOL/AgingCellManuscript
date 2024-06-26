# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(cowplot)
library(SCPA)
library(circlize)
library(magrittr)
library(msigdbr)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS(file = "../out_large/Belfast_scRNAseq_ECRetinaMouse/object/data.combined_annotated_norm_fix_regressCC_DoubletSinglet.rds")

# generate the table for the plots ----------------------------------------
# get the metadata from the other object
meta <- data.combined@meta.data %>%
  rownames_to_column(var = "barcodes") %>%
  mutate(treat = factor(treat,levels = c("Young","Adult","Old")))

# get the coordinates
UMAP1_df <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column(var = "barcodes")

# plot A ------------------------------------------------------------------
# following Reinhold request I am fixing the annotation:
df_plot_UMAP <- purrr::reduce(list(meta,UMAP1_df),left_join, by="barcodes") %>% 
  mutate(cell_id_reinhold = case_when(cell_id %in% c("MAC_inf","IMM_13")~"MG",
                                      cell_id %in% c("IMM_14")~"IMM",
                                      str_detect(cell_id,pattern = "BP_")~"BP",
                                      T~cell_id))

# 
df_plot_UMAP %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2,col=cell_id_reinhold))+geom_point(size=0.3,alpha=0.6)+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))
ggsave("out/image/plot_A01.pdf",width = 5,height = 4)

# adjust following reinhold request
df_plot_UMAP %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2,col=cell_id_reinhold))+geom_point(size=0.3,alpha=0.6)+theme_cowplot(line_size = 2)+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))+
  theme(strip.background = element_blank(),
        #axis.text.x = element_text(hjust = 1,angle = 45),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"))
ggsave("out/image/plot_A01_reinhold.pdf",width = 5,height = 4)

# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
# Alk1|Alk2|Alk5|ActRIIa|ActRI
# Acvrl1|Acvr1|Tgfbr1|Acvr2a|Acvr1
# str_subset(rownames(data.combined),pattern = c("Acvr|Acvrl1|Acvr1|Tgfbr1|Acvr2a"))
# GOI <- c("Acvr2a","Acvr1","Tgfbr1","Acvr2b","Acvrl1","Acvr1b")
# df_GOI <- read_csv("data/gene_BMP9.csv")
# GOI <- df_GOI %>% 
#   pull(mouse_gene_symbol) %>% 
#   unique()
GOI <- c("Irf7","Isg15","Ifi44")

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

# generate the dataset for mapping the data in the umamp
dim(UMAP1_df)
dim(df_exp)
dim(meta)

df_tot <- purrr::reduce(list(df_plot_UMAP,df_exp),left_join, by="barcodes")
df_tot_avg <- df_tot %>% group_by(cell_id_reinhold) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

dim(df_tot)

# plot the average expression per sample use the variable cell tyep per donor as grouping
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_id)
cell_id_reinhold <- data.combined@meta.data %>%
  mutate(cell_id_reinhold = case_when(cell_id %in% c("MAC_inf","IMM_13")~"MG",
                                      cell_id %in% c("IMM_14")~"IMM",
                                      str_detect(cell_id,pattern = "BP_")~"BP",
                                      T~cell_id)) %>% 
  pull(cell_id_reinhold)

data.combined$group <- paste0(data.combined$treat,".",cell_id_reinhold)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_GOI <- AverageExpression(data.combined,features = GOI)

# plot basic violins ------------------------------------------------------
# fix the order of the metadata
data.combined$treat <- factor(data.combined$treat,levels = c("Young","Adult","Old"))

Idents(data.combined) <- "cell_id"

# basic seurat graphic
VlnPlot(data.combined, features = GOI,group.by = "treat")
ggsave("out/image/plot_B02.pdf",width = 6,height = 4)

# # use the ggplot graphic
# list_plot <- lapply(GOI, function(x){
#   test <- VlnPlot(data.combined, features = x, group.by = "treat",split.by = "cell_id", combine = FALSE,slot = "data")
# })
# 
# # extract the data from the plot 
# df_violin <- lapply(list_plot,function(x){
#   df <- x[[1]]$data
#   
#   # extract the name of the gene
#   gene <- colnames(df)[1]
#   
#   df %>%
#     mutate(gene = gene) %>%
#     setNames(c("exp","treat","cell_id","gene"))
# }) %>%
#   bind_rows() %>%
#   tibble() %>%
#   mutate(cell_id = as.character(cell_id))
# 
# 
# df_violin %>%
#   mutate(group = case_when(str_detect(cell_id,pattern = "BP")~"BP",
#                            T~cell_id)) %>%
#   filter(group %in% c("Endo","BP","Rod")) %>%
#   ggplot(aes(x = treat, y = exp)) +
#   geom_violin(scale = "width")+
#   #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) +
#   geom_jitter(alpha = 0.05) +
#   facet_grid(group~gene,scales = "free") +
#   theme_bw() +
#   theme(strip.background = element_blank(), 
#         panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45))
# # ggsave("out/image/violin_data.combined_annotated_norm_fix_regressCC_DoubletSinglet_GOI2_split.pdf",width = 6,height = 4)

# general UMAP ------------------------------------------------------------
# build the plot using both info
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=cell_id_reinhold),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = cell_id_reinhold),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_cowplot()
ggsave("out/image/plot_A02.pdf",width = 5,height = 4)
# facet_wrap(~infection)
# ggsave("images/UMAP_data.combined_annotated_norm_fix_subset_AC_CCA_celltype.pdf",width = 6,height = 4)

# split the sample
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=cell_id_reinhold),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = cell_id_reinhold),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_cowplot() +
  facet_wrap(~treat)+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
  # )
  theme(strip.background = element_blank())
ggsave("out/image/plot_A03.pdf",width = 10,height = 4)

# adjust following reinhold
ggplot(label= TRUE) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=cell_id_reinhold),size=0.3) +
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))+
  theme_cowplot(line_size = 2)+
  facet_wrap(~treat)+
  theme(strip.background = element_blank(),
        #axis.text.x = element_text(hjust = 1,angle = 45),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"))
ggsave("out/image/plot_A03_reinhold.pdf",width = 10,height = 4)

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
  separate(group,into = c("treat","cell_id"),sep = "\\.",remove = F) %>%
  mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  # mutate(treat = case_when(treat == "ctrl"~"control",
  #                          T~"chronic active")) %>%
  # mutate(treat = factor(treat,levels = c("control","chronic active"))) %>%
  ggplot(aes(x = cell_id,y = avg,col = treat)) +
  # geom_boxplot(outlier.shape = NA) +
  # geom_point(position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.2),alpha=1,size=2)+ theme_bw()+
  geom_point(position=position_dodge(width=0.9),alpha=0.8,size=2)+ theme_bw()+
  facet_wrap(~gene,scales = "free_y",ncol = 1)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
  )+
  scale_color_manual(values = c("#4662D7FF","#FABA39FF","#7A0403FF"))
ggsave("out/image/plot_B03.pdf",width = 5,height = 4)
# ggsave("out/image/dotplot_data.combined_annotated_norm_fix_regressCC_DoubletSinglet_GOI2_split.pdf",width = 6,height = 4)

# expression of the GOI on UMAP -------------------------------------------
# # by counts
# df_tot %>%
#   mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
#   arrange(exp_cat) %>%
#   ggplot(aes(x = UMAP_1, y = UMAP_2,col = count)) + geom_point(alpha = 0.5,size = 0.2) +
#   # facet_grid(gene~NMDA_time) + 
#   # facet_grid(treat~gene) + 
#   facet_grid(~gene) + 
#   theme_bw() + scale_color_gradient(low = "gray",high = "blue")+
#   theme(strip.background = element_blank(), 
#         panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
#   )
# # ggsave("images/UMAP_data.combined_annotated_norm_fix_subset_AC_CCA_SPP1_CD44_C1QB_count.pdf",width = 10,height = 6)

# # by min max normalized counts
# df_tot %>%
#   mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
#   # arrange(exp_cat) %>%
#   arrange(norm_min_max) %>%
#   ggplot(aes(x = UMAP_1, y = UMAP_2,col = norm_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
#   # facet_grid(gene~NMDA_time) + 
#   # facet_grid(gene~treat) + 
#   theme_bw() + 
#   scale_color_gradient(low = "gray",high = "blue") +
#   theme(strip.background = element_blank(), 
#         panel.border = element_rect(colour = "black", fill = NA))
# # ggsave("images/UMAP_data.combined_annotated_norm_fix_subset_AC_CCA_SPP1_CD44_count_scaled.pdf",width = 8,height = 6)

df_tot %>%
  mutate(treat_fix = case_when(treat %in% c("Young")~"03 m",
                               treat %in% c("Adult")~"12 m",
                               treat %in% c("Old")~"23 m")) %>% 
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  # arrange(exp_cat) %>%
  arrange(norm_min_max) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = norm_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
  # facet_grid(gene~NMDA_time) + 
  # facet_grid(gene~treat) + 
  facet_grid(gene~treat_fix) + 
  theme_cowplot() + 
  # viridis::scale_color_viridis(option = "H")+
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
ggsave("out/image/plot_B01.pdf",width = 8.5,height = 5)

# modify reinhold
df_tot %>%
  mutate(treat_fix = case_when(treat %in% c("Young")~"03 m",
                               treat %in% c("Adult")~"12 m",
                               treat %in% c("Old")~"23 m")) %>% 
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  # arrange(exp_cat) %>%
  arrange(norm_min_max) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = norm_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
  # facet_grid(gene~NMDA_time) + 
  # facet_grid(gene~treat) + 
  facet_grid(gene~treat_fix) + 
  theme_cowplot(line_size = 2)+
  viridis::scale_color_viridis(option = "H")+
  # scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank(),
        #axis.text.x = element_text(hjust = 1,angle = 45),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"))
ggsave("out/image/plot_B01_reinhold1.pdf",width = 8.5,height = 7.5)

df_tot %>%
  mutate(treat_fix = case_when(treat %in% c("Young")~"03 m",
                               treat %in% c("Adult")~"12 m",
                               treat %in% c("Old")~"23 m")) %>% 
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  # arrange(exp_cat) %>%
  arrange(norm_min_max) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = norm_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
  # facet_grid(gene~NMDA_time) + 
  # facet_grid(gene~treat) + 
  facet_grid(gene~treat_fix) + 
  theme_cowplot(line_size = 2)+
  # viridis::scale_color_viridis(option = "H")+
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank(),
        #axis.text.x = element_text(hjust = 1,angle = 45),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"))
ggsave("out/image/plot_B01_reinhold2.pdf",width = 8.5,height = 7.5)

# positive
df_tot %>%
  mutate(treat_fix = case_when(treat %in% c("Young")~"03 m",
                               treat %in% c("Adult")~"12 m",
                               treat %in% c("Old")~"23 m")) %>% 
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  # arrange(exp_cat) %>%
  arrange(norm_min_max) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.2) +
  # facet_grid(gene~NMDA_time) + 
  # facet_grid(gene~treat) + 
  facet_grid(gene~treat_fix) + 
  theme_cowplot(line_size = 2)+
  # viridis::scale_color_viridis(option = "H")+
  # scale_color_gradient(low = "gray",high = "blue") +
  scale_color_manual(values = c("gray","blue")) +
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank(),
        #axis.text.x = element_text(hjust = 1,angle = 45),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"))
ggsave("out/image/plot_B01_reinhold1.pdf",width = 8.5,height = 7.5)

# plot degs ---------------------------------------------------------------
gene_id <- read_tsv("out/table/response_Old_vs_Young_data.combined_fix_regressCC_DoubletSinglet.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no")) %>%
  filter(DE_cat != "no",
         cell_id == "Endo")


# cell_id_reinhold <- data.combined@meta.data %>%
#   mutate(cell_id_reinhold = case_when(cell_id %in% c("MAC_inf","IMM_13")~"MG",
#                                       cell_id %in% c("IMM_14")~"IMM",
#                                       str_detect(cell_id,pattern = "BP_")~"BP",
#                                       T~cell_id)) %>% 
#   pull(cell_id_reinhold)

data.combined$group <- paste0(data.combined$treat,"_",cell_id_reinhold)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_DE <- AverageExpression(data.combined,features = unique(gene_id$gene))

average_DE$RNA %>% 
  data.frame() %>%
  rownames_to_column("gene") %>%
  select(gene,contains("Endo")) %>%
  pivot_longer(names_to = "sample",values_to = "avg_exp",-gene) %>%
  separate(sample,into = c("treat","cell_id"),sep = "_",remove = F) %>%
  mutate(treat_fix = case_when(treat %in% c("Young")~"03 m",
                               treat %in% c("Adult")~"12 m",
                               treat %in% c("Old")~"23 m")) %>% 
  mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  left_join(gene_id,by = "gene") %>%
  ggplot(aes(x=treat_fix,y=avg_exp,group=gene,col=DE_cat))+geom_line()+facet_wrap(~gene,scales = "free")+theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)) +
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/Plot_C01_ENDO.pdf",width = 10,height = 10)

average_DE$RNA %>% 
  data.frame() %>%
  rownames_to_column("gene") %>%
  select(gene,contains("Rod")) %>%
  pivot_longer(names_to = "sample",values_to = "avg_exp",-gene) %>%
  separate(sample,into = c("treat","cell_id"),sep = "_",remove = F) %>%
  mutate(treat_fix = case_when(treat %in% c("Young")~"03 m",
                               treat %in% c("Adult")~"12 m",
                               treat %in% c("Old")~"23 m")) %>% 
  mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  left_join(gene_id,by = "gene") %>%
  ggplot(aes(x=treat_fix,y=avg_exp,group=gene,col=DE_cat))+geom_line()+facet_wrap(~gene,scales = "free")+theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)) +
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/Plot_C01_ROD.pdf",width = 10,height = 10)

average_DE$RNA %>% 
  data.frame() %>%
  rownames_to_column("gene") %>%
  select(gene,contains(c("Endo","Rod"))) %>%
  pivot_longer(names_to = "sample",values_to = "avg_exp",-gene) %>%
  separate(sample,into = c("treat","cell_id"),sep = "_",remove = F) %>%
  mutate(treat_fix = case_when(treat %in% c("Young")~"03 m",
                               treat %in% c("Adult")~"12 m",
                               treat %in% c("Old")~"23 m")) %>% 
  mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  # left_join(gene_id,by = c("gene","cell_id")) %>%
  ggplot(aes(x=treat_fix,y=avg_exp,col=cell_id,group=cell_id))+geom_line()+facet_wrap(~gene,scales = "free")+theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)) +
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/Plot_C01.pdf",width = 10,height = 10)

# adjust following reinhold
average_DE$RNA %>% 
  data.frame() %>%
  rownames_to_column("gene") %>%
  dplyr::filter(gene %in% c("Bst2", "Irf7", "Ifit1","Ifi44","Isg15","Vcam1", "Tfrc", "Hspa1a", "Fosb")) %>% 
  select(gene,contains(c("Endo","Rod","Cone","BP"))) %>%
  pivot_longer(names_to = "sample",values_to = "avg_exp",-gene) %>%
  separate(sample,into = c("treat","cell_id"),sep = "_",remove = F) %>%
  mutate(treat_fix = case_when(treat %in% c("Young")~"03 m",
                               treat %in% c("Adult")~"12 m",
                               treat %in% c("Old")~"23 m")) %>% 
  mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  # left_join(gene_id,by = c("gene","cell_id")) %>%
  ggplot(aes(x=treat_fix,y=avg_exp,col=cell_id,group=cell_id))+
  geom_line(linewidth = 2,alpha=0.7)+
  facet_wrap(~gene,scales = "free")+
  theme_bw(base_rect_size = 2)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.ticks = element_line(colour = "black",linewidth = 2)
        #axis.ticks.length = unit(.25, "cm")
        ) +
  scale_color_manual(values = c("green","yellow","blue","red")) +
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/Plot_C01_reinhold.pdf",width = 6,height = 6)

average_DE$RNA %>% 
  data.frame() %>%
  rownames_to_column("gene") %>%
  dplyr::filter(gene %in% c("Bst2", "Irf7", "Ifit1","Ifi44","Isg15","Vcam1", "Tfrc", "Hspa1a", "Fosb")) %>% 
  select(gene,contains(c("Endo","Rod","Cone","Muller"))) %>%
  pivot_longer(names_to = "sample",values_to = "avg_exp",-gene) %>%
  separate(sample,into = c("treat","cell_id"),sep = "_",remove = F) %>%
  mutate(treat_fix = case_when(treat %in% c("Young")~"03 m",
                               treat %in% c("Adult")~"12 m",
                               treat %in% c("Old")~"23 m")) %>% 
  mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  # left_join(gene_id,by = c("gene","cell_id")) %>%
  ggplot(aes(x=treat_fix,y=avg_exp,col=cell_id,group=cell_id))+
  geom_line(linewidth = 2,alpha=0.7)+
  facet_wrap(~gene,scales = "free")+
  theme_bw(base_rect_size = 2)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.ticks = element_line(colour = "black",linewidth = 2)
        #axis.ticks.length = unit(.25, "cm")
  ) +
  scale_color_manual(values = c("green","blue","yellow","red")) +
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/Plot_C01_reinhold2.pdf",width = 6,height = 6)

# plot D ------------------------------------------------------------------
# # read in the dataset
# data.combined <- readRDS(file = "../out_large/Belfast_scRNAseq_ECRetinaMouse/object/data.combined_annotated_norm_fix_regressCC_DoubletSinglet.rds")

cell_id_reinhold <- data.combined@meta.data %>%
  mutate(cell_id_reinhold = case_when(cell_id %in% c("MAC_inf","IMM_13")~"MG",
                                      cell_id %in% c("IMM_14")~"IMM",
                                      str_detect(cell_id,pattern = "BP_")~"BP",
                                      T~cell_id)) %>%
  pull(cell_id_reinhold)


# split the dataset based on the cell id
data.combined$cell_id_reinhold <- cell_id_reinhold
list.data.combined <- SplitObject(data.combined,split.by = "cell_id_reinhold")

# for this dataset skip RPE and IMM_14 beacuse there are too few cells
lapply(list.data.combined,function(x){
  table(x@meta.data$treat)
})

# list.data.combined_filter <- list.data.combined[c("Endo","BP_7","BP_9","Cone","BP_6","Peri","Muller","MAC_inf","Rod","IMM_13","BP_10")]
list.data.combined_filter <- list.data.combined

# run SCPA on tailored signatures -----------------------------------------
# read in the senescence pathways
gs_senescence <- readRDS("data/signature_new/senescence_pathways_mouse.rds") %>% 
  bind_rows() %>%
  # remove the ximerakis
  dplyr::filter(str_detect(Pathway,pattern = "XIMERAKIS",negate = T)) %>% 
  mutate(class = "senescence") %>% 
  # dplyr::filter(Pathway %in% c("EPC_TEAM_SENESCENCE_UP_NO_IF","FRIDMAN_SENESCENCE_UP","senmayo","CellAge_Induces","EPC_TEAM_SENESCENCE_DOWN_NO_IF","CellAge_Inhibits")) %>% 
  split(f = .$Pathway)

# run it on all the files the KEGG sigantures
pmap(list(list.data.combined_filter,names(list.data.combined_filter)),function(x,y){
  # which is the dataset in process
  print(y)
  
  # subset the dataset for the comparison
  control <- seurat_extract(x,meta1 = "treat", value_meta1 = "Young")
  treat <- seurat_extract(x,meta1 = "treat", value_meta1 = "Old")
  
  # print(paste("comparing", i))
  scpa_out <- compare_pathways(list(control, treat), gs_senescence) %>%
    # select(Pathway, qval) %>%
    mutate(cluster = y)
  
  # save the output
  write_tsv(scpa_out,paste0("out/table/subset_SCPA_reinhold_FinalPlot/SCPA_out_fix_normalized_data_regressCC_DoubletSinglet_SENshortlist_Old_",y,".tsv"))
})

# inflammation pathways ---------------------------------------------------
# read in the senescence pathways
gs_inflammation <- readRDS("data/signature_new/inflammation_pathways_mouse.rds") %>% 
  bind_rows() %>% 
  mutate(class = "inflammation") %>% 
  split(f = .$Pathway)

# run it on all the files the KEGG sigantures
pmap(list(list.data.combined_filter,names(list.data.combined_filter)),function(x,y){
  # which is the dataset in process
  print(y)
  
  # subset the dataset for the comparison
  control <- seurat_extract(x,meta1 = "treat", value_meta1 = "Young")
  treat <- seurat_extract(x,meta1 = "treat", value_meta1 = "Old")
  
  # print(paste("comparing", i))
  scpa_out <- compare_pathways(list(control, treat), gs_inflammation) %>%
    # select(Pathway, qval) %>%
    mutate(cluster = y)
  
  # save the output
  write_tsv(scpa_out,paste0("out/table/subset_SCPA_reinhold_FinalPlot/SCPA_out_fix_normalized_data_regressCC_DoubletSinglet_INFshortlist_Old_",y,".tsv"))
})

# save all the tables in one summary --------------------------------------
file <- dir("out/table/subset_SCPA_reinhold_FinalPlot/")%>%
  str_subset(pattern = "_Old_")

list_res_tot <- lapply(file,function(x){
  name <- paste0("out/table/subset_SCPA_reinhold_FinalPlot/",x)
  read_tsv(name)
}) %>%
  setNames(file)

df_summary_tot <- list_res_tot %>%
  bind_rows(.id = "file") %>%
  mutate(annotation = str_extract(file,pattern = "SEN|INF"))

df_summary_tot %>% 
  write_tsv("out/table/summary_SCPA_Old_Reinhold_shortlist_plot.tsv")

# plot all the ones shortlisted by reinhold
df_summary_tot <- read_tsv("out/table/summary_SCPA_Old_Reinhold_shortlist_plot.tsv")

# table(df_summary_tot$Pathway)
wide_df_SEN <- df_summary_tot %>%
  # filter(annotation == "SEN") %>%
  #
  dplyr::filter(Pathway %in% c("MOSERLE_IFNA_RESPONSE","FRIDMAN_SENESCENCE_UP","ISM_SCORE","EPC_TEAM_SENESCENCE_UP_NO_IF","CellAge_Induces","senmayo")) %>% 
  select(Pathway,qval,cluster) %>%
  # mutate(cluster = paste0("cluster_",cluster)) %>%
  pivot_wider(values_from = qval,
              names_from = cluster) %>%
  column_to_rownames("Pathway")

col_hm_SEN <- colorRamp2(colors = c("blue", "white", "red"), breaks = c(0, 1, 2))

ht <- Heatmap(wide_df_SEN,
              col = col_hm_SEN,
              name = "Qval",
              show_row_names = T,
              column_names_gp = gpar(fontsize = 8),
              border = T,
              column_km = 2,
              row_km = 3)

pdf("out/image/Plot_D01.pdf",width = 10,height = 5)
draw(ht, heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 50), "mm"))
dev.off()

# adjust following reinhold
ht_reinhold <- Heatmap(wide_df_SEN,
                       col = col_hm_SEN,
                       name = "Qval",
                       show_row_names = T,
                       column_names_gp = gpar(fontsize = 8),
                       border = T,
                       # column_km = 2,
                       # row_km = 3,
                       show_row_dend = F,
                       rect_gp = gpar(col = "white", lwd = 2))

pdf("out/image/Plot_D01_reinhold.pdf",width = 7,height = 3)
draw(ht_reinhold, heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 50), "mm"))
dev.off()

wide_df_SEN_fix <- wide_df_SEN[c("ISM_SCORE","MOSERLE_IFNA_RESPONSE","FRIDMAN_SENESCENCE_UP","CellAge_Induces","senmayo","EPC_TEAM_SENESCENCE_UP_NO_IF"),]

# adjust following reinhold
ht_reinhold2 <- Heatmap(wide_df_SEN_fix,
                        col = col_hm_SEN,
                        name = "Qval",
                        show_row_names = T,
                        column_names_gp = gpar(fontsize = 8),
                        border = T,
                        # column_km = 2,
                        # row_km = 3,
                        show_row_dend = F,
                        rect_gp = gpar(col = "white", lwd = 2),cluster_rows = F)

pdf("out/image/Plot_D01_reinhold2.pdf",width = 7,height = 3)
draw(ht_reinhold2, heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 50), "mm"))
dev.off()

# # plot INF ----------------------------------------------------------------
# df_summary_tot <- read_tsv("out/table/summary_SCPA_Old_Reinhold_shortlist.tsv")
# 
# wide_df_INF <- df_summary_tot %>%
#   filter(annotation == "INF") %>%
#   select(Pathway,qval,cluster) %>%
#   # mutate(cluster = paste0("cluster_",cluster)) %>%
#   pivot_wider(values_from = qval,
#               names_from = cluster) %>%
#   column_to_rownames("Pathway")
# 
# col_hm_INF <- colorRamp2(colors = c("blue", "white", "red"), breaks = c(0, 1.5, 3))
# 
# ht <- Heatmap(wide_df_INF,
#               col = col_hm_INF,
#               name = "Qval",
#               show_row_names = T,
#               column_names_gp = gpar(fontsize = 8),
#               border = T,
#               column_km = 2,
#               row_km = 3)
# 
# pdf("out/image/SCPA_heatmap_INF_Old_reinhold_shortlist.pdf",width = 10,height = 5)
# draw(ht, heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 50), "mm"))
# dev.off()
# 
# 

# -------------------------------------------------------------------------
# make the dotplot of the same panel for the split object
data.combined <- readRDS(file = "../out_large/Belfast_scRNAseq_ECRetinaMouse/object/data.combined_annotated_norm_fix_regressCC_DoubletSinglet.rds")

cell_id_reinhold <- data.combined@meta.data %>%
  mutate(cell_id_reinhold = case_when(cell_id %in% c("MAC_inf","IMM_13")~"MG",
                                      cell_id %in% c("IMM_14")~"IMM",
                                      str_detect(cell_id,pattern = "BP_")~"BP",
                                      T~cell_id)) %>%
  pull(cell_id_reinhold)

# notice that this is done only on the subset of the young (control) cells
data.combined$cell_id2 <- cell_id_reinhold
#
Idents(data.combined) <- "cell_id2"
DimPlot(data.combined)

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
ggsave("out/image/Dotplot_clusters_annotated_big_panle_integration_norm_fix_long_regressCC_DoubletSinglet4.pdf",width = 16,height = 8)

Idents(data.combined) <- "cell_id2"
test_long2 <- DotPlot(data.combined, features = shortlist_features_list_long[1:7], dot.scale = 8,cluster.idents = T) +
  RotatedAxis()

test_long3 <- DotPlot(data.combined, features = shortlist_features_list_long[8:14], dot.scale = 8,cluster.idents = T) +
  RotatedAxis()

test_long2/test_long3

ggsave("out/image/Dotplot_clusters_annotated_big_panle_integration_norm_fix_long_regressCC_DoubletSinglet5.pdf",width = 16,height = 8)

Idents(data.combined) <- "cell_id2"
test_long22 <- DotPlot(data.combined, features = shortlist_features_list_long[c("ENDO","ROD","MULLER")], dot.scale = 8,cluster.idents = T) +
  RotatedAxis()

df_test22 <- lapply(shortlist_features_list_long[c("ENDO","ROD","MULLER")],function(x){
  test_long22$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

plot22 <- df_test22 %>%
  mutate(cell_type = factor(cell_type,c("ENDO","ROD","MULLER"))) %>%
  mutate(id = factor(id,c("MG","BP","IMM","Cone","RPE","Peri","Muller","Rod","Endo"))) %>%
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_size(range = c(0, 6),breaks = c(0,25,50,75,100),limits = c(0,100)) +
  # facet_grid(~cell_type,scales = "free",space="free")+
  facet_wrap(~cell_type,scales = "free_x",nrow=1)+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")

test_long32 <- DotPlot(data.combined, features = shortlist_features_list_long[c("BP","PERiCYTE","CONE")], dot.scale = 8,cluster.idents = T) +
  RotatedAxis()

df_test32 <- lapply(shortlist_features_list_long[c("BP","PERiCYTE","CONE")],function(x){
  test_long32$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

plot32 <- df_test32 %>%
  mutate(cell_type = factor(cell_type,c("BP","PERiCYTE","CONE"))) %>% 
  mutate(id = factor(id,c("MG","BP","IMM","Cone","RPE","Peri","Muller","Rod","Endo"))) %>%
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_size(range = c(0, 6),breaks = c(0,25,50,75,100),limits = c(0,100)) +
  # facet_grid(~cell_type,scales = "free",space="free")+
  facet_wrap(~cell_type,scales = "free_x",nrow=1)+
  theme_cowplot()+
  # theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90),legend.position = "null")+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")

plot22/plot32

ggsave("out/image/Dotplot_clusters_annotated_big_panle_integration_norm_fix_long_regressCC_DoubletSinglet6.pdf",width = 10,height = 7)



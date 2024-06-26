# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(limma)

# aim ---------------------------------------------------------------------
# count the immune cells CD83 or RUNX1 pos

# read in the data --------------------------------------------------------
# read in the annotated object
data.combined <- readRDS("../out_large/Belfast_scRNAseq_ECFC_senescence/object/data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.rds")

# compare the clusters per treatment --------------------------------------
# make suer the correct defalult dataset il loaded should be RNA
DefaultAssay(data.combined)

# check that in the RNA slot the data object is indeed loaded with normalized values
data.combined@assays$RNA@data[1:10,1:10]

# define the grouping variables for the comparison of CSF vs CTRL
head(data.combined@meta.data)

# define the new grouping
# data.combined$infection.treat.cluster <- paste(data.combined$infection, data.combined$treat,data.combined$seurat_clusters, sep = "_")
# head(data.combined@meta.data)

# update the idents of the object
Idents(data.combined) <- "orig.ident"

# avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
# run the differential expression over the WT
clusters_id <- as.character(sort(unique(data.combined$orig.ident))) %>%
  str_subset(pattern = "P09",negate = T)

# Idents(data.combined)
# for the WT the cluster 8 does not exist, therefore I will cropp it out
# x<-0
# id_1 <- paste0("WT_CSF_",x)
# id_2 <- paste0("WT_CTRL_",x)
# response <- FindMarkers(data.combined, ident.1 = id_1, ident.2 = id_2, verbose = FALSE,logfc.threshold = 0)
# response %>%
#   rownames_to_column("gene") %>%
#   mutate(id_1 = id_1,
#          id_2 = id_2) %>%
#   mutate(cluster = x)


list_test_vs_young <- lapply(clusters_id,function(x){
  id_1 <- paste0(x)
  id_2 <- paste0("P09")
  response <- FindMarkers(data.combined, ident.1 = id_1, ident.2 = id_2, verbose = T,logfc.threshold = 0)
  response %>%
    rownames_to_column("gene") %>%
    mutate(id_1 = id_1,
           id_2 = id_2) %>%
    mutate(cluster = x)
})

list_test_vs_young %>%
  setNames(paste0(clusters_id,"_vs_P09")) %>%
  bind_rows(.id = "cluster") %>%
  write_tsv("out/table/response_treat_vs_P09_regressCC_DoubletSinglet.tsv")

# reinhold also wanted to pull the genes for the comparisoin P26 vs P08
response <- FindMarkers(data.combined, ident.1 = "P26", ident.2 = "P08", verbose = T,logfc.threshold = 0)
response %>%
  rownames_to_column("gene") %>%
  mutate(id_1 = "P26",
         id_2 = "P08") %>%
  mutate(cluster = "P26_vs_P08") %>%
  write_tsv("out/table/response_treat_P26_vs_P08_regressCC_DoubletSinglet.tsv")

# -------------------------------------------------------------------------
#
test <- read_tsv("out/table/response_treat_vs_P09_regressCC_DoubletSinglet.tsv")
# check the expression of some of the genes on interes
GOI <- c("IRF7","DDX58","CD34","PROCR")
test %>%
  filter(gene %in% GOI)

# -------------------------------------------------------------------------
# plot the volcanos per cluster
test_plot <- read_tsv("out/table/response_treat_vs_P09_regressCC_DoubletSinglet.tsv")

# test_plot <- list_test_vs_young %>%
#   setNames(paste0(clusters_id,"_vs_P09")) %>%
#   bind_rows(.id = "cluster")

# test_plot %>%
#   ggplot(aes(x=avg_log2FC))+geom_histogram()

# show the distribution of the pvalues
test_plot %>%
  ggplot(aes(x = p_val)) + geom_histogram()+theme_bw()+facet_wrap(~cluster)
ggsave("out/image/dist_p_value_regressCC_DoubletSinglet.pdf",width = 6,height = 3)

test %>%
  ggplot(aes(x = p_val)) + geom_histogram()+theme_bw()+facet_wrap(~cluster)

# render all of them as a volcano plot
test_significant <- test_plot %>%
  mutate(threshold = case_when(abs(avg_log2FC) > 1 & p_val_adj<0.05~1,
                               T~0)) %>%
  filter(threshold == 1)

library(ggrepel)
test %>%
  # filter(symbol %in% setdiff(GOI_SLC,GOI)) %>%
  ggplot(aes(x = avg_log2FC,y = -log(p_val_adj))) + 
  geom_point(alpha = 0.01) +
  geom_point(data = test_significant,aes(x = avg_log2FC,y = -log(p_val_adj)),col="red",alpha = 0.5) +
  geom_vline(xintercept = c(-1,1),col="gray",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="gray",linetype="dashed")+
  geom_text_repel(data = test_significant,aes(x = avg_log2FC,y = -log(p_val_adj),label = gene)) +
  facet_wrap(~cluster) +
  theme_bw()
ggsave("out/image/volcano_regressCC_DoubletSinglet.pdf",width = 12,height = 6)

# -------------------------------------------------------------------------
#
# check the expresison of the panel alos for P26 vs P08 comparison
test <- read_tsv("out/table/response_treat_P26_vs_P08_regressCC_DoubletSinglet.tsv")
# check the expression of some of the genes on interes
GOI <- c("IRF7","DDX58","CD34")
test %>%
  filter(gene %in% GOI)

# -------------------------------------------------------------------------
# plot the volcanos per cluster
test_plot <- read_tsv("out/table/response_treat_P26_vs_P08_regressCC_DoubletSinglet.tsv")

# test_plot %>%
#   ggplot(aes(x=avg_log2FC))+geom_histogram()

# show the distribution of the pvalues
test_plot %>%
  ggplot(aes(x = p_val)) + geom_histogram()+theme_bw()+facet_wrap(~cluster)
ggsave("out/image/dist_p_value_P26_vs_P08_regressCC_DoubletSinglet.pdf",width = 3,height = 3)

# render all of them as a volcano plot
test_significant <- test_plot %>%
  mutate(threshold = case_when(abs(avg_log2FC) > 1 & p_val_adj<0.05~1,
                               T~0)) %>%
  filter(threshold == 1)

library(ggrepel)
test %>%
  # filter(symbol %in% setdiff(GOI_SLC,GOI)) %>%
  ggplot(aes(x = avg_log2FC,y = -log(p_val_adj))) + 
  geom_point(alpha = 0.01) +
  geom_point(data = test_significant,aes(x = avg_log2FC,y = -log(p_val_adj)),col="red",alpha = 0.5) +
  geom_vline(xintercept = c(-1,1),col="gray",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="gray",linetype="dashed")+
  geom_text_repel(data = test_significant,aes(x = avg_log2FC,y = -log(p_val_adj),label = gene)) +
  facet_wrap(~cluster) +
  theme_bw()
ggsave("out/image/volcano_P26_vs_08_regressCC_DoubletSinglet.pdf",width = 4,height = 6)

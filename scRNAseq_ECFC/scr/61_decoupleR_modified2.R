# https://www.bioconductor.org/packages/release/bioc/vignettes/decoupleR/inst/doc/tf_sc.html
# Transcription factor activity inference from scRNA-seq
# Pau Badia-i-Mompel1

# scRNA-seq yield many molecular readouts that are hard to interpret by themselves. One way of summarizing this information is by inferring transcription factor (TF) activities from prior knowledge.

# In this notebook we showcase how to use decoupleR for transcription factor activity inference with a down-sampled PBMCs 10X data-set. The data consists of 160 PBMCs from a Healthy Donor. The original data is freely available from 10x Genomics here from this webpage.
# https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

# 1Loading packages -------------------------------------------------------
# First, we need to load the relevant packages, Seurat to handle scRNA-seq data and decoupleR to use statistical methods.

## We load the required packages
library(Seurat)
library(decoupleR)
library(tidyverse)
library(cowplot)

# # Only needed for data handling and plotting
# library(dplyr)
# library(tibble)
# library(tidyr)
# library(patchwork)
# library(ggplot2)
# library(pheatmap)
library(ComplexHeatmap)

# 2Loading the data-set ---------------------------------------------------
# read in the data
data <- readRDS("../out_large/Belfast_scRNAseq_ECFC_senescence/object/data.combined_NOT_annotated_norm_fix_regressCC_DoubletSinglet.rds")

# We can observe that we have different cell types:
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# 3DoRothEA network -------------------------------------------------------
# DoRothEA is a comprehensive resource containing a curated collection of TFs and their transcriptional targets. Since these regulons were gathered from different types of evidence, interactions in DoRothEA are classified in different confidence levels, ranging from A (highest confidence) to D (lowest confidence). 
# Moreover, each interaction is weighted by its confidence level and the sign of its mode of regulation (activation or inhibition). The direction is coded in the mor column of the dataset.

# For this example we will use the human version (mouse is also available) and we will use the confidence levels ABC. We can use decoupleR to retrieve it from OmniPath:
net <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))
net

# is IRF7 in the dataset?
net %>% 
  dplyr::filter(source == "IRF7")

# add the curated list from reinhold
read_tsv("data/IRF7_targets_partial_babelgene_manual.txt")
  
net2 <- read_tsv("data/TFLink_targets_of_IRF7.txt") %>% 
  pull(Name.Target) %>% 
  str_split(pattern = ";") %>% 
  unlist() %>% 
  unique() %>% 
  data.frame(human_symbol=.) %>% 
  dplyr::filter(!is.na(human_symbol)) %>% 
  mutate(source = "IRF7",
         confidence = "D",
         target = human_symbol,
         mor = 1) %>% 
  # shape it as the net object
  dplyr::select(source,
                confidence,
                target,
                mor) %>% 
  # some genes are alteady present in net, remove them
  dplyr::filter(!target %in% c("CXCL10","IFNB1","IFNA13","IFNA2"))

# # confirm there are no duplicated
# read_tsv("data/IRF7_targets_partial_babelgene_manual.txt") %>% 
#   dplyr::filter(!is.na(human_symbol)) %>% 
#   mutate(source = "IRF7",
#          confidence = "D",
#          target = human_symbol,
#          mor = 1) %>% 
#   # shape it as the net object
#   dplyr::select(source,
#                 confidence,
#                 target,
#                 mor) %>% 
#   pull(target) %>% table()

# add this novel information to the final table
net_final <- bind_rows(list(net,net2))

# confirm the addition
net_final %>% 
  dplyr::filter(source == "IRF7") %>% 
  group_by(target) %>%
  summarise(n = n()) %>% arrange(desc(n))

# the number of targets is small, make sure to change the threshold for the analysis and include it
net_final %>% 
  dplyr::filter(source == "IRF7") %>% 
  mutate(present = target %in% rownames(data@assays$RNA@data)) %>% 
  summarise(present_tot = sum(present))

# 18 targets genes are present in the dataset

# 4Activity inference with Weighted Mean ----------------------------------
# To infer activities we will run the Weighted Mean method (wmean). It infers regulator activities by first multiplying each target feature by its associated weight which then are summed to an enrichment score wmean. Furthermore, permutations of random target features can be performed to obtain a null distribution that can be used to compute a z-score norm_wmean, or a corrected estimate corr_wmean by multiplying wmean by the minus log10 of the obtained empirical p-value.

# In this example we use wmean but we could have used any other. To see what methods are available use show_methods().

# To run decoupleR methods, we need an input matrix (mat), an input prior knowledge network/resource (net), and the name of the columns of net that we want to use.

# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA@data)
mat[1:10,1:10]

# Run wmean
acts <- run_wmean(mat=mat, net=net_final, .source='source', .target='target',
                  .mor='mor', times = 100, minsize = 3)
acts
# save the table of of the results
acts %>% 
  write_tsv("out/table/decoupleR_NOT_annotated_norm_fix_regressCC_DoubletSinglet_modified2.tsv")

# confrim IRF7 made the cut
acts %>% 
  dplyr::filter(str_detect(source,pattern = "IRF7"))

# robust TF signatures ----------------------------------------------------
# read in the result from run_wmean
acts <- read_tsv("out/table/decoupleR_NOT_annotated_norm_fix_regressCC_DoubletSinglet_modified2.tsv")

# Extract norm_wmean and store it in tfswmean in pbmc
data[['tfswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "tfswmean"

# Scale the data
data <- ScaleData(data)
data@assays$tfswmean@data <- data@assays$tfswmean@scale.data

# unbiased exploration of the main TFs ------------------------------------
# We can also see what is the mean activity per group of the top 20 more variable TFs:
n_tfs <- 25

# get the meta of the dataset
LUT_sample <- data@meta.data %>% 
  data.frame() %>% 
  rownames_to_column("barcodes")

# Extract activities from object as a long dataframe
df_acts <- t(as.matrix(data@assays$tfswmean@data)) %>%
  as.data.frame() %>%
  rownames_to_column("barcodes") %>% 
  # mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -barcodes, names_to = "source", values_to = "score") %>%
  # add all the meta
  left_join(LUT_sample,"barcodes")

# summarise per cluster
df_acts_summ_cluster <- df_acts %>%
  group_by(seurat_clusters, source) %>%
  summarise(mean = mean(score))

# summarise per condition
df_acts_summ_cond <- df_acts %>%
  group_by(orig.ident, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs_cluster <- df_acts_summ_cluster %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Get top tfs with more variable means across condition
tfs_cond <- df_acts_summ_cond %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix cluster
top_acts_mat_cluster <- df_acts_summ_cluster %>%
  dplyr::filter(source %in% tfs_cluster) %>%
  pivot_wider(id_cols = 'seurat_clusters', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('seurat_clusters') %>%
  as.matrix()

# Subset long data frame to top tfs and transform to wide matrix sample
top_acts_mat_cond <- df_acts_summ_cond %>%
  dplyr::filter(source %in% tfs_cond) %>%
  pivot_wider(id_cols = 'orig.ident', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('orig.ident') %>%
  as.matrix()

# Choose color palette
palette_length <- 100
my_color <- colorRampPalette(c("Darkblue", "white","red"))(palette_length)

# my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
#                seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
# pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
# implement it in complexheatmap
# plot TF per cluster
pdf("out/image/decoupleR_NOT_annotated_norm_fix_regressCC_DoubletSinglet_topTF_cluster_modified2.pdf",width = 6,height = 3)
Heatmap(top_acts_mat_cluster)
dev.off()

# plot TF per sample
pdf("out/image/decoupleR_NOT_annotated_norm_fix_regressCC_DoubletSinglet_topTF_condition_modified2.pdf",width = 6,height = 3)
Heatmap(top_acts_mat_cond)
dev.off()

# visualization -----------------------------------------------------------
# confirm the activity of specific factors
# from act chsck MYC
DefaultAssay(object = data)

p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('clusters')
p2 <- (FeaturePlot(data, features = c("IRF7")) & scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('IRF7 activity')

DefaultAssay(object = data) <- "RNA"
p3 <- FeaturePlot(data, features = c("IRF7")) + ggtitle('IRF7 expression')
DefaultAssay(object = data) <- "tfswmean"
p1 | p2 | p3
ggsave("out/image/decoupleR_robust_NOT_annotated_norm_fix_regressCC_DoubletSinglet_IRF7_activity_modified2.pdf",width = 10,height = 3)

# plot the activity score for the sample activity
df_UMAP <- data@reductions$umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column("barcodes")

df_meta <- data@meta.data %>% 
  data.frame() %>% 
  rownames_to_column("barcodes")

df_activity <- data@assays$tfswmean@data %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(names_to = "barcodes",values_to = "activity",-gene) %>% 
  mutate(barcodes = str_replace(barcodes,pattern = "\\.",replacement = "-")) %>% 
  filter(gene %in% c("IRF7"))

# join the three dataset
df_all <- purrr::reduce(list(df_meta,df_UMAP,df_activity),left_join,by="barcodes")
df_all %>% 
  arrange(desc(activity)) %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2,col=activity)) +
  geom_point(size = 0.5) +
  # scale_color_viridis_c(option = "viridis") +
  scale_color_gradient2(low = 'blue',
                        high = 'red',
                        mid = 'white',
                        midpoint = 0
                        # limit = c(-1, 1)
  ) +
  theme_cowplot() +
  facet_wrap(~orig.ident,nrow = 1) +
  theme(strip.background = element_blank())
ggsave("out/image/decoupleR_robust_NOT_annotated_norm_fix_regressCC_DoubletSinglet_IRF7_activity_split_modified2.pdf",width = 10,height = 3)

# # visualize specific TFs activity -----------------------------------------
# # IRF7
# DefaultAssay(object = data)
# 
# p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
#   NoLegend() + ggtitle('clusters')
# p2 <- (FeaturePlot(data, features = c("IRF7")) & scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
#   ggtitle('IRF7 activity')
# 
# DefaultAssay(object = data) <- "RNA"
# p3 <- FeaturePlot(data, features = c("IRF7")) + ggtitle('IRF7 expression')
# DefaultAssay(object = data) <- "tfswmean"
# p1 | p2 | p3
# ggsave("out/image/decoupleR_OneGene_NOT_annotated_norm_fix_regressCC_DoubletSinglet_IRF7_activity.pdf",width = 10,height = 3)

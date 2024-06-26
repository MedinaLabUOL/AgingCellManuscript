# this script is run to explore the dataset individually and define a consensus threshold for the mito and reads counts
# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)
library(DoubletFinder)

# READ IN DATA ------------------------------------------------------------
# file <- c("ECFC_P09")
id_sample <- dir("../data/cellranger/") %>%
  str_subset(pattern = "PN")

# load the LUT
LUT <- read_csv("data/LUT_sample.csv")

# do the preprocessing over all the dataset and save the objects
# /media/edo/INTENSO/RNAseq/scRNASeq/raw_data/Belfast_scRNAseq_ECFC_senescence/cellranger
# x <- "PN0371_0001"
list_datasc <- lapply(id_sample,function(x){
  data <- Read10X(data.dir = paste0("../data/cellranger/",x,"/outs/filtered_feature_bc_matrix/"))
  
  datasc <- CreateSeuratObject(counts = data, project = LUT %>%
                                 filter(sample_id == x) %>%
                                 pull(sample), min.cells = 20, min.features = 200)
  
  # datasc <- CreateSeuratObject(counts = data, min.cells = 20, min.features = 200)
  
  # datasc@meta.data
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^mt-")
  
  # label the cells based on the mt reads content
  datasc$mt_bin <- datasc@meta.data %>%
    mutate(test = case_when(percent.mt < 1~"low",
                            percent.mt < 20~"mid",
                            T ~ "high")) %>%
    pull(test)
  
  # add the meta for the treatmnett
  datasc$treat <- LUT %>%
    filter(sample_id == x) %>%
    pull(treat)
  
  # add the filtering variable based on the fixed threshold
  datasc$test <- datasc@meta.data %>%
    # mutate(test = percent.mt < 20 & nFeature_RNA > 700 & nFeature_RNA < 9000) %>%
    # mutate(test = percent.mt > 1 & percent.mt < 10 & nFeature_RNA > 1000 & nFeature_RNA < 9000) %>%
    mutate(test = percent.mt > 1 & percent.mt < 20 & nFeature_RNA > 500 & nFeature_RNA < 6000) %>%
    pull(test)
  
  # add the filtering variable based on the
  stats <- cbind(log10(datasc@meta.data$nCount_RNA), log10(datasc@meta.data$nFeature_RNA),
                 datasc@meta.data$nFeature_RNA)
  
  # library(robustbase)
  outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
  #library(scater)
  multi.outlier <- isOutlier(outlying, type = "higher")
  # summary(multi.outlier)
  
  datasc$not_outlier <- !as.vector(multi.outlier)
  
  datasc
}) %>%
  setNames(id_sample)

# plot QC -----------------------------------------------------------------
# run and plot QC
# extract the metadata from each dataset
meta_total <- lapply(list_datasc, function(x){
  x@meta.data
}) %>%
  bind_rows(.id = "dataset") %>%
  rownames_to_column("barcode")

meta_total %>%
  write_tsv("out/table/meta_datasc_fix_filter_norm_total.tsv")

# how many cells are considered outliers
meta_total %>%
  dplyr::count(orig.ident,test)

meta_total %>%
  dplyr::count(orig.ident,not_outlier)

# fixed threshold scatter nFeature vs percent.mt
meta_total %>%
  ggplot(aes(y = percent.mt,x = nFeature_RNA,col=test)) + geom_point(alpha=0.3)+ facet_wrap(~orig.ident) + 
  theme_bw() +
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
# save the plot
ggsave("out/image/fixed_scatter_feature_mito.pdf",width = 9,height = 3)

meta_total %>%
  ggplot(aes(y = percent.mt,x = nFeature_RNA,col=not_outlier)) + geom_point(alpha=0.3)+ facet_wrap(~orig.ident) +
  theme_bw() +
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
# save the plot
ggsave("out/image/adaptive_scatter_feature_mito.pdf",width = 9,height = 3)

#
meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt)) %>%
  ggplot(aes(x=orig.ident,y=value))+geom_violin()+geom_jitter(width = 0.2,alpha=0.01) + facet_wrap(~var,scales = "free") + theme_bw() + 
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45,hjust = 1))
# save the plot
ggsave("out/image/fixed_boxplot_reads.pdf",width = 9,height = 3)

#
meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt)) %>%
  filter(var == "percent.mt") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) + facet_wrap(orig.ident~var,scales = "free")+theme_bw()+scale_x_log10() + geom_vline(xintercept = c(1,20),col="red")+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
# save the plot
ggsave("out/image/fixed_histo_mito.pdf",width = 9,height = 3)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt)) %>%
  filter(var == "nFeature_RNA") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) + facet_wrap(orig.ident~var,scales = "free")+theme_bw()+scale_x_log10() + geom_vline(xintercept = c(500,6000),col="red")+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
# save the plot
ggsave("out/image/fixed_histo_features.pdf",width = 9,height = 3)

#
meta_total %>%
  dplyr::count(mt_bin,orig.ident)
# we are going to select only the test T
meta_total %>%
  dplyr::count(orig.ident,mt_bin,test)

# color the bins for the amount of reads
meta_total %>%
  ggplot(aes(x = nCount_RNA,y = nFeature_RNA,col=mt_bin)) + geom_point(alpha=0.3) + facet_grid(orig.ident~mt_bin,scales = "free_y")+theme_bw() + 
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45,hjust = 1))
# save the plot
ggsave("out/image/fixed_scatter_mito.pdf",width = 10,height = 9)

# filtering of the cells --------------------------------------------------
# perform the filtering based on the fixed threshold
list_datasc_fixed <- lapply(list_datasc,function(x){
  datasc_filter <- subset(x, subset = test == 1)
})

# Pre-process the data for DoubletFinder ----------------------------------
# use the Normalize function. for the doublet removla I need to run the full preprocessing steps
list_datasc_fixed_norm <- lapply(list_datasc_fixed, function(x){
  datasc_filter <- x
  datasc_filter <- NormalizeData(datasc_filter)
  datasc_filter <- FindVariableFeatures(datasc_filter, selection.method = "vst", nfeatures = 2000)
  datasc_filter <- ScaleData(datasc_filter) 
  datasc_filter <- RunPCA(datasc_filter) 
  datasc_filter <- FindNeighbors(datasc_filter,dims = 1:30) 
  datasc_filter <- FindClusters(datasc_filter) 
  datasc_filter <- RunUMAP(datasc_filter,dims = 1:30) 
  
  datasc_filter
})

# DoubletFinder -----------------------------------------------------------
# determine the number of cell recovered
df_doublet <- read_csv("data/doublets.csv")

list_nExp <- lapply(list_datasc_fixed_norm,function(x){
  recover_cell <- round(dim(x)[2]/1000,digits = 0)*1000
  recover_cell
  ifelse(test = recover_cell > 10000,
         yes = 0.076,
         no = df_doublet[df_doublet$CellRecovered == recover_cell,][["MultipletRate"]]
  )
})

# runt he simulation for doublet finder
list_datasc_fixed_norm_doublet <- pmap(list(list_datasc_fixed_norm,list_nExp),function(x,y){
  
  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep_v3(x, PCs = 1:30, sct = FALSE) 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats) 
  # plot the results to justify the pK choice 
  # bcmvn %>% 
  #   ggplot(aes(pK,BCmetric,group=1)) + 
  #   geom_point()+ 
  #   geom_line() 
  # save the max BCmetrics 
  pK <- bcmvn %>% 
    filter(BCmetric == max(BCmetric)) %>% 
    pull(pK) %>% 
    as.character() %>% 
    as.numeric() 
  
  # homotipic doublet proportion estimate ----------------------------------- 
  annotation <- x$seurat_clusters 
  homotipic.prop <- modelHomotypic(annotation) 
  
  # 0.076 is the expected nbumber of doublets based on the number of cells recovered 
  nExp.poi <- round(y*nrow(x@meta.data)) 
  nExp.poi.adj <- round(nExp.poi*(1-homotipic.prop)) 
  
  # run the DoubletFinder --------------------------------------------------- 
  pbmc.seurat.filteres <- doubletFinder_v3(x,PCs = 1:30,pN = 0.25,pK = pK,nExp = nExp.poi.adj,reuse.pANN = F,sct = F)
  pbmc.seurat.filteres
})

# list_datasc_fixed_norm_doublet$PN0372_0001@meta.data
# c("orig.ident","nCount_RNA", "nFeature_RNA", "percent.mt", "mt_bin", "treat" "test", "not_outlier", "RNA_snn_res.0.8", "seurat_clusters", "pANN", "DF")
#
lapply(list_datasc_fixed_norm_doublet, function(x){
  meta <- x@meta.data
  colnames(meta) <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","mt_bin","treat","test","not_outlier","RNA_snn_res.0.8","seurat_clusters","pANN","DF")
  x@meta.data <- meta
  head(x@meta.data)
  dim(x@meta.data)[1]
})

# save the individula filtered and normzlized objects
pmap(list(names(list_datasc_fixed_norm_doublet),list_datasc_fixed_norm_doublet),function(x,y){
  # fix the meta and save the object
  meta <- y@meta.data
  colnames(meta) <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","mt_bin","treat","test","not_outlier","RNA_snn_res.0.8","seurat_clusters","pANN","DF")
  y@meta.data <- meta
  # head(y@meta.data)
  #
  saveRDS(object = y,file = paste0("../out_large/Belfast_scRNAseq_ECRetinaMouse/object/datasc_fix_filter_norm_doublet_",x,".rds"))
})

# save the taola meta with the doublet imputation
meta_total_doublet <- lapply(list_datasc_fixed_norm_doublet, function(x){
  meta <- x@meta.data
  colnames(meta) <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","mt_bin","treat","test","not_outlier","RNA_snn_res.0.8","seurat_clusters","pANN","DF")
  x@meta.data <- meta
  # head(x@meta.data)
  x@meta.data
}) %>%
  bind_rows(.id = "dataset") %>%
  rownames_to_column("barcode")

# save the total meta
meta_total_doublet %>%
  write_tsv("out/table/meta_datasc_fix_filter_norm_doublet.tsv")

# subset the datasets after removal of the doublets befreo the integration
pmap(list(names(list_datasc_fixed_norm_doublet),list_datasc_fixed_norm_doublet),function(x,y){
  # fix the meta
  meta <- y@meta.data
  colnames(meta) <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","mt_bin","treat","test","not_outlier","RNA_snn_res.0.8","seurat_clusters","pANN","DF")
  y@meta.data <- meta
  # head(y@meta.data)
  # filter only the singlets
  datasc_filter <- subset(y, subset = DF == "Singlet")
  saveRDS(object = datasc_filter,file = paste0("../out_large/Belfast_scRNAseq_ECRetinaMouse/object/datasc_fix_filter_norm_doublet_Singlet_",x,".rds"))
})

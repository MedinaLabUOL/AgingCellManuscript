# libraries ---------------------------------------------------------------
library("tidyverse")
library("DESeq2")
library("limma")
library("AnnotationDbi") 
library("AnnotationHub")
library("ashr")

# read in the data --------------------------------------------------------
# save the filtered object
dds2_filtered <- readRDS("out/object/dds2_filtered.rds")
design <- readRDS("out/object/design.rds")

# print the contrast
resultsNames(dds2_filtered)

# save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast <- makeContrasts(RIGIvsNORIP = RIP_RIGI,
                          SNRP70vsNORIP = RIP_SNRP70,
                          RIGIvsSNRP70 = RIP_RIGI-RIP_SNRP70,
                          levels = design)

res_RIGIvsNORIP <- results(dds2_filtered, contrast=contrast[,"RIGIvsNORIP"],alpha = 0.05)
res_RIGIvsNORIP %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  write_tsv("out/table/res_RIGIvsNORIP.tsv")
res_RIGIvsNORIP_shr <- lfcShrink(dds2_filtered, res = res_RIGIvsNORIP, type = "ashr")
res_RIGIvsNORIP_shr %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  write_tsv("out/table/res_RIGIvsNORIP_shr.tsv")

summary(res_RIGIvsNORIP)
summary(res_RIGIvsNORIP_shr)

res_RIGIvsSNRP70 <- results(dds2_filtered, contrast=contrast[,"RIGIvsSNRP70"],alpha = 0.05)
res_RIGIvsSNRP70 %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  write_tsv("out/table/res_RIGIvsSNRP70.tsv")
res_RIGIvsSNRP70_shr <- lfcShrink(dds2_filtered, res = res_RIGIvsSNRP70, type = "ashr")
res_RIGIvsSNRP70_shr %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  write_tsv("out/table/res_RIGIvsSNRP70_shr.tsv")
summary(res_RIGIvsSNRP70)
summary(res_RIGIvsSNRP70_shr)

res_SNRP70vsNORIP <- results(dds2_filtered, contrast=contrast[,"SNRP70vsNORIP"],alpha = 0.05)
res_SNRP70vsNORIP %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  write_tsv("out/table/res_SNRP70vsNORIP.tsv")
res_SNRP70vsNORIP_shr <- lfcShrink(dds2_filtered, res = res_SNRP70vsNORIP, type = "ashr")
res_SNRP70vsNORIP_shr %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  write_tsv("out/table/res_SNRP70vsNORIP_shr.tsv")
summary(res_SNRP70vsNORIP)
summary(res_SNRP70vsNORIP_shr)

# check some trends of interest
res_SNRP70vsNORIP %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  arrange(padj) %>% dplyr::filter(str_detect(gene,pattern = "SNR"))

res_RIGIvsSNRP70 %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  arrange(padj) %>% dplyr::filter(str_detect(gene,pattern = "SNR"))

# add annotation ----------------------------------------------------------
# read in the annotation fo the genes
df_LUT <- readRDS("out/table/df_LUT_gtf_genes.tsv")
# how many genes are ribosomal in the dataset
df_LUT %>% 
  dplyr::filter(biotype=="rRNA")

# loop all the non corrected files
folder <- "out/table/"
file_raw <- dir(folder) %>% 
  str_subset("^res_") %>% 
  str_subset("_shr.tsv",negate = T)

file_shr <- dir(folder) %>% 
  str_subset("^res_") %>% 
  str_subset("_shr.tsv",negate = F)
  
list_raw <- lapply(file_raw, function(x){
  df <- read_tsv(paste0(folder,x)) %>% 
    left_join(df_LUT,by = c("gene"="gene_ids")) %>% 
    arrange(padj)
  
  return(df)
}) %>% 
  setNames(str_remove_all(file_raw,pattern = "res_|.tsv"))
saveRDS(list_raw,"out/object/res_raw.rds")

list_shr <- lapply(file_shr, function(x){
  df <- read_tsv(paste0(folder,x)) %>% 
    left_join(df_LUT,by = c("gene"="gene_ids")) %>% 
    arrange(padj)
  
  return(df)
}) %>% 
  setNames(str_remove_all(file_shr,pattern = "res_|.tsv"))
saveRDS(list_shr,"out/object/res_shr.rds")

# explore the dataset -----------------------------------------------------
# what are the biotype in the dataset
df_raw <- list_raw %>% 
  bind_rows(.id = "comparison") %>% 
  mutate(sig = case_when(padj<0.05~"sig",
                         T~"non_sig"))

df_shr <- list_shr %>% 
  bind_rows(.id = "comparison") %>% 
  mutate(sig = case_when(padj<0.05~"sig",
                         T~"non_sig"))

table(df_raw$biotype,df_raw$comparison)
table(df_shr$biotype,df_shr$comparison)

# how many ribosomal are in the dataset
df_raw %>% 
  dplyr::filter(biotype=="rRNA") %>% 
  group_by(gene,biotype) %>% 
  summarise() %>% 
  print(n=40)

# save the barplot of the number of genes in the dataset
df_raw %>% 
  group_by(gene,biotype) %>% 
  summarise() %>% 
  ungroup() %>% 
  group_by(biotype) %>% 
  summarise(n_gene = n()) %>% 
  dplyr::filter(!is.na(biotype)) %>% 
  mutate(biotype = fct_reorder(biotype,-n_gene)) %>% 
  ggplot(aes(x=biotype,y=n_gene))+geom_col()+theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))

# do the same for the number of counts
test <- read_tsv("out/table/featureCounts_Multimapping_fraction_output_ALL_unfiltered_norm.tsv")
test %>% 
  pivot_longer(names_to = "sample",values_to = "norm_counts",-gene) %>% 
  left_join(df_LUT,by = c("gene"="gene_ids")) %>% 
  group_by(biotype) %>% 
  summarise(tot = sum(norm_counts)) %>% 
  arrange(desc(tot)) %>% 
  # dplyr::filter(!is.na(GENEBIOTYPE)) %>% 
  mutate(biotype = fct_reorder(biotype,-tot)) %>% 
  ggplot(aes(x=biotype,y=tot))+geom_col()+theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))

# split the normalized counts per condition
lut <- dds2_filtered@colData %>% data.frame()
test %>% 
  pivot_longer(names_to = "sample",values_to = "norm_counts",-gene) %>% 
  left_join(df_LUT,by = c("gene"="gene_ids")) %>% 
  left_join(lut,by = c("sample")) %>% 
  group_by(condition,biotype) %>% 
  summarise(tot_type = sum(norm_counts)) %>% 
  arrange(desc(tot_type)) %>% 
  ungroup() %>% 
  group_by(condition) %>% 
  mutate(tot = sum(tot_type),
         prop = tot_type/tot) %>% 
  # dplyr::filter(!is.na(GENEBIOTYPE)) %>% 
  mutate(biotype = fct_reorder(biotype,-prop)) %>% 
  ggplot(aes(x=biotype,y=prop,fill=condition))+geom_col(position = "dodge")+theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90),strip.background = element_blank())
facet_wrap(~condition)

# maplot of the datasets
df_raw %>% 
  ggplot() +
  geom_point(data = df_raw %>% dplyr::filter(sig=="non_sig"),aes(x=baseMean,y=log2FoldChange),alpha=0.1,col="gray") +
  geom_point(data = df_raw %>% dplyr::filter(sig=="sig"),aes(x=baseMean,y=log2FoldChange),alpha=0.7,col="red") +
  theme_bw()+facet_wrap(~comparison)+theme(strip.background = element_blank())+scale_x_log10()+geom_hline(yintercept = 0,col="black",linetype="dashed")
ggsave("out/image/MAplot_res_raw.pdf",width = 18,height = 4)

df_shr %>% 
  ggplot() +
  geom_point(data = df_shr %>% dplyr::filter(sig=="non_sig"),aes(x=baseMean,y=log2FoldChange),alpha=0.1,col="gray") +
  geom_point(data = df_shr %>% dplyr::filter(sig=="sig"),aes(x=baseMean,y=log2FoldChange),alpha=0.7,col="red") +
  theme_bw()+facet_wrap(~comparison)+theme(strip.background = element_blank())+scale_x_log10()+geom_hline(yintercept = 0,col="black",linetype="dashed")
ggsave("out/image/MAplot_res_shr.pdf",width = 18,height = 4)

# only one condition
df_shr %>% 
  filter(comparison %in% c("RIGIvsSNRP70_shr")) %>% 
  ggplot() +
  geom_point(data = df_shr %>% dplyr::filter(sig=="non_sig") %>% filter(comparison %in% c("RIGIvsSNRP70_shr")),aes(x=baseMean,y=log2FoldChange),alpha=0.1,col="gray") +
  geom_point(data = df_shr %>% dplyr::filter(sig=="sig") %>% filter(comparison %in% c("RIGIvsSNRP70_shr")),aes(x=baseMean,y=log2FoldChange),alpha=0.7,col="red") +
  theme_bw()+facet_wrap(~comparison)+theme(strip.background = element_blank())+scale_x_log10()+geom_hline(yintercept = 0,col="black",linetype="dashed")
ggsave("out/image/MAplot_res_shr2.pdf",width = 4,height = 4)

# highlight the genes in the plot
df_shr %>% 
  ggplot() +
  geom_point(data = df_shr %>% dplyr::filter(sig=="non_sig"),aes(x=baseMean,y=log2FoldChange),alpha=0.1,col="gray") +
  geom_point(data = df_shr %>% dplyr::filter(sig=="sig"),aes(x=baseMean,y=log2FoldChange),alpha=0.7,col="red") +
  ggrepel::geom_text_repel(data = df_shr %>% dplyr::filter(sig=="sig"),aes(x=baseMean,y=log2FoldChange,label = gene))+
  theme_bw()+facet_wrap(~comparison)+theme(strip.background = element_blank())+scale_x_log10()+geom_hline(yintercept = 0,col="black",linetype="dashed")
ggsave("out/image/MAplot_res_shr_label.pdf",width = 25,height = 8)

#
df_shr %>% 
  filter(comparison %in% c("RIGIvsSNRP70_shr")) %>% 
  ggplot() +
  geom_point(data = df_shr %>% dplyr::filter(sig=="non_sig") %>% filter(comparison %in% c("RIGIvsSNRP70_shr")),aes(x=baseMean,y=log2FoldChange),alpha=0.1,col="gray") +
  geom_point(data = df_shr %>% dplyr::filter(sig=="sig") %>% filter(comparison %in% c("RIGIvsSNRP70_shr")),aes(x=baseMean,y=log2FoldChange),alpha=0.7,col="red") +
  ggrepel::geom_text_repel(data = df_shr %>% dplyr::filter(sig=="sig") %>% filter(comparison %in% c("RIGIvsSNRP70_shr")),aes(x=baseMean,y=log2FoldChange,label = gene))+
  theme_bw()+facet_wrap(~comparison)+theme(strip.background = element_blank())+scale_x_log10()+geom_hline(yintercept = 0,col="black",linetype="dashed")
ggsave("out/image/MAplot_res_shr_label2.pdf",width = 8,height = 8)

# highlight the ribosomal genes
df_shr %>% dplyr::filter(biotype=="rRNA") %>% 
  print(n=120)
df_shr %>% 
  ggplot() +
  geom_point(data = df_shr,aes(x=baseMean,y=log2FoldChange),alpha=0.1,col="gray") +
  # geom_point(data = df_shr %>% dplyr::filter(sig=="non_sig"),aes(x=baseMean,y=log2FoldChange),alpha=0.1,col="gray") +
  # geom_point(data = df_shr %>% dplyr::filter(sig=="sig"),aes(x=baseMean,y=log2FoldChange),alpha=0.7,col="red") +
  ggrepel::geom_text_repel(data = df_shr %>% dplyr::filter(biotype=="rRNA"),aes(x=baseMean,y=log2FoldChange,label = gene),max.overlaps = 20,segment.alpha=0.3)+
  geom_point(data = df_shr %>% dplyr::filter(biotype=="rRNA",sig=="sig"),aes(x=baseMean,y=log2FoldChange),col="red",alpha=0.7) +
  geom_point(data = df_shr %>% dplyr::filter(biotype=="rRNA",sig=="non_sig"),aes(x=baseMean,y=log2FoldChange),col="black",alpha=0.7) +
  theme_bw()+facet_wrap(~comparison)+theme(strip.background = element_blank())+scale_x_log10()+geom_hline(yintercept = 0,col="black",linetype="dashed")
ggsave("out/image/MAplot_res_shr_label_rRNA.pdf",width = 30,height = 10)

#
df_shr %>%
  filter(comparison %in% c("RIGIvsSNRP70_shr")) %>% 
  ggplot() +
  geom_point(data = df_shr %>% filter(comparison %in% c("RIGIvsSNRP70_shr")),aes(x=baseMean,y=log2FoldChange),alpha=0.1,col="gray") +
  # geom_point(data = df_shr %>% dplyr::filter(sig=="non_sig"),aes(x=baseMean,y=log2FoldChange),alpha=0.1,col="gray") +
  # geom_point(data = df_shr %>% dplyr::filter(sig=="sig"),aes(x=baseMean,y=log2FoldChange),alpha=0.7,col="red") +
  ggrepel::geom_text_repel(data = df_shr %>% dplyr::filter(biotype=="rRNA") %>% filter(comparison %in% c("RIGIvsSNRP70_shr")),aes(x=baseMean,y=log2FoldChange,label = gene),max.overlaps = 20,segment.alpha=0.3)+
  geom_point(data = df_shr %>% dplyr::filter(biotype=="rRNA",sig=="sig") %>% filter(comparison %in% c("RIGIvsSNRP70_shr")),aes(x=baseMean,y=log2FoldChange),col="red",alpha=0.7) +
  geom_point(data = df_shr %>% dplyr::filter(biotype=="rRNA",sig=="non_sig") %>% filter(comparison %in% c("RIGIvsSNRP70_shr")),aes(x=baseMean,y=log2FoldChange),col="black",alpha=0.7) +
  theme_bw()+facet_wrap(~comparison)+theme(strip.background = element_blank())+scale_x_log10()+geom_hline(yintercept = 0,col="black",linetype="dashed")
ggsave("out/image/MAplot_res_shr_label_rRNA2.pdf",width = 8,height = 8)

df_raw %>% 
  ggplot() +
  geom_point(data = df_raw,aes(x=baseMean,y=log2FoldChange),alpha=0.1,col="gray") +
  # geom_point(data = df_shr %>% dplyr::filter(sig=="non_sig"),aes(x=baseMean,y=log2FoldChange),alpha=0.1,col="gray") +
  # geom_point(data = df_shr %>% dplyr::filter(sig=="sig"),aes(x=baseMean,y=log2FoldChange),alpha=0.7,col="red") +
  ggrepel::geom_text_repel(data = df_raw %>% dplyr::filter(biotype=="rRNA"),aes(x=baseMean,y=log2FoldChange,label = gene),max.overlaps = 20,segment.alpha=0.3)+
  geom_point(data = df_raw %>% dplyr::filter(biotype=="rRNA",sig=="sig"),aes(x=baseMean,y=log2FoldChange),col="red",alpha=0.7) +
  geom_point(data = df_raw %>% dplyr::filter(biotype=="rRNA",sig=="non_sig"),aes(x=baseMean,y=log2FoldChange),col="black",alpha=0.7) +
  theme_bw()+facet_wrap(~comparison)+theme(strip.background = element_blank())+scale_x_log10()+geom_hline(yintercept = 0,col="black",linetype="dashed")
ggsave("out/image/MAplot_res_raw_label_rRNA.pdf",width = 30,height = 10)

# plot some genes ---------------------------------------------------------
GOI <- c("RNR2","RNR1")
lut <- dds2_filtered@colData %>% data.frame()
counts(dds2_filtered,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  filter(symbol %in% GOI) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  # left_join(MR,by = "sample") %>%
  left_join(lut,by = "sample") %>%
  mutate(count_norm_adj = count + 0.5)%>%
  ggplot(aes(x=condition,y = count_norm_adj,label=sample))+geom_point(position = position_jitter(width = 0.1,seed = 1),alpha=0.6)+
  ggrepel::geom_text_repel(position = position_jitter(seed = 1))+
  facet_wrap(~symbol,scales = "free")+scale_y_log10()+ theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("out/image/scatterplot_RNR2_RNR1.pdf",width = 12,height = 6)


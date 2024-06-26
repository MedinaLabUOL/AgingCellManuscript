# libraries ---------------------------------------------------------------
library("tidyverse")
library("DESeq2")

# sample input ------------------------------------------------------------
# library("pasilla")
# pasCts <- system.file("extdata",
#                       "pasilla_gene_counts.tsv",
#                       package="pasilla", mustWork=TRUE)
# pasAnno <- system.file("extdata",
#                        "pasilla_sample_annotation.csv",
#                        package="pasilla", mustWork=TRUE)
# 
# cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
# 
# cts[1:10,1:3]

# read in the data --------------------------------------------------------
# read in all the files
folder <- "../data/05_featureCounts_out_50/"
file <- dir(folder) %>% 
  str_subset(pattern = "featureCounts_Multimapping_fraction_output.txt$")

sample_id <- file %>% 
  str_extract_all(pattern = "PN0401_\\d+") %>% 
  unlist()

# # sample read one file
# Cts <-"../data/05_featureCounts_out_50/PN0401_0001_featureCounts_Multimapping_fraction_output.txt"
# test <- read.delim(Cts,skip = 1)
# colnames(test) <- c("Geneid","Chr","Start","End","Strand","Length","Count")
# # in case of fractional count I should round them to the closest integer https://support.bioconductor.org/p/75848/
# test_fix <- test %>%
#   mutate(Count_rounded = round(Count))
# 
# head(test_fix)
# 
# df_out <- test_fix %>%
#   dplyr::select(Geneid,Count_rounded)
# 
# 
# 
# # are there any replicates in the geneid ?
# df_out %>%
#   group_by(Geneid) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))
# 
# test_fix[1:10,]

# read in all the files
list_reads <- lapply(sample_id,function(x){
  Cts <-paste0("../data/05_featureCounts_out_50/",x,"_featureCounts_Multimapping_fraction_output.txt")
  test <- read.delim(Cts,skip = 1)
  colnames(test) <- c("Geneid","Chr","Start","End","Strand","Length","Count")
  # in case of fractional count I should round them to the closest integer https://support.bioconductor.org/p/75848/
  test_fix <- test %>% 
    mutate(Count_rounded = round(Count))
  
  df_out <- test_fix %>% 
    dplyr::select(Geneid,Count_rounded)
  
  return(df_out)
}) %>% 
  setNames(sample_id)

# confirm the dimension of all the files is the same
lapply(list_reads,function(x){dim(x)})

# confirm no dataset have replicated gene names
lapply(list_reads,function(x){x %>% group_by(Geneid) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n))
  })

# make the list as a single matrix
df_exp <- pmap(list(list_reads,names(list_reads)),function(x,name){
  x %>%
    dplyr::rename(!!name := "Count_rounded")
}) %>% 
  purrr::reduce(left_join,by=c("Geneid"))

# save the table
df_exp %>% 
  write_tsv("out/table/featureCounts_Multimapping_fraction_output_ALL_raw.tsv")

mat_exp <- df_exp %>% 
  column_to_rownames("Geneid") %>% 
  as.matrix()

# tot counts per sample
colSums(mat_exp)

# build the annotation ----------------------------------------------------
coldata <- data.frame(sample = colnames(mat_exp)) %>% 
  mutate(condition = case_when(sample %in% c("PN0401_0001","PN0401_0002","PN0401_0003","PN0401_0004")~"RIP_RIGI",
                               sample %in% c("PN0401_0005")~"RIP_SNRP70",
                               sample %in% c("PN0401_0006")~"NO_RIP")) %>% 
  mutate(rowname = sample) %>% 
  column_to_rownames("rowname")

# We examine the count matrix and column data to see if they are consistent in terms of sample order.
head(mat_exp,2)

coldata

# Note that these are not in the same order with respect to samples!
# It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order.

# As they are not in the correct order as given, we need to re-arrange one or the other so that they are consistent in terms of sample order (if we do not, later functions would produce an error). We additionally need to chop off the "fb" of the row names of coldata, so the naming is consistent.

rownames(coldata) <- rownames(coldata)
all(rownames(coldata) %in% colnames(mat_exp))

all(rownames(coldata) == colnames(mat_exp))

mat_exp <- mat_exp[, rownames(coldata)]
all(rownames(coldata) == colnames(mat_exp))

# build the design --------------------------------------------------------
design <- model.matrix(~ coldata$condition)
colnames(design) <- c("intercept","RIP_RIGI","RIP_SNRP70")

# save the disegn
saveRDS(design,"out/object/design.rds")

# build the deseq2 objects ------------------------------------------------
# library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = mat_exp,
                              colData = coldata,
                              design = design)
dds

# remove potential non infirmative genes
dds_filter <- dds[rowSums(counts(dds)) > 10, ]

nrow(dds)
nrow(dds_filter)

# generate the filtered and unfiltered processed data
dds2 <- DESeq(dds)
dds2_filter <- DESeq(dds_filter)

# save the outputs --------------------------------------------------------
saveRDS(dds2,"out/object/dds2_unfiltered.rds")
saveRDS(dds2_filter,"out/object/dds2_filtered.rds")

# scale the dataset
vds_filter <- vst(dds2_filter, blind = F)
saveRDS(vds_filter,file = "out/object/vds_filter.rds")

counts(dds2,normalized=T) %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>%
  write_tsv("out/table/featureCounts_Multimapping_fraction_output_ALL_unfiltered_norm.tsv")

counts(dds2_filter,normalized=T) %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>%
  write_tsv("out/table/featureCounts_Multimapping_fraction_output_ALL_filtered_norm.tsv")
  

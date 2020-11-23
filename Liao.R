## analysis Liao
library(tidyverse)

## Cells info
meta <- read_delim("data/Liao/meta.tsv", delim = "\t")

meta <- filter(meta, celltype == "Macrophages")

## Preprocess Liao data
## single cell data downloaded From 

## Merge cells to recreate bulk
# library(Seurat)
# s <- readRDS("nCoV.rds")
# 
# exp <- GetAssayData(s, slot = "counts")
# 
# merge_counts <- function(ids){
#   data.frame(gene = rownames(exp), counts = Matrix::rowSums(exp[, ids]))
# }
# 
# merged <- meta %>% 
#   group_by(sample_new) %>%
#   summarise(merge_counts(ID))
# 
# merged <- pivot_wider(merged, names_from = sample_new, values_from = counts)
# 
# write_delim(merged, path = "data/Liao/Liao_merged.tsv", delim = "\t")



## Differential analysis
library(limma)
library(edgeR)

trs <- read_delim("data/Liao/Liao_merged.tsv", delim = "\t")

dge <- DGEList(counts = as.matrix(column_to_rownames(trs, "gene")))

coldata <- colnames(trs)[-1] %>% 
  enframe() %>% 
  mutate(group = str_sub(value, start = 1, end = -2)) %>% 
  select(-name) %>% 
  as.data.frame() %>% 
  column_to_rownames("value")

coldata <- coldata[colnames(dge), , drop = FALSE]

design <- model.matrix(~ -1 + group, data = coldata)
colnames(design) <- make.names(colnames(design))

keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)
# convert to lo2cpm
v1 <- voom(dge, design = design, plot = TRUE)
res_fit <- lmFit(v1)

write_delim(rownames_to_column(as.data.frame(v1$E), "gene"), delim = "\t", path = "data/Liao/Liao_norm.tsv")


## Comparison Severe vs Healthy
cont_S_HC <- makeContrasts(groupS-groupHC, levels=design)
res_cont_S_HC <- contrasts.fit(res_fit, cont_S_HC)
res_cont_S_HC <- eBayes(res_cont_S_HC)

res_export_S_HC <- topTable(res_cont_S_HC, p.value = 1, number=Inf)
res_export_S_HC <- res_export_S_HC %>% 
  as.data.frame() %>% 
  rownames_to_column("genes")
res_export <- res_export[order(res_export$P.Value), c("genes", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]

## Comparison Moderate vs Healthy
cont_HC_M <- makeContrasts(groupM - groupHC, levels=design)
res_cont_HC_M <- contrasts.fit(res_fit, cont_HC_M)
res_cont_HC_M <- eBayes(res_cont_HC_M)

res_export_HC_M <- topTable(res_cont_HC_M, p.value = 1, number=Inf)
res_export_HC_M <- res_export_HC_M %>% 
  as.data.frame() %>% 
  rownames_to_column("genes")
res_export_HC_M <- res_export_HC_M[order(res_export_HC_M$P.Value), c("genes", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]


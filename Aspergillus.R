library(limma)
library(edgeR)
library(tidyverse)
library(plyr)



genes_anno <- read.csv("data/Aspergillus/tableannot.csv",
                       row.names = 1,
                       stringsAsFactors = FALSE)

samples_desc <- read_csv("data/Aspergillus/SraRunTable.txt")
samples_desc <- samples_desc %>% 
  select(barcode = `GEO_Accession (exp)`, infection_status) %>% 
  unique() %>% 
  as.data.frame() %>% 
  column_to_rownames("barcode")



counts <- read_csv("data/Aspergillus/tablecounts_raw.csv")
counts <- column_to_rownames(as.data.frame(counts), "X1")

#focus on protein coding genes
allcor <- read.table("data/hgnc_complete_set.txt",
                     header = TRUE,
                     sep = "\t",
                     fill = T,
                     comment.char = "~",
                     quote = "")

allcor_id <- as.character(allcor[which(allcor$locus_group == "protein-coding gene"), "ensembl_gene_id"])

allcor1 <- allcor[which(allcor$name %in% grep("tumor", grep("ribosom|mitochon", allcor$name, value = TRUE), value = TRUE, invert = TRUE) & allcor$locus_Group == "protein-coding gene"),]

allcor_id <- allcor_id[which(!allcor_id %in% allcor1$ensembl_gene_id)]

allcor2 <- allcor[which(as.character(allcor$ensembl_gene_id) %in% allcor_id), ]
allcor2$chr <- sapply(strsplit(gsub("p", "q", allcor2$location), "q"), "[", 1)
allcor2 <- allcor2[which(nchar(allcor2$chr) < 3), ]
allcor2 <- allcor2[which(!allcor2$chr %in% c("un", "Y")), ]

allcor_id <- allcor_id[which(allcor_id %in% as.character(allcor2$ensembl_gene_id))]

counts <- counts[intersect(rownames(counts), allcor_id), ]
counts <- counts[, rownames(samples_desc)]


# Normalizing raw counts
dge <- DGEList(counts = counts, genes = rownames(counts), samples = samples_desc)

design <- model.matrix( ~ -1 + infection_status, data = dge$samples)
colnames(design) <- make.names(colnames(design))

keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

# convert to lo2cpm
v1 <- voomWithQualityWeights(dge, design = design, plot = TRUE)

v1t <- v1$E

write.table(v1t, "results/norm_counts_Aspergillus.txt", row.names = TRUE, col.names = TRUE, sep = "\t",
            quote = FALSE)

res_fit <- lmFit(v1, design = design)

cont <- makeContrasts(inf_2_wt = infection_status2h.post.infection-infection_statusuninfected, 
                      inf_6_wt = infection_status6h.post.infection-infection_statusuninfected, 
                      inf_6_inf_2 = infection_status6h.post.infection-infection_status2h.post.infection,
                      levels = design)

res_cont <- contrasts.fit(res_fit, cont)
res_cont <- eBayes(res_cont)


prep_res <- function(tab){
  tab$genes <- mapvalues(rownames(tab),
                         genes_anno$gene_id,
                         genes_anno$gene_name,
                         warn_missing = FALSE)
  tab <- tab[order(tab$P.Value), c("genes", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
}

table_cont_inf_2_wt <- prep_res(topTable(res_cont, p.value = 1, number = Inf, coef = "inf_2_wt"))
table_cont_inf_6_wt <- prep_res(topTable(res_cont, p.value = 1, number = Inf, coef = "inf_6_wt"))
table_cont_inf_6_inf_2 <- prep_res(topTable(res_cont, p.value = 1, number = Inf, coef = "inf_6_inf_2"))


write_csv(table_cont_inf_2_wt, path = "results/genes_Aspergillus_2h.csv")
write_csv(table_cont_inf_6_wt, path = "results/genes_Aspergillus_6h.csv")





library(limma)
library(edgeR)
library(tidyverse)
library(plyr)



genes_anno <- read.csv("data/Leishmania/tableannot.csv",
                       row.names = 1,
                       stringsAsFactors = FALSE)

samples_desc <- read_csv("data/Leishmania/SraRunTable.txt")
samples_desc <- filter(samples_desc, str_detect(treatment, "Infected with Leishmania major") | str_detect(treatment, "Uninfected"))

samples_desc <- separate(samples_desc, treatment, into = c("treat", "time"), remove = FALSE, sep = ",")

samples_desc <- mutate(samples_desc, treat = if_else(str_detect(treat, "Uninfected"), "Uninfected", "Infected"), time = str_remove_all(time, " "))
samples_desc <- unite(samples_desc, col = treat, treat, time)

samples_desc <- samples_desc %>% 
  select(barcode = Run, treatment, treat) %>% 
  unique() %>% 
  as.data.frame() %>% 
  column_to_rownames("barcode")



counts <- read_csv("data/Leishmania/tablecounts_raw.csv")
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

samples_desc <- samples_desc[intersect(colnames(counts), rownames(samples_desc)),,drop = FALSE]

counts <- counts[intersect(rownames(counts), allcor_id), ]
counts <- counts[, rownames(samples_desc)]


# Normalizing raw counts
dge <- DGEList(counts = counts, genes = rownames(counts), samples = samples_desc)

design <- model.matrix( ~ -1 + treat, data = dge$samples)
colnames(design) <- make.names(colnames(design))

keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

# convert to lo2cpm
v1 <- voom(dge, design = design, plot = TRUE)

v1t <- v1$E

write.table(v1t, "results/norm_counts_Leishmania.txt", row.names = TRUE, col.names = TRUE, sep = "\t",
            quote = FALSE)

res_fit <- lmFit(v1, design = design)

cont <- makeContrasts(inf_4h = treatInfected_4hpi - treatUninfected_4hpi,
                      inf_24h = treatInfected_24hpi - treatUninfected_24hpi, 
                      inf_48h = treatInfected_48hpi - treatUninfected_48hpi, 
                      inf_72h = treatInfected_72hpi - treatUninfected_72hpi,
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

table_cont_inf_4h <- prep_res(topTable(res_cont, p.value = 1, number = Inf, coef = "inf_4h"))
table_cont_inf_24h <- prep_res(topTable(res_cont, p.value = 1, number = Inf, coef = "inf_24h"))
table_cont_inf_48h <- prep_res(topTable(res_cont, p.value = 1, number = Inf, coef = "inf_48h"))
table_cont_inf_72h <- prep_res(topTable(res_cont, p.value = 1, number = Inf, coef = "inf_72h"))



write_csv(table_cont_inf_4h, path = "results/genes_Leishmania_4h.csv")
write_csv(table_cont_inf_24h, path = "results/genes_Leishmania_24h.csv")




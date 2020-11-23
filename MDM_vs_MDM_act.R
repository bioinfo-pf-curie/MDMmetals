## Analysis MDM vs activated MDM

library(tidyverse)
library(plyr)
library(readxl)
library(edgeR)
library(limma)

## ----samples--------------------------------------------------------------
desc <- read_excel("data/samples_desc.xlsx")

## ----norm---------------------------------------------------
genes_anno <- read.csv("data/tableannot.csv",
                       row.names = 1,
                       stringsAsFactors = FALSE)

counts <- read.csv("data/tablecounts_raw.csv", row.names = 1)

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

coldata <- data.frame("Condition" = desc$Condition,
                      Group = desc$Group,
                      "Blood" = paste("Blood", desc$Donor, sep = ""))
rownames(coldata) <- desc$SampleID

counts <- counts[, rownames(coldata)]

# Normalizing raw counts
coldata$Group <- factor(coldata$Group, levels = c("WT", "Activated"))
# coldata <- coldata[colnames(dge),]
dge <- DGEList(counts = counts, genes = rownames(counts), samples = coldata)

design <- model.matrix( ~ -1 + Group, data = dge$samples)
colnames(design) <- make.names(colnames(design))

keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)


# estimate the correlation between the patients
cor_fit <- duplicateCorrelation(voom(dge), block = dge$targets$Blood)

v2 <- voom(dge,
           design = design,
           plot = TRUE,
           block = dge$samples$Blood,
           correlation = cor_fit$consensus.correlation
)

cor_fit <- duplicateCorrelation(v2, block = coldata$Blood)

res_fit <- lmFit(v2,
                 design = design,
                 block = coldata$Blood,
                 correlation = cor_fit$consensus.correlation
)


padj_th <- 0.05

#"GM-CSF+LPS+IFN","GM-CSF
contrast.matrix <- makeContrasts(GroupActivated-GroupWT, levels = design)
res_cont_act_wt <- contrasts.fit(res_fit, contrast.matrix)
res_cont_act_wt <- eBayes(res_cont_act_wt)

table_cont_act_wt <- topTable(res_cont_act_wt, p.value = 1, number = Inf)
table_cont_act_wt$genes <- mapvalues(rownames(table_cont_act_wt),
            genes_anno$gene_id,
            genes_anno$gene_name,
            warn_missing = FALSE)
table_cont_act_wt <- table_cont_act_wt[order(table_cont_act_wt$P.Value), c("genes", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
write_csv(table_cont_act_wt, path = "results/genes_GMCSF_activated_vs_WT.csv")



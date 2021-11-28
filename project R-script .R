##Library 
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(apeglm)
library (dplyr)
library(biomaRt)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(GO.db)
library(GOstats)
library(ggplot2)
library("pheatmap")

#read meta data [file used: CLP = 85,87.htseq.out | CD8 = 17,20.htseq.out]
metadata <- read_csv("metadata.csv")
file <- paste0(metadata$Sample, ".htseq.out")
file.exists(file)

#create new datatable
sample_df <- data.frame(sampleName = metadata$Sample,
                        fileName = file,
                        cell_type = metadata$`Cell type`)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sample_df,
                                       design = ~ cell_type)

#set control condition: ref = 'CLP'
ddsHTSeq$cell_type <- relevel(ddsHTSeq$cell_type, ref = "CLP")


#pre-filtering of lowly expressed genes
keep <- rowSums(counts(ddsHTSeq)) >= 10
deseq <- ddsHTSeq[keep, ]

ddsHTSeq$cell_type
deseq

#DESeq modelling 
dds <- DESeq(deseq)

##sample clustering with log-transformed data
#PCA plot
rld <- rlog(dds)
plotPCA(rld,intgroup = "cell_type" )

#heatmap
sample_dists <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dists)
colnames(sample_dist_matrix) <- NULL

colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists, 
         col = colours)

#differentially expressed genes

resultsNames(dds)

res1 <- results(dds, name = "cell_type_CD8_vs_CLP")

res2 <- as.data.frame(results(dds, name = "cell_type_CD8_vs_CLP")) %>%
  rownames_to_column("ENSEMBL")

significant_res <- subset(res1, padj < 0.05)

up_in_CD8 <- subset(res1, padj < 0.05 & log2FoldChange > 0 )
down_in_CD8 <- subset(res1, padj < 0.05 & log2FoldChange < 0 )


##result visulization 

#logfold change - MA plot
resLFC <- lfcShrink(dds, coef = "cell_type_CD8_vs_CLP", type = "apeglm" )
plotMA(resLFC, ylim=c(-10,10))

#logfold change - volcanol plot
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
threshold <- res2$padj < padj.cutoff & abs(res2$log2FoldChange) > lfc.cutoff

res2 %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj),colour = threshold)) + geom_point()


#plot the gene with the lowest p.adj value between cell types
plotCounts(dds, gene=which.min(res1$padj), intgroup="cell_type")

#plot specific gene with gene ID
plotCounts(dds, "ENSMUSG00000005474",intgroup = c("cell_type"))

#gene expression heatmap 
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized =T )),
                decreasing = T)

df <- as.data.frame(colData(dds)[,c("cell_type")])
pheatmap(assay(ntd)[select,], cluster_rows = F, show_rownames = F, cluster_cols = F, annotation_colors = df )


#adding annotation to genes 

res2 <- as.data.frame(results(dds, name = "cell_type_CD8_vs_CLP")) %>%
  rownames_to_column("ENSEMBL")

anno <- AnnotationDbi::select(org.Mm.eg.db,keys=res2$ENSEMBL,
                              columns=c("ENSEMBL","SYMBOL","GENENAME","ENTREZID"),
                              keytype="ENSEMBL") %>% filter(!duplicated(ENSEMBL))

dim(res2)
dim(anno)

res <- left_join(res2, anno,by="ENSEMBL")
head(res)

#filtering annotated differentially expressed gene for export
significant_res <- subset(res, padj < 0.05) %>% arrange (, padj)
head(significant_res)

more_upregulated_in_CD8 <- subset(res, padj < 0.05 & log2FoldChange > 0 ) %>% arrange (, padj)
more_downregulated_in_CD8 <- subset(res, padj < 0.05 & log2FoldChange < 0 ) %>% arrange (, padj)

head(more_upregulated_in_CD8)
head(more_downregulated_in_CD8)

#data export
write.csv(as.data.frame(significant_res), file = "sorted cell_type_CD8_vs_CLP.csv")

write.csv(as.data.frame(more_upregulated_in_CD8), file = "sorted upregulated in CD8 vs CLP.csv")

write.csv(as.data.frame(more_downregulated_in_CD8), file = "sorted downregulated in CD8 vs CLP.csv")



##GO Term enrichment analysis
significant_res1 <- subset(res1, padj < 0.05)

annotated_significant_results <- significant_res1

annotated_significant_results$symbol <- mapIds(
  org.Mm.eg.db,
  keys = rownames(annotated_significant_results),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)


annotated_significant_results$entrez <- mapIds(
  org.Mm.eg.db,
  keys = rownames(annotated_significant_results),
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)


all_genes <- annotated_significant_results %>% 
  as.data.frame() %>% 
  pull(entrez) %>% 
  unique()


#top10 GO terms for upregulated genes in CD8
genes_upregulated <- annotated_significant_results %>% 
  as.data.frame() %>% 
  filter(log2FoldChange > 4) %>% 
  pull(entrez) %>% 
  unique()

go_bp_upregulated <- hyperGTest(new("GOHyperGParams",
                                    geneIds = genes_upregulated,
                                    universeGeneIds = all_genes,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.01,
                                    conditional = FALSE,
                                    testDirection = "over"))

top10_enriched_GO_terms_up <- go_bp_upregulated %>% summary() %>% head(10)

top10_enriched_GO_terms_up

#top10 GO terms fordownregulated genes in CD8

genes_downregulated <- annotated_significant_results %>% 
  as.data.frame() %>% 
  filter(log2FoldChange < 4) %>% 
  pull(entrez) %>% 
  unique()

go_bp_downregulated <- hyperGTest(new("GOHyperGParams",
                                    geneIds = genes_downregulated,
                                    universeGeneIds = all_genes,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.01,
                                    conditional = FALSE,
                                    testDirection = "over"))

top10_enriched_GO_terms_down <- go_bp_downregulated %>% summary() %>% head(10)

top10_enriched_GO_terms_down


#For annotated peaks --------------------------------------------------
#read peak annotated files
cd8_ann <- read_csv("cd8_annotated.csv")
clp_cd8_ann <- read_csv("clp_cd8_annotated.csv")
clp_ann <- read_csv("clp_nocd8_annotated.csv")

#rename Entrez ID columns to match in tables
names(more_upregulated_in_CD8)[names(more_upregulated_in_CD8) == "ENTREZID"] <- "Entrez ID"
names(more_downregulated_in_CD8)[names(more_downregulated_in_CD8) == "ENTREZID"] <- "Entrez ID"

#return joined rows with matching entrez ids
cd8_upreg = merge(x=cd8_ann, y=more_upregulated_in_CD8, by="Entrez ID") %>% 
  arrange (, padj)
cd8_upreg = subset(cd8_upreg, select = (c("Chr", "Start", "End", "Peak Score")))
write.csv(cd8_upreg, file= "~/Desktop/MICB405/MICB405Project/cd8_upreg.csv", row.names=FALSE)

cd8_downreg = merge(x=cd8_ann, y=more_downregulated_in_CD8, by="Entrez ID") %>% 
  arrange (, padj) 
cd8_downreg = subset(cd8_upreg, select = (c("Chr", "Start", "End", "Peak Score")))
write.csv(cd8_downreg, file= "~/Desktop/MICB405/MICB405Project/cd8_upreg.csv", row.names=FALSE)

clp_downreg = merge(x=clp_ann, y=more_upregulated_in_CD8, by="Entrez ID") %>%
  arrange (, padj)
clp_downreg = subset(clp_downreg, select = (c("Chr", "Start", "End", "Peak Score")))
write.csv(clp_downreg, file= "~/Desktop/MICB405/MICB405Project/clp_downreg.csv", row.names=FALSE)

clp_upreg = merge(x=clp_ann, y=more_downregulated_in_CD8, by="Entrez ID") %>%
  arrange (, padj)
clp_upreg = subset(clp_upreg, select = (c("Chr", "Start", "End", "Peak Score")))
write.csv(clp_upreg, file= "~/Desktop/MICB405/MICB405Project/clp_upreg.csv", row.names=FALSE)





 







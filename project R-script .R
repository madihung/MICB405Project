##Library 
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(apeglm)

#read meta data

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

ddsHTSeq
deseq

#DESeq modelling 
dds <- DESeq(deseq)

#sample clustering with log-transformed data
rld <- rlog(dds)
plotPCA(rld,intgroup = "cell_type" )

sample_dists <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dists)
colnames(sample_dist_matrix) <- NULL

colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists, 
         col = colours)



#filtering differentially expressed gene
resultsNames(dds)

res <- results(dds, name = "cell_type_CD8_vs_CLP")
head(res)
summary(res)

significant_res <- subset(res, padj < 0.05)
head(significant_res)

more_upregulated_in_CD8 <- subset(res, padj < 0.05 & log2FoldChange > 0 )
more_downregulated_in_CD8 <- subset(res, padj < 0.05 & log2FoldChange < 0 )


##result visulization 

#logfold change
resLFC <- lfcShrink(dds, coef = "cell_type_CD8_vs_CLP", type = "apeglm" )
resLFC
plotMA(resLFC, ylim=c(-10,10))

#plot the gene with the lowest p.adj value between cell types
plotCounts(dds, gene=which.min(res$padj), intgroup="cell_type")


#data export
write.csv(as.data.frame(significant_res), file = "cell_type_CD8_vs_CLP.csv")

write.csv(as.data.frame(more_upregulated_in_CD8), file = "upregulated in CD8 vs CLP.csv")

write.csv(as.data.frame(more_downregulated_in_CD8), file = "downregulated in CD8 vs CLP.csv")

##GO Term enrichment analysis
library(AnnotationDbi)
library(org.Mm.eg.db)
library(GO.db)
library(GOstats)

annotated_significant_results <- significant_res

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

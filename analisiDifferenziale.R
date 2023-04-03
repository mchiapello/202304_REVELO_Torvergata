library(tximeta)
library(DESeq2)

coldata <- read.csv("sample_table.csv", 
                    row.names=1, stringsAsFactors=FALSE)

coldata$names <- coldata$Run

files <- fs::dir_ls("quants_Gencode/", 
                    recurse = TRUE, regexp = "quant.sf")

coldata$files <- files

rownames(coldata) <- coldata$Run

se <- tximeta(coldata)

rowRanges(se)

################################################################################
library(airway)
library(DESeq2)

data(gse)

gse

metadata(gse)


str(gse, max.level = 2)


rowData(gse)
colData(gse)
gse

assay(gse)
assays(gse)
assays(gse)$abundance

## Pre-filtering the dataset
nrow(gse)
keep <- rowSums(assay(gse)) > 1
dds <- gse[keep,]
nrow(dds)

data <- data.frame(assay(dds))

library(ggplot2)

data |> 
  ggplot(aes(x = SRR1039508)) +
  geom_histogram(bins = 200)

data |> 
  tidyr::pivot_longer(cols = everything(),
                      names_to = "Sample",
                      values_to = "Count") |> 
  ggplot(aes(x = Count)) +
  geom_histogram(bins = 200) +
  facet_wrap(vars(Sample))

dds <- DESeqDataSet(dds, design = ~ donor + condition)

vs <- vst(dds, blind = FALSE)

head(assay(dds))
head(assay(vs))


sampleDist <- dist(t(assay(vs)))
sampleDistMatrix <- as.matrix(sampleDist)

rownames(sampleDistMatrix) <- paste(vs$donor, vs$condition, 
                                    sep = " - ")

library(pheatmap)
library(RColorBrewer)

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDist,
         # col = colors,
         cutree_rows = 2,
         cutree_cols = 2)

plotPCA(vs, intgroup = c("donor", "condition"))
plotPCA(vs, intgroup = c("donor"))
plotPCA(vs, intgroup = c("condition"))

pcaData <- plotPCA(vs, intgroup = c("donor", "condition"),
                   returnData = TRUE)

str(pcaData)

percentVar <- round(attr(pcaData, "percentVar") * 100)

pcaData |> 
  ggplot(aes(x = PC1, y = PC2, 
             color = condition,
             shape = donor)) +
  geom_point(size = 5) +
  labs(
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%"),
    title = "PCA with VST data"
  ) +
  theme(plot.title = element_text(hjust = .5))

dds <- DESeq(dds)

res <- results(dds)

summary(res)

table(res$padj < 0.1)

table(res$padj < 0.1 & res$log2FoldChange > 0)


resDF <- res |> 
  data.frame()
  
resDF |> 
  tibble::as_tibble() |> 
  dplyr::mutate(reg = dplyr::case_when(padj < 0.1 & log2FoldChange > 0 ~ "up",
                                padj < 0.1 & log2FoldChange < 0 ~ "down",
                                TRUE ~ "NR")) |> 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             color = reg)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.1))

res.a1 <- results(dds, alpha = 0.0001)

summary(res.a1)

resFC1 <- results(dds, lfcThreshold = 1)
summary(resFC1)



plotMA(res, ylim = c(-5, 5))

library("apeglm")

resultsNames(dds)
resS <- lfcShrink(dds, coef = "condition_Dexamethasone_vs_Untreated",
                  type = "apeglm", res = res)

plotMA(resS, ylim = c(-5, 5))


write.csv(resS, "pippo.csv")






































































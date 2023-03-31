################################################################################
# Script for RNASeq analysis
################################################################################

# Load libraries
library(tximeta)
library(DESeq2)

################################################################################
# Only the teacher
################################################################################
# Import from Salmon quantification
## Read sample table
coldata <- read.csv("sample_table.csv", row.names=1, stringsAsFactors=FALSE)
## Modify sample table to fit the data
coldata$names <- coldata$Run
files <- fs::dir_ls("quants_Gencode/", recurse = TRUE, regexp = "quant.sf")
coldata$files <- files
rownames(coldata) <- coldata$Run
## Check the presence of the files
file.exists(coldata$files)
## Read the data in (transcript-level)
setr <- tximeta(coldata)
# summarize the transcript-level quantifications to the gene level 
sege <- summarizeToGene(setr)
## Explore the RangedSummarizedExperiment
rowRanges(setr)
rowData(setr)
assay(se)
colData(se)
metadata(se)

################################################################################
# Exploratory analysis and visualization
################################################################################
# Load libraries
library(ggplot2)
library(dplyr)
library(airway)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(DESeq2)
# Load data
data(airway)
data(gse)

# RangedSummarizedExperiment modifications
gse$cell <- gse$donor
gse$dex <- gse$condition
gse <- DESeqDataSet(gse, design = ~ cell + dex)

## Pre-filtering the dataset
nrow(gse)
keep <- rowSums(assay(gse)) > 1
dds <- gse[keep,]
nrow(dds)

# have a look at the distribution of raw read counts
data <- data.frame(assay(dds))
ggplot(data) +
  geom_histogram(aes(x = SRR1039508), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

data |> 
  tidyr::pivot_longer(cols = everything(),
                      names_to = "Sample",
                      values_to = "Count") |> 
  ggplot() +
    geom_histogram(aes(x = Count), stat = "bin", bins = 200) +
    xlab("Raw expression counts") +
    ylab("Number of genes") +
  facet_wrap(~Sample)

################################################################################
# The variance stabilizing transformation and the rlog

# variance stabilizing transformation (VST) for negative binomial data
vsd <- vst(dds, blind = FALSE)
head(assay(dds), 3)
head(assay(vsd), 3)

# regularized-logarithm transformation
rld <- rlog(dds, blind = FALSE)
head(assay(dds), 3)
head(assay(rld), 3)

#' we specified blind = FALSE, which means that differences between cell lines 
#' and treatment (the variables in the design) will not contribute to the 
#' expected variance-mean trend of the experiment

df <- bind_rows(
  as_tibble(assay(dds)[, 1:2]) |> 
    mutate(transformation = "original"),
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("original", "log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

df |> 
  filter(transformation != "original") |> 
ggplot(aes(x = x, y = y)) + 
  geom_hex(bins = 80) +
  # coord_fixed() + 
  facet_grid( . ~ transformation)  

################################################################################
# Sample distances
## Compute the distance
sampleDists <- dist(t(assay(vsd)))
sampleDists
## Matrix transformation
sampleDistMatrix <- as.matrix( sampleDists )
## Change ros and column names
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
# colnames(sampleDistMatrix) <- NULL
## Create color palette
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
## Plot the heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

## Plot PCA
plotPCA(vsd, intgroup = c("cell"))
plotPCA(vsd, intgroup = c("dex"))
plotPCA(vsd, intgroup = c( "dex", "cell"))

## Optional
pcaData <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData |> 
  ggplot(aes(x = PC1, y = PC2, color = dex, shape = cell)) +
    geom_point(size =3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    ggtitle("PCA with VST data")

################################################################################
# Differential expression analysis
#' we run the differential expression pipeline on the raw counts with a single call to the function DESeq
dds <- DESeq(dds)

# Result
res <- results(dds)
res
# Get information on each column in results
mcols(res, use.names = TRUE)
# Summarise results
summary(res)

table(res$padj < 0.1)
table(res$padj < 0.1 & res$log2FoldChange > 0)
table(res$padj < 0.1 & res$log2FoldChange < 0)
2373+1949

## Volcano plot
resDF <- res %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>% 
  tibble::as_tibble() |> 
  mutate(reg = case_when(padj < 0.1 & log2FoldChange > 0 ~ "up",
                         padj < 0.1 & log2FoldChange < 0 ~ "down",
                         TRUE ~ "NR"))

resDF |> count(reg)

ggplot(resDF) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = reg)) +
  ggtitle("Volcano Plot") +
  geom_hline(yintercept = -log10(0.1)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


#' There are two ways to be more strict about which set of genes are considered significant:
#' - lower the false discovery rate threshold (the threshold on padj in the results table)
#' - raise the log2 fold change threshold from 0 using the lfcThreshold argument of results

# 1. lower the false discovery rate
res.a1 <- results(dds, alpha = 0.0001)
table(res.a1$padj < 0.0001)

resDF.a1 <- res.a1 %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>% 
  tibble::as_tibble() |> 
  mutate(reg = case_when(padj < 0.0001 & log2FoldChange > 0 ~ "up",
                         padj < 0.0001 & log2FoldChange < 0 ~ "down",
                         TRUE ~ "NR"))

ggplot(resDF.a1) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = reg)) +
  ggtitle("Volcano Plot") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept = -log10(0.0001)) +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

# 2. raise the log2 fold change

resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

resDF.lfc1 <- resLFC1 %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>% 
  tibble::as_tibble() |> 
  mutate(reg = case_when(padj < 0.1 & log2FoldChange > 1 ~ "up",
                         padj < 0.1 & log2FoldChange < -1 ~ "down",
                         TRUE ~ "NR"))

resDF.lfc1 |> count(reg)

ggplot(resDF.lfc1) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = reg)) +
  ggtitle("Volcano Plot") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept = -log10(0.1)) +
  geom_vline(xintercept = c(-1, 1)) +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

################################################################################
# Plottign results

## Counts plot
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("dex"))

geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"),
                        returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  geom_point()

# Count down 
topGenesD <- resDF |> 
  filter(reg == "down") |> 
  arrange(padj) |> 
  slice(1) |> 
  pull(gene)

geneCounts <- plotCounts(dds, gene = topGenesD, intgroup = c("dex","cell"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  geom_point()

# Count top 
topGenesD <- resDF |> 
  filter(reg == "up") |> 
  arrange(padj) |> 
  slice(1) |> 
  pull(gene)

geneCounts <- plotCounts(dds, gene = topGenesD, intgroup = c("dex","cell"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  geom_point()

# More then one
topGene <- resDF |> 
  # filter(reg == "down") |> 
  arrange(padj) |> 
  slice(1:5) |> 
  pull(gene)

as.data.frame(assay(vsd)[topGene,]) |> 
  tibble::rownames_to_column("gene") |> 
  tidyr::pivot_longer(cols = starts_with("SRR"), 
                      names_to = "names", 
                      values_to = "count") |> 
  dplyr::left_join(as.data.frame(colData(vsd))) |> 
  ggplot(aes(dex, count, color = cell)) + 
  geom_point() + 
  facet_wrap(~gene)

################################################################################
# MA-plot
#' An MA-plot (Dudoit et al. 2002) provides a useful overview for the distribution 
#' of the estimated coefficients in the mode
#' 
#' On the y-axis, the “M” stands for “minus” – subtraction of log values is 
#' equivalent to the log of the ratio 
#' On the x-axis, the “A” stands for “average”
plotMA(res, ylim = c(-5, 5))
plotMA(resLFC1, ylim = c(-5, 5))

#' We need to use the lfcShrink function to shrink the log2 fold changes 
#' for the comparison of dex treated vs untreated samples
#' Shrinking the log2 fold changes will not change the total number of genes that are identified as significantly differentially expressed. The shrinkage of fold change is to help with downstream assessment of results. For example, if you wanted to subset your significant genes based on fold change for further evaluation, you may want to use shruken values. Additionally, for functional analysis tools such as GSEA which require fold change values as input you would want to provide shrunken values.
library("apeglm")
resultsNames(dds)
resS <- lfcShrink(dds, coef="dex_Dexamethasone_vs_Untreated", type="apeglm",
                  res = res)
plotMA(resS, ylim = c(-5, 5))

resL <- lfcShrink(dds, coef="dex_Dexamethasone_vs_Untreated", type="apeglm",
                  res = resLFC1)
plotMA(resL, ylim = c(-5, 5))

# Gene clustering
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
pheatmap(mat, annotation_col = anno)


# Gene clustering of up-regulated
topVarGenes <- resDF |> 
  filter(reg == "up") |> 
  arrange(padj) |> 
  dplyr::slice(1:50) |> 
  pull(gene)

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
pheatmap(mat, annotation_col = anno)

# Exporting results
write.csv(resDF, file = "results.csv")

library("ReportingTools")
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")
publish(resDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)

################################################################################
# Session information
sessionInfo()









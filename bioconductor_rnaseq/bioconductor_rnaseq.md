# Bioconductor-RNA-seq: RNA-seq workflow: gene-level exploratory analysis and differential expression

In this tutorial, we use a [bioconductor docker images](https://www.bioconductor.org/help/docker/) that includes [asciinema](https://asciinema.org/) functionality to run the [RNA-seq tutorial](https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html). Here we walk through an end-to-end gene-level RNA-seq differential expression workflow using Bioconductor packages. We will start from the FASTQ files, show how these were quantified to the reference transcripts, and prepare gene-level count datasets for downstream analysis. We will perform exploratory data analysis (EDA) for quality assessment and to explore the relationship between samples, perform differential gene expression analysis, and visually explore the results.

In the examples below, `$` indicates the command line prompt within the container.

<!-- blank line -->
----
<!-- blank line -->

## Learning Objectives:
 - pull a [bioconductor docker image](https://hub.docker.com/r/bioconductor/bioconductor_docker) from DockerHub
 - run the bioconductor/bioconductor_docker container via docker
 - perform exploratory data analysis (EDA) on RNA-seq data

## Assignment 
1. Complete the assignment described below.
2. Upload a link to your dockerhub account.
3. Upload a link with screen-cast.

### Prerequisites
* create an asciinema account using email at [asccinema.org](https://asciinema.org/login/new) 
* navigate to the directory: ~/assignments/bioconductor_asciinema
<!-- blank line -->
----
<!-- blank line -->

 ### Assignment Points
|  Rubric        | Points | 
|----------------|-------|
| Screencast     |  -/5  |
| On Time        |  -/5  |
*Total Points: -/10*

## Getting Started

### 1. open docker teminal
<!-- blank line -->
----
<!-- blank line -->

### 2. pull a docker image from DockerHub
```
docker pull dominicklemas/bioconductor_rnaseq:04_2023
```

## 3. boot into container as bash
```
docker run -it -v [path-to-working-directory]:/data dominicklemas/bioconductor_rnaseq:04_2023 bash

example: docker run -it -v c:/Users/djlemas/Documents/rna-seq:/data dominicklemas/bioconductor_rnaseq:04_2023 bash
```
<!-- blank line -->
----
<!-- blank line -->

### 4. link your container to your asciinema.org account by opening the URL in a web browser 
```
$ asciinema auth
```
<!-- blank line -->
----
<!-- blank line -->

### 5. add screen-cast headers and check your docker container
```
$ asciinema rec
$ # Name: 
$ # Date: 
$ # bioconductor:: bioconductor_rnaseq
```
<!-- blank line -->
----
<!-- blank line -->

### 6. start R 
```
R
```
<!-- blank line -->
----
<!-- blank line -->

### 7. Start the Tutorial at 2.3 [Reading in data with tximeta](https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)
``` 
library("airway")
dir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(dir)
list.files(file.path(dir, "quants"))
csvfile <- file.path(dir, "sample_table.csv")
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
coldata
coldata <- coldata[1:2,]
coldata$names <- coldata$Run
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf.gz")
file.exists(coldata$files)
library("tximeta")
se <- tximeta(coldata)
dim(se)
head(rownames(se))
gse <- summarizeToGene(se)
dim(gse)
head(rownames(gse))
```
<!-- blank line -->
----
<!-- blank line -->

### 8. (2.5) SummarizedExperiment
```
library("DESeq2")
data(gse)
gse
assayNames(gse)
head(assay(gse), 3)
colSums(assay(gse))
rowRanges(gse)
seqinfo(rowRanges(gse))
colData(gse)
```
<!-- blank line -->
----
<!-- blank line -->

### 9. (3) The DESeqDataSet object, sample information and the design formula
```
gse$donor
gse$condition
gse$cell <- gse$donor
gse$dex <- gse$condition
levels(gse$dex)
levels(gse$dex) <- c("untrt", "trt")
library("magrittr")
gse$dex %<>% relevel("untrt")
gse$dex
gse$dex <- relevel(gse$dex, "untrt")
library("DESeq2")
dds <- DESeqDataSet(gse, design = ~ cell + dex)
countdata <- round(assays(gse)[["counts"]])
head(countdata, 3)
coldata <- colData(gse)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ cell + dex)
```
<!-- blank line -->
----
<!-- blank line -->
### 10. (4) Exploratory analysis and visualization
```
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)
keep <- rowSums(counts(dds) >= 10) >= 3
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)
log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
library("dplyr")
library("ggplot2")
dds <- estimateSizeFactors(dds)
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
colnames(df)[1:2] <- c("x", "y")  
lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  
  sampleDists <- dist(t(assay(vsd)))
sampleDists
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$dex, dds$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
plotPCA(vsd, intgroup = c("dex", "cell"))
pcaData <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$dex <- dds$dex
gpca.dat$cell <- dds$cell
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = dex, shape = cell)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
mds <- as.data.frame(colData(vsd))  %>%
         cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
mdsPois <- as.data.frame(colData(dds)) %>%
   cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")
```
<!-- blank line -->
----
<!-- blank line -->

### 11. (5)  Differential expression analysis
```
dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, contrast=c("dex","trt","untrt"))
mcols(res, use.names = TRUE)
summary(res)
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
results(dds, contrast = c("cell", "N061011", "N61311"))
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
sum(res$padj < 0.1, na.rm=TRUE)
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])
```
<!-- blank line -->
----
<!-- blank line -->

### 12. (6) Plotting Results
```
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("dex"))
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()
library("apeglm")
resultsNames(dds)
res <- lfcShrink(dds, coef="dex_trt_vs_untrt", type="apeglm")
plotMA(res, ylim = c(-5, 5))
res.noshr <- results(dds, name="dex_trt_vs_untrt")
plotMA(res.noshr, ylim = c(-5, 5))
plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
	 library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
pheatmap(mat, annotation_col = anno)
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
                          mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
                     ylab = "fraction of small p values")
library("IHW")
res.ihw <- results(dds, filterFun=ihw)					 
```
<!-- blank line -->
----
<!-- blank line -->

### 13. (6) Plotting Results
```
library("AnnotationDbi")
library("org.Hs.eg.db")	
columns(org.Hs.eg.db)
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
resOrdered <- res[order(res$pvalue),]
head(resOrdered)	
resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "results.csv")		
library("ReportingTools")
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)	
resGR <- lfcShrink(dds, coef="dex_trt_vs_untrt", type="apeglm", format="GRanges")
resGR	
ens.str <- substr(names(resGR), 1, 15)
resGR$symbol <- mapIds(org.Hs.eg.db, ens.str, "SYMBOL", "ENSEMBL") 
library("Gviz")
window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)
status <- factor(ifelse(resGRsub$padj < 0.05 & !is.na(resGRsub$padj),
                        "sig", "notsig"))
options(ucscChromosomeNames = FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name = "gene ranges", feature = status)
d <- DataTrack(resGRsub, data = "log2FoldChange", baseline = 0,
               type = "h", name = "log2 fold change", strand = "+")
plotTracks(list(g, d, a), groupAnnotation = "group",
           notsig = "grey", sig = "hotpink")
```
<!-- blank line -->
----
<!-- blank line -->

### 14. (8) Removing hidden batch effects
```
library("sva")
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ dex, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
svseq$sv
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$cell, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
 }
 ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + dex
library("RUVSeq")
set <- newSeqExpressionSet(counts(dds))
idx  <- rowSums(counts(set) > 5) >= 2
set  <- set[idx, ]
set <- betweenLaneNormalization(set, which="upper")
not.sig <- rownames(res)[which(res$pvalue > .1)]
empirical <- rownames(set)[ rownames(set) %in% not.sig ]
set <- RUVg(set, empirical, k=2)
pData(set)
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(pData(set)[, i] ~ dds$cell, vertical = TRUE, main = paste0("W", i))
  abline(h = 0)
 }
 ddsruv <- dds
ddsruv$W1 <- set$W_1
ddsruv$W2 <- set$W_2
design(ddsruv) <- ~ W1 + W2 + dex
```
<!-- blank line -->
----
<!-- blank line -->

## 15. (9) Time course experiments
```
library("fission")
data("fission")
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),], 4)
fiss <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("minute","strain"), returnData = TRUE)
fiss$minute <- as.numeric(as.character(fiss$minute))
ggplot(fiss,
  aes(x = minute, y = count, color = strain, group = strain)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()
 resultsNames(ddsTC)
res30 <- results(ddsTC, name="strainmut.minute30", test="Wald")
res30[which.min(res30$padj),]
betas <- coef(ddsTC)
colnames(betas)
topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)
```
<!-- blank line -->
----
<!-- blank line -->

## 16. (10) Session Info
```
sessionInfo()
```
<!-- blank line -->
----
<!-- blank line -->
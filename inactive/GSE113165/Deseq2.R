# ref: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
setwd("/Volumes/HDD24TB/MetaAnalysisProject_Apr2021/inactive/GSE113165")

#### helpers ####
library("stringr")
library("dyplyr")
library("tidyr")
library("magrittr")
library("ggplot2")

#### Importing Data ####
counts <- read.table("countmatrix.sickle.txt", sep = "\t", header = TRUE)
head(counts)
dim(counts)  # 60664    56

# samples <- read.table("pdata.txt", sep = "\t", header = TRUE)
samples <- read.table("sampletable.txt", sep = "\t", header = TRUE)
head(samples, 3)
#                   run age sample.name sex subject susceptibility time  color
# SRR7007949 SRR7007949 old  GSM3098323   m   S3005            low  pre   blue
# SRR7007950 SRR7007950 old  GSM3098324   m   S3005            low post orange
# SRR7007951 SRR7007951 old  GSM3098325   m   S3008            low  pre   blue


#### Constructing DESeqDataSet (DESeqDataSetFromMatrix) ####
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples, 
                              design = ~ subject + time)  # ref = "pre"
dds
# class: DESeqDataSet 
# dim: 60664 56 
# metadata(1): version
# assays(1): counts
# rownames(60664): 1 2 ... 60663 60664
# rowData names(0):
# colnames(56): SRR7007949 SRR7007950 ... SRR7008003 SRR7008004
# colData names(7): run age ... susceptibility times
dim(dds)   # 60664    56
colData(dds)
colnames(dds)
head(assay(dds))
head(counts(dds))
summary(assay(dds))


#### Filtering low counts data ####
# minimal filtering
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
dim(dds)  # 47412    56
dim(colData(dds))  # 56  8

# cf. at least 3 samples with a count of 10 or higher
# keep2 <- rowSums(counts(dds) >= 10) >= 3
# table(keep2)



#### Exploratory analysis and visualization ####
# 1. vst(): Variance Stabilizing Transformation ----
# recommended for n > 30 dataset
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
# blind = FALSE: means that differences between cell lines and treatment 
# (the variables in the design) will not contribute to the expected variance-mean trend of the experiment

# these will be used soon

# cf. rlog()
# n < 30 dataset might go well with rlog()
# rld <- rlog(dds, blind = FALSE)
# head(assay(rld), 3)



# 2. Normalization ---
par(mfrow = c(1, 2))
dds <- estimateSizeFactors(dds)

# 2-1. box-plot ----
boxplot(log2(counts(dds) + 1), col = dds$color, cex.axis = 0.7, 
        las = 1, xlab = "log2(counts)", horizontal = TRUE, main = "Raw counts")
boxplot(log2(counts(dds, normalized = TRUE) + 1), col = dds$color, cex.axis = 0.7, 
        las = 1, xlab = "log2(normalized counts)", horizontal = TRUE, main = "Normalized counts") 

# 2.2 density-plot ----
library(affy)
plotDensity(log2(counts(dds) + 1),  col = dds$color, 
            xlab = "log2(counts)", cex.lab = 0.7, panel.first = grid()) 
plotDensity(log2(counts(dds, normalized = TRUE) + 1), col = dds$color, 
            xlab = "log2(normalized counts)", cex.lab = 0.7, panel.first = grid()) 



# 3. PCA plot ----
# by plotPCA() from DESeq2
plotPCA(vsd, intgroup = c("time"))
plotPCA(vsd, intgroup = c("age"))
plotPCA(vsd, intgroup = c("sex"))

# by ggplot2
pcaData <- plotPCA(vsd, intgroup = c("time", "sex"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
### pca plot rewite
ggplot(pcaData, aes(x = PC1, y = PC2, color = time, shape = age)) +
  geom_point(size = 3) +
  facet_wrap(. ~ sex) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")



# 4. MDS plot ----
sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$time, vsd$subject, sep = " - " )
colnames(sampleDistMatrix) <- NULL

mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = time, shape = sex)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")




#### DESeq2 ####
dds$time <- relevel(as.factor(dds$time), ref = "pre")
dds <- DESeq(dds)
res <- results(dds)  # results(dds, contrast = c("Time", "post", "pre"))
res
# mcols(res, use.names = TRUE)  # metadata

summary(res)  # FDR < 0.1
# out of 47412 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2209, 4.7%
# LFC < 0 (down)     : 2380, 5%
# outliers [1]       : 0, 0%
# low counts [2]     : 23900, 50%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# lower the FDR to 0.05
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
# FALSE  TRUE 
# 16227  3609

sum(res$padj < 0.1, na.rm=TRUE)  # 4589

# subset the res, sort it by the log2 fold change:
resSig <- subset(res, padj < 0.1)
# strongest down-regulation
head(resSig[order(resSig$log2FoldChange),])
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   ENSG00000227496  99.68188       -3.08841  0.254037 -12.15732 5.24500e-34 1.37023e-30
# ENSG00000174502  19.08919       -2.54539  0.312625  -8.14199 3.88844e-16 1.54958e-13
# ENSG00000162763   3.68283       -2.36392  0.449604  -5.25778 1.45801e-07 7.21700e-06
# ENSG00000129988  25.67004       -2.23460  0.333673  -6.69697 2.12782e-11 2.96032e-09
# ENSG00000286986  13.96839       -2.17008  0.271501  -7.99289 1.31809e-15 4.55750e-13
# ENSG00000101251   4.09929       -1.79228  0.342578  -5.23174 1.67924e-07 8.12963e-06

# strongest upregulation
head(resSig[order(-resSig$log2FoldChange),])
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   ENSG00000170525 2633.17262        1.60167  0.208617   7.67755 1.62160e-14 4.53894e-12
# ENSG00000100362    5.48105        1.59309  0.441939   3.60477 3.12426e-04 4.38436e-03
# ENSG00000253619    6.92388        1.58036  0.267468   5.90860 3.45022e-09 2.73136e-07
# ENSG00000152785    3.59613        1.46046  0.613140   2.38193 1.72222e-02 9.14060e-02
# ENSG00000214514   11.17850        1.45929  0.310552   4.69904 2.61392e-06 8.24945e-05
# ENSG00000205215   41.99741        1.43161  0.297997   4.80409 1.55458e-06 5.31269e-05


#### Plotting results ####
# gene that has minimum padj by pre vs post ----
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup = c("time", "sex"))


# Normalized counts with lines connecting cell lines ----
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("time", "sex"),
                         returnData = TRUE)

ggplot(geneCounts, aes(x = time, y = sex, color = time)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = time, y = subject, color = subject, group = sex)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


# Histgram: p values
# with mean normalized count larger than 1
par(mfrow = c(1,1))
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")



#### Annotation ####
library("AnnotationDbi")
library("Homo.sapiens")
# columns(Homo.sapiens)
res$symbol <- mapIds(Homo.sapiens,
                         keys = row.names(res),
                         column = "SYMBOL",
                         keytype = "ENSEMBL" # change
                         )
head(res)
dim(res)  # [1] 47412     7
res <- na.omit(res)
dim(res)  # [1] 17823     7

keep <- res$symbol[!duplicated(res$symbol)]
dup_symbol <- res$symbol[duplicated(res$symbol)]
length(keep)  # [1] 17752
length(dup_symbol)  # [1] 248
res <- subset(res,!duplicated(res$symbol))
dim(res)  # [1] 17752     7

write.table(res, "deseqRes_GSE113165M.txt", sep = "\t")
write.table(res, "deseqRes_GSE113165F.txt", sep = "\t")
dm <- read.table("deseqRes_GSE113165M.txt", sep = "\t", header = TRUE)
df <- read.table("deseqRes_GSE113165F.txt", sep = "\t", header = TRUE)
head(dm)
head(df)


tab = data.frame(log2FC = res$log2FoldChange,
                 log10pVal = -log10(res$pvalue),
                 row.names = rownames(res))
head(tab)

par(mar = c(5, 4, 4, 6), mfrow = c(1,1))
plot(tab,
     pch = 16, cex = 0.6,
     xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue))

## Log2 fold change and adj.p-value cutoffs
lfc = 2
pval = 0.01

signGenes = (abs(tab$log2FC) > lfc & tab$log10pVal > -log10(pval))
## Identifying the selected genes
points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red")
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),
      cex = 0.8, line = 0.5)


################ Explore Results ################
#table(df_cor$p < 0.1)
dfSig <- subset(res, pvalue < 0.01)
dfSigOrdered <- dfSig[order(dfSig$pvalue),]
head(dfSigOrdered)
dim(dfSigOrdered)  # [1] 3262    7

dfSigadj5 <- subset(res, padj < 0.05)
dfSig5Ordered <- dfSigadj5[order(dfSigadj5$padj),]
head(dfSig5Ordered)
dim(dfSig5Ordered)  # [1] 3018    7

dfSigadj01 <- subset(res, padj < 0.01)
dfSig01Ordered <- dfSigadj01[order(dfSigadj01$padj),]
head(dfSig01Ordered)
dim(dfSig01Ordered)  # [1] 1856    7
###################################################
###################################################



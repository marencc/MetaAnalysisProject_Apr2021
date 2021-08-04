# ref: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
setwd("/Users/Emma/Documents/Bioinformatics/DEG/MetaAnalysisProject_Apr2021/inactive/GSE113165")

#### helpers ####
library("stringr")
library("dyplyr")
library("tidyr")
library("magrittr")
library("ggplot2")

#### Importing Data ####
counts <- read.table("counts.txt", sep = "\t", header = TRUE)
head(counts)
colnames(counts) <- sub(".sort.bam", "", colnames(counts))

# samples <- read.table("pdata.txt", sep = "\t", header = TRUE)
samples <- read.table("SraRunTable.txt", sep = ",", header = TRUE)
samples <- samples[,c("Run", "age", "Sample.Name", "sex", "Subject", "susceptibility", "Time")]
head(samples, 3)
# samples$Time <- ifelse(str_detect(samples$Time, "pre"), "pre", "post")

# tiding up sample information
samples$age <- relevel(as.factor(samples$age), ref = "young")
samples$sex <- str_sub(samples$sex, 1, 1)
samples$sex <- relevel(as.factor(samples$sex), ref = "m")
samples$Subject <- paste0("S", samples$Subject)
samples$susceptibility <- str_replace(samples$susceptibility, "Not used in Susceptibility Study", "NA")
samples$susceptibility <- as.factor(samples$susceptibility)
samples$Time <- ifelse(str_detect(samples$Time, "pre"), "pre", "post")
samples$Time <- relevel(as.factor(samples$Time), ref = "pre")
summary(samples)
#  Run               age     Sample.Name        sex      Subject          susceptibility   Time   
# Length:56          young:18   Length:56          m:26   Length:56          high:24        pre :28  
# Class :character   old  :38   Class :character   f:30   Class :character   low :28        post:28  
# Mode  :character              Mode  :character          Mode  :character   NA  : 4 

head(samples, 3)
# Run age Sample.Name sex Subject susceptibility Time
# 1 SRR7007949 old  GSM3098323   m   S3005            low  pre
# 2 SRR7007950 old  GSM3098324   m   S3005            low post
# 3 SRR7007951 old  GSM3098325   m   S3008            low  pre

table(samples$sex, samples$Time)
#   pre post
# m  13   13
# f  15   15


#### Constructing DESeqDataSet (DESeqDataSetFromMatrix) ####
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples, 
                              design = ~ Subject + Time)  # ref = "pre"
dds
dim(dds)
colData(dds)
colnames(dds)
head(assay(dds))
summary(assay(dds))


#### Filtering low counts data ####
# minimal filtering
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
colData(dds) <- colData(dds)[,keep]
dim(dds)
dim(colData(dds))

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
boxplot(log2(counts(dds) + 1),  col = Time, cex.axis = 0.7, 
        las = 1, xlab = "log2(counts)", horizontal = TRUE, main = "Raw counts")
boxplot(log2(counts(dds, normalized = TRUE) + 1),  col = Time, cex.axis = 0.7, 
        las = 1, xlab = "log2(normalized counts)", horizontal = TRUE, main = "Normalized counts") 

# 2.2 density-plot ----
plotDensity(log2(counts(dds) + 1),  col = Time, 
            xlab = "log2(counts)", cex.lab = 0.7, panel.first = grid()) 
plotDensity(log2(counts(dds, normalized = TRUE) + 1), col = Time, 
            xlab = "log2(normalized counts)", cex.lab = 0.7, panel.first = grid()) 



# 3. PCA plot ----
# by plotPCA() from DESeq2
plotPCA(vsd, intgroup = c( "Time", "Subject"))

# by ggplot2
pcaData <- plotPCA(vsd, intgroup = c( "Time", "Subject"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = Time, shape = Subject)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")



# 4. MDS plot ----
sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Time, vsd$Subject, sep = " - " )
colnames(sampleDistMatrix) <- NULL

mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Time, shape = Subject)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")




#### DESeq2 ####
dds$Time <- relevel(as.factor(dds$Time), ref = "pre")
dds <- DESeq(dds)
res <- results(dds)  # results(dds, contrast = c("Time", "post", "pre"))
res
mcols(res, use.names = TRUE)  # metadata

summary(res)  # FDR < 0.1

# lower the FDR to 0.05
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

sum(res$padj < 0.1, na.rm=TRUE)

# subset the res, sort it by the log2 fold change:
resSig <- subset(res, padj < 0.1)
# strongest down-regulation
head(resSig[ order( resSig$log2FoldChange ), ])
# strongest upregulation
head(resSig[ order( -resSig$log2FoldChange ), ])



#### Plotting results ####
# gene that has minimum padj by pre vs post ----
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup = c("Time"))


# Normalized counts with lines connecting cell lines ----
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("Time", "Subject"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Time, y = Subject, color = Subjects)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = Time, y = Subject, color = Subject, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


# Histgram: p values
# with mean normalized count larger than 1
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")



#### Annotation ####
library("AnnotationDbi")
library("Homo.sapiens")
# columns(Homo.sapiens)
res$symbol <- mapIds(Homo.sapiens,
                         keys = row.names(res),
                         column = "SYMBOL",
                         keytype = "SYMBOL" # change
                         )

head(res)
dim(res)  # [1] 54583     4
res <- na.omit(res)
dim(res)  # [1] 31020     4

keep <- res$gene_name[!duplicated(res$gene_name)]
dup_gene_names <- res$gene_name[duplicated(res$gene_name)]
length(keep)  # [1] 30772
length(dup_gene_names)  # [1] 248
res <- res %>% 
    filter(!duplicated(res$gene_name))
dim(res)  # [1] 30772     4

final_res <- data.frame(res[1:3], row.names = res$gene_name)
head(final_res)

write.table(final_res, "adults_vs_fetus_deseq_res.txt", sep = "\t")
d <- read.table("adults_vs_fetus_deseq_res.txt", sep = "\t", header = TRUE)
head(d)


tab = data.frame(log2FC = final_res$log2FoldChange,
                 log10pVal = -log10(final_res$pvalue),
                 row.names = rownames(final_res))
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
dfSig <- subset(final_res, pvalue < 0.01)
dfSigOrdered <- dfSig[order(dfSig$pvalue),]
head(dfSigOrdered)
dim(dfSigOrdered)  # [1] 10325     3

dfSigadj5 <- subset(final_res, padj < 0.05)
dfSig5Ordered <- dfSigadj5[order(dfSigadj5$padj),]
head(dfSig5Ordered)
dim(dfSig5Ordered)  # [1] 11450     3

dfSigadj01 <- subset(final_res, padj < 0.01)
dfSig01Ordered <- dfSigadj01[order(dfSigadj01$padj),]
head(dfSig01Ordered)
dim(dfSig01Ordered)  # [1] 8563    3
###################################################
###################################################











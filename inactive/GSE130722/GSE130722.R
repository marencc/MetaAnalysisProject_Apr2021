setwd("/Users/Emma/Documents/Bioinformatics/DEG/MetaAnalysisProject_Apr2021/inactive/GSE130722/")

counts <- read.table("counts.txt", sep = "\t", header = TRUE)
head(counts)
colnames(counts) <- sub(".sort.bam", "", colnames(counts))

sampleInfo <- read.table("SraRunTable.txt", sep = ",", header = TRUE)
head(sampleInfo) # no gender or id information

# keep columns needed
sampleInfo <- sampleInfo[,c("Run", "Sample.Name", "Time_point")]
head(sampleInfo)

library(GEOquery)
geo <- getGEO("GSE130722", getGPL = FALSE)
length(geo)
gse <- geo[[1]]
head(pData(gse))

sampleInfo <- sampleInfo[match(colnames(counts), sampleInfo$Run),]
rownames(sampleInfo) <- colnames(counts)
head(sampleInfo)  # to check if rowname and Run are the same

sampleInfo$Group <- as.factor(sampleInfo$Group)
col.group <- c("adult"="#4169e1","fetus"="#f08080") # adult: midnight blue, fetus: lightcoral
sampleInfo$color <- col.group[as.vector(sampleInfo$Group)]
head(sampleInfo)
str(sampleInfo)





# https://dputhier.github.io/ASG/practicals/rnaseq_diff_Snf2/rnaseq_diff_Snf2.html
stats.per.sample <- data.frame(t(do.call(cbind, lapply(counts, summary))))
head(stats.per.sample)
dim(stats.per.sample)

## We can now add some columns to the stats per sample
stats.per.sample$libsum <- apply(counts, 2, sum) ## libsum
# Add some percentiles
stats.per.sample$perc05 <- apply(counts, 2, quantile, 0.05)
stats.per.sample$perc10 <- apply(counts, 2, quantile, 0.10)
stats.per.sample$perc90 <- apply(counts, 2, quantile, 0.90)
stats.per.sample$perc95 <- apply(counts, 2, quantile, 0.95)
stats.per.sample$zeros <- apply(counts==0, 2, sum)
stats.per.sample$percent.zeros <- 100*stats.per.sample$zeros/nrow(counts)

# View(stats.per.sample)
stats.per.sample
#             Min. X1st.Qu. Median      Mean X3rd.Qu.    Max.   libsum perc05 perc10 perc90 perc95 zeros percent.zeros
# SRR1554535    0        0      2  657.6492      117 1624827 35896468      0      0 1255.0 2539.0 23215      42.53156
# SRR1554541    0        0      3 1102.3105      185  747136 60167413      0      0 2411.0 5244.8 20057      36.74587
# SRR1554556    0        0      3  830.0649      162 1380381 45307431      0      0 1680.8 3415.0 21041      38.54863
# SRR1554561    0        0      2  684.0431      110 1350225 37337122      0      0 1245.0 2746.9 23685      43.39263
# SRR1554567    0        0      4  981.0084      204  724042 53546382      0      0 2394.8 4842.8 19496      35.71808
# SRR1554568    0        0      2  761.1721      142  534835 41547057      0      0 1822.0 3782.0 22113      40.51261





### Bar plot
par(mar = c(5,6,4,2) + 0.1)
barplot(colSums(counts)/1000000, 
        main = "Total number of reads per sample (million)",
        col = sampleInfo$color, 
        las = 1,  horiz = TRUE,
        cex.names = 0.8,
        xlab = "Million counts")
par(mar = c(5,4,4,2) + 0.1)


### Box plot
# 1.
# boxplot(counts, col=sampleInfo$color, pch=".", 
#         cex.axis=0.5,
#         ylab="Samples", xlab="log2(Counts +1)")
boxplot(log2(counts + 1 ), col=sampleInfo$color, pch=".", 
        cex.axis=0.8,
        ylab="log2(Counts +1)", xlab="Samples")
legend("topright", legend = names(col.group), col = col.group, lwd = 1)

# 2.
library(reshape)
pseudoCount <- log2(counts + 1)
df = melt(pseudoCount, variable_name = "Samples")
head(df)
sampleGroup <- ifelse(df$Samples %in% c("SRR1554535", "SRR1554556", "SRR1554561"), "adult", "fetus")
df = data.frame(df, Condition = sampleGroup)
df$Color <- col.group[as.vector(df$Condition)]

dim(df)

library(ggplot2)
ggplot(df, aes(x = Samples, y = value, fill = Condition)) + geom_boxplot() + xlab("") +
    ylab(expression(log[2](count + 1))) +
    scale_fill_manual(values = sampleInfo$color) +
    labs(fill = "Group") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position = "top")
### title


### Density plot
# 1.
library(affy)
plotDensity(log2(counts + 1), lty = 1, col = sampleInfo$color, lwd = 1)
grid()
legend("topright", legend = names(col.group), col = col.group, lwd = 1)
### title

# 2.
ggplot(df, aes(x = value, colour = Samples, fill = Samples)) + ylim(c(0, 0.25)) +
    geom_density(alpha = 0.2, size = 1) + facet_wrap(~ Condition) +
    theme(legend.position = "top") + xlab(expression(log[2](count + 1)))
### title



### PCA plot (before creating deseq data)
annot = AnnotatedDataFrame(sampleInfo)
expSet = new("ExpressionSet", exprs = as.matrix(pseudoCount), phenoData = annot)

exp_raw <- exprs(expSet)
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])


# PCA plot of the log
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Individual = pData(expSet)$Run,
                     Group = pData(expSet)$Group,
                     RIN = pData(expSet)$RIN,
                     Col = pData(expSet)$color,
                     AGE = pData(expSet)$AGE,
                     Bases = as.character(pData(expSet)$Bases)
)

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(shape = Individual, colour = Group)) +
    ggtitle("PCA plot of the log-transformed expression data") +
    geom_text(aes(label = paste0("age: ", AGE)), nudge_y = 5, size = 3) +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = sampleInfo$color)







# Eliminating undetected genes
# All genes from genome the S. cerevisiae where quantified.
# However only a fraction of them where expressed and some of them
# where to weakly expressed to be detected in any of the sample.
# As a result the count table may contain rows with only zero values (null counts).
# 
# What is the percentage of gene having null counts per sample. Draw a barplot.
# Some genes were not detected in any of the sample. Count their number,
# and delete them from the count.table data.frame.
prop.null <- apply(counts, 2, function(x) 100*mean(x==0))
print(head(prop.null))
# SRR1554535 SRR1554541 SRR1554556 SRR1554561 SRR1554567 SRR1554568 
# 42.53156   36.74587   38.54863   43.39263   35.71808   40.51261


par(mar = c(6,4,4,2) + 0.1)
barplot(prop.null, main="Percentage of null counts per sample", 
        cex.names=0.8, las = 2,
        col=sampleInfo$color, ylab='% of null counts')
par(mar = c(5,4,4,2) + 0.1)




## Some genes were not detected at all in these samples. We will discard them.
counts <- counts[rowSums(counts) > 0, ]
dim(counts)
# [1] 41735     6


library(DESeq2)
deseq_data <- DESeqDataSetFromMatrix(countData = counts, colData = sampleInfo, design = ~ Group)  # ref = "adult"
colData(deseq_data)
colnames(deseq_data)
deseq_data
head(assay(deseq_data))
summary(assay(deseq_data))
# SRR1554535          SRR1554541       SRR1554556          SRR1554561        SRR1554567       SRR1554568      
# Min.   :      0.0   Min.   :     0   Min.   :      0.0   Min.   :      0   Min.   :     0   Min.   :     0.0  
# 1st Qu.:      0.0   1st Qu.:     0   1st Qu.:      0.0   1st Qu.:      0   1st Qu.:     0   1st Qu.:     0.0  
# Median :      2.0   Median :     3   Median :      3.0   Median :      2   Median :     4   Median :     2.0  
# Mean   :    657.6   Mean   :  1102   Mean   :    830.1   Mean   :    684   Mean   :   981   Mean   :   761.2  
# 3rd Qu.:    117.0   3rd Qu.:   185   3rd Qu.:    162.0   3rd Qu.:    110   3rd Qu.:   204   3rd Qu.:   142.0  
# Max.   :1624827.0   Max.   :747136   Max.   :1380381.0   Max.   :1350225   Max.   :724042   Max.   :534835.0


# assay(deseq_data) <- assay(deseq_data) + 1
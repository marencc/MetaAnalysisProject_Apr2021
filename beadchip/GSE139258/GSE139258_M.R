#### GSE139258: ####
#### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139258
#### https://physoc.onlinelibrary.wiley.com/doi/full/10.14814/phy2.14416

# Summary:
# short-term aerobic training -> improvement in skeletal muscle mitochondrial capacity
# -> changes in Tribbles 1 expression
# purpose: understand the molecular mechanisms involved in the adaptation to exercise
# method: performed multiple measurements of mitochondrial capacity both
#        in vivo and ex vivo in lean or overweight individuals.
#        18-day aerobic exercise training.
# results: mitochondrial oxidative respiratory capacity
#          without an appreciable increase in mitochondrial content
#          (robust transcriptome changes)
#          Tribbles pseudokinase 1, TRIB1, might be a potential mediator of
#          the exercise response in human skeletal muscle

# Illumina HumanHT-12 V4.0 expression beadchip
# https://physoc.onlinelibrary.wiley.com/doi/full/10.14814/phy2.14416
# 
# Workflow: http://www.bioconductor.org/packages/release/data/experiment/vignettes/BeadArrayUseCases/inst/doc/BeadArrayUseCases.pdf
# rawdata preprocesinng: https://academic.oup.com/bioinformatics/article/23/16/2183/198810
# GSE139258

#### List of packages required for the workflow ####
# General Bioconductor packages
library(Biobase)
library(oligoClasses)

# # Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

# Analysis and statistics packages
library(limma)
# library(topGO)
# library(ReactomePA)
# library(clusterProfiler)

# Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
# library(RColorBrewer)
# library(pheatmap)

# Formatting/documentation packages
library(rmarkdown)
library(BiocStyle)
library(dplyr)
library(tidyr)

# Helpers:
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
library(devtools)


#### Import raw data as “ExpressionSet” ####
library(GEOquery)
gse <- getGEO(GEO = "GSE139258")[["GSE139258_series_matrix.txt.gz"]]
dim(gse)
gse  # expressionset
head(pData(gse))
exprs(gse)[1:5 , 1:2]

# Data deposited in the GEO database may be either raw or normalized.
head(pData(gse)$data_processing)

# rawdt <- gse
# head(pData(rawdt)$data_processing)


# filter out probes with poor annotation
library(illuminaHumanv4.db)
ids3 <- featureNames(gse)
qual2 <- unlist(mget(ids3, illuminaHumanv4PROBEQUALITY, ifnotfound = NA))
table(qual2)
rem2 <- qual2 == "No match" | qual2 == "Bad" | is.na(qual2)
table(rem2)

# Remove any differences between labs using the removeBatchEffect
# gse.batchcorrect <- removeBatchEffect(log2(exprs(gse)), batch = sites)
# no batchEffect this time
gse_rem2 <- gse[!rem2,]


# if is.na() === TRUE: remove NAs
bad.sample <- colMeans(is.na(exprs(gse_rem2))) > 0.8
bad.gene <- rowMeans(is.na(exprs(gse_rem2))) > 0.5
gse_rem2 <- gse_rem2[!bad.gene,!bad.sample]

# Take a look to make a subset
gse_rem2
varLabels(gse_rem2)


#. Subselect the columns of interest ----
dat <- gse_rem2
head(Biobase::pData(dat))
tail(Biobase::pData(dat))
Biobase::pData(dat) <- Biobase::pData(dat)[, c("title", # individual
                                               "Sex:ch1", # gender
                                               "groups:ch1", # group type
                                               "group:ch1", # group type
                                               "time:ch1")]       # timepoint

#. Rename ----
dat <- dat[,1:12]  # excluding LA & OBS
names(pData(dat)) <- c("individual", "gender", "group-short", "group", "timepoint")
pData(dat)$individual <- factor(
    c("A", "A", "B", "B", "C", "C", "D", "D", "E", "E", "F", "F"))  # male: A,F
                                                                    # female: B,C,D,E

pData(dat)$id <- paste(pData(dat)$individual, pData(dat)$timepoint, sep = ".")

head(pData(dat))
tail(pData(dat))
dim(pData(dat))



# workflow ----
# (also checking if the data is preprocessed/normalized)
boxplot(log2(exprs(dat)), las = 2, cex.names = 0.5, ylab = expression(log [2](intensity)),
        outline = FALSE, main = "Before batch correction")



#### Quality control of the raw data ####
# Take a look
Biobase::exprs(dat)[1:5, 1:5]
head(exprs(dat)[,1:2])


# log2
exp_raw <- log2(Biobase::exprs(dat))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

#. PCA plot of the log ----
pData(dat)
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     individual = pData(dat)$individual,
                     gender = pData(dat)$gender,
                     timepoint = pData(dat)$timepoint,
                     id = pData(dat)$id
)


# all gender × exercise.type
ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(gender ~ .)


#. Create subsets ----
# male
male <- dat[,dat$gender == "Male"]
pData(male)
dim(male)

#. Intensity boxplots of the log2 ----
oligo::boxplot(male, target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               las = 2,
               cex.axis=0.7)

# Perform a differential expression analysis
pData(male)
individual <- factor(Biobase::pData(male)$individual)
timepoint <- factor(Biobase::pData(male)$timepoint)
id <- factor(Biobase::pData(male)$id)

design <- model.matrix(~ 0 + timepoint + individual)
colnames(design) <- c(levels(timepoint), levels(individual)[-1])
rownames(design) <- individual

fit <- lmFit(male, design)
contrast <- makeContrasts(POST-PRE, levels = design)
contr.fit <- eBayes(contrasts.fit(fit, contrast))
# ?eBayes

table <- topTable(contr.fit, coef = 1, number = Inf)
?topTable

volcanoplot(contr.fit, main = "PRE - POST")

# workflow p35: differential expression analysis ----
ids4 <- rownames(male)
chr2 <- mget(ids4, illuminaHumanv4CHR, ifnotfound = NA)
chrloc2 <- mget(ids4, illuminaHumanv4CHRLOC, ifnotfound = NA)
refseq2 <- mget(ids4, illuminaHumanv4REFSEQ, ifnotfound = NA)
entrezid2 <- mget(ids4, illuminaHumanv4ENTREZID, ifnotfound = NA)
symbols2 <- mget(ids4, illuminaHumanv4SYMBOL, ifnotfound = NA)
genename2 <- mget(ids4, illuminaHumanv4GENENAME, ifnotfound = NA)
anno2 <- data.frame(Ill_ID = ids4 , Chr = as.character(chr2),
                    Loc = as.character(chrloc2), RefSeq = as.character(refseq2),
                    Symbol = as.character(symbols2), Name = as.character(genename2),
                    EntrezID = as.numeric(entrezid2))
contr.fit$genes <- anno2
write.csv(contr.fit, file = "res_139258_M.csv")


# Exploring Results
head(table)
nrow(subset(table, P.Value < 0.001))
nrow(subset(table, adj.P.Val < 0.05))

table(table$adj.P.Val < 0.1)
sum(table$adj.P.Val < 0.1, na.rm=TRUE)

resSig <- subset(table, adj.P.Val < 0.1)
head(resSig[order(-resSig$adj.P.Val), ])
head(resSig[order(resSig$adj.P.Val), ])

resOrdered <- table[order(table$adj.P.Val),]
head(resOrdered)

# GSE8479
# Resistance Exercise Reverses Aging in Human Skeletal Muscle
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE8479
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000465

## Summary: ##
# Overall design
# A muscle biopsy was taken from the vastus lateralis muscle of the right or left leg (randomized)
# before exercise or immobilization (young people, N = 26 total) and
# before (N = 25), and after (N = 14), the training period in older adults,
# ~20 cm proximal to the knee joint using a 5mm Bergström biopsy needle.
# The muscle was dissected of fat and connective tissue,
# immediately frozen in liquid nitrogen, and stored at -80°C for subsequent analysis.
# All subjects were required to abstain from strenuous physical activity
# for 48 h prior to the muscle biopsy.

# A muscle biopsy was taken from the vastus lateralis muscle of the right or left leg (randomized)
# before exercise or immobilization (young people, N = 26 total) and
# before (N = 25), and after (N = 14), the training period in older adults
# ref: subjects biopsied after training (N = 14) ... ./supplement/Table_S1.xls

# Sentrix HumanRef-8 Expression BeadChip
# ref: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL2700
# (Sentrix HumanRef-8 Expression BeadChip)

# Annotation:
# method: https://support.bioconductor.org/p/119296/
# another ref: https://www.illumina.com/documents/products/techbulletins/techbulletin_whole_genome_expression.pdf


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
## for convenience, renamed "gse" to "raw"
raw <- getGEO(GEO = "GSE8479")[["GSE8479_series_matrix.txt.gz"]]
gse <- raw

dim(gse)  ## gse == raw[[1]]
# Features  Samples 
# 24354       65

gse
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 24354 features, 65 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM210289 GSM210290 ... GSM210353 (65 total)
# varLabels: title geo_accession ... Sample Group:ch1 (44 total)
# varMetadata: labelDescription
# featureData
# featureNames: GI_10047089-S GI_10047091-S ... trpF (24354 total)
# fvarLabels: ID SequenceSource GB_ACC SPOT_ID
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 17520024 
# Annotation: GPL2700  

head(pData(gse))
exprs(gse)[1:5 , 1:2]


head(fData(gse))  # Will use GB_ACC for annotation
#                       ID         SequenceSource GB_ACC SPOT_ID
# GI_10047089-S GI_10047089-S         RefSeq NM_014332.1      NA
# GI_10047091-S GI_10047091-S         RefSeq NM_013259.1      NA
# GI_10047093-S GI_10047093-S         RefSeq NM_016299.1      NA
# GI_10047099-S GI_10047099-S         RefSeq NM_016303.1      NA
# GI_10047103-S GI_10047103-S         RefSeq NM_016305.1      NA
# GI_10047105-S GI_10047105-S         RefSeq NM_016352.1      NA


### Annotation ###
library(org.Hs.eg.db)
# @ dbfile(org.Hs.eg.db)

# ref: https://support.bioconductor.org/p/119296/
refseq <- sapply(strsplit(fData(gse)[,3], "\\."), "[", 1)
head(refseq)

refseq[is.na(refseq)] <- "NA"
# columns(org.Hs.eg.db)
anno <- lapply(c("ENTREZID","SYMBOL","GENENAME"), 
               function(x) mapIds(org.Hs.eg.db, refseq, x, "ACCNUM"))

annodf <- data.frame(PROBEID = fData(gse)[,1], 
                     ACCNUM = refseq, ENTREZID = anno[[1]],
                     SYMBOL = anno[[2]], GENENAME = anno[[3]],
                     row.names = fData(gse)[,1])

head(annodf)
fData(gse) <- annodf
gse
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 24354 features, 65 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM210289 GSM210290 ... GSM210353 (65 total)
# varLabels: title geo_accession ... Sample Group:ch1 (44 total)
# varMetadata: labelDescription
# featureData
# featureNames: GI_10047089-S GI_10047091-S ... trpF (24354 total)
# fvarLabels: PROBEID ACCNUM ... GENENAME (5 total)
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 17520024 
# Annotation: GPL2700

View(fData(gse))
# fData(gse)

# Data deposited in the GEO database may be either raw or normalized.
head(pData(gse)$data_processing, 1)
# [1] "Raw data (see related publication for normalization procedure)"


# filter out probes with poor annotation (skipped) ----
ids3 <- featureNames(gse)
head(featureNames(gse))
# [1] "GI_10047089-S" "GI_10047091-S" "GI_10047093-S" "GI_10047099-S" "GI_10047103-S"
# [6] "GI_10047105-S"

# library(illuminaHumanv1.db)
    # memo: CHIPNAME: HumanWG6v1
# qual2 <- unlist(mget(ids3, illuminaHumanv1PROBEQUALITY, ifnotfound = NA))

# length(unique(ids3))  # [1] 24354
# length(qual2)  # [1] 24361
    # skipped because this length difference may be causing error for "gse[!rem2,]" below
# table(qual2)
# qual2
#       Bad        Good     Good***    Good****    No match     Perfect  Perfect***  Perfect**** 
#     4040         767          24         160          88       17548         657          1073

# rem2 <- qual2 == "No match" | qual2 == "Bad" | is.na(qual2)
# table(rem2)
# rem2
# FALSE  TRUE 
# 20229  4132  # sum == 24361

# gse_rem2 <- gse[!rem2,]


# if is.na() === TRUE: remove NAs
# bad.sample <- colMeans(is.na(exprs(gse_rem2))) > 0.8
# bad.gene <- rowMeans(is.na(exprs(gse_rem2))) > 0.5
# gse_rem2 <- gse_rem2[!bad.gene,!bad.sample]

# Take a look to make a subset
# gse_rem2
# varLabels(gse_rem2)


#. Subselect the columns of interest ----
dat <- gse
head(Biobase::pData(dat))
tail(Biobase::pData(dat))
Biobase::pData(dat) <- Biobase::pData(dat)[, c("title", # individual
                                               "Gender:ch1")] # gender

#. Rename ----
names(pData(dat)) <- c("individual", "gender")
pData(dat)$individual <- gsub("Subject code ", "", pData(dat)$individual)
pData(dat)$id <- gsub("EB", "", pData(dat)$individual) # individual
pData(dat)$gender <- str_sub(pData(dat)$gender, 0, 1) # gender
pData(dat)$timepoint <- ifelse(str_detect(pData(dat)$individual, "EB"), "pre", "post") # timepoint

library(readxl)
af <- read_xls("./supplements/Table_S1.xls", sheet = 2)  ## who did the biopsy after the training
after <- data.frame(af[1])

keep <- pData(dat) %>% 
    filter(individual %in% after$Subject.Code. | str_detect(individual, "EB"))
keep$individual


nrow(keep)  # extracted only elderly subjects, who did biopsy
            # both before and after exercise training

kp <- pData(dat)$individual %in% keep$individual
dt <- dat[,kp]  # creating the dataset containing only elderly subjects mentioned above
pData(dt) <- keep


head(pData(dt))
tail(pData(dt))
dim(pData(dt))
dt
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 24354 features, 28 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM210315 GSM210316 ... GSM210353 (28 total)
# varLabels: individual gender id timepoint
# varMetadata: labelDescription
# featureData
# featureNames: GI_10047089-S GI_10047091-S ... trpF (24354 total)
# fvarLabels: PROBEID ACCNUM ... GENENAME (5 total)
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 17520024 
# Annotation: GPL2700


# quantile normalization
library(beadarray)
norm <- normaliseIllumina(dt)  
# ?normaliseIllumina

# workflow ----
# (also checking if the dta is preprocessed/normalized)
boxplot(exprs(dt), las = 2, cex.names = 0.5, ylab = expression(log [2](intensity)),
        outline = FALSE, main = "raw")
boxplot(exprs(norm), las = 2, cex.names = 0.5, ylab = expression(log [2](intensity)),
        outline = FALSE, main = "normalized")



#### Quality control of the raw dta ####
# Take a look
Biobase::exprs(dt)[1:5, 1:5]
Biobase::exprs(norm)[1:5, 1:5]


# log2
exp_raw <- log2(Biobase::exprs(norm))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

#. PCA plot of the log ----
pData(norm)
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     individual = pData(norm)$id,
                     gender = pData(norm)$gender,
                     timepoint = pData(norm)$timepoint
)


# all gender × exercise.type
ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(gender ~ .)


#. Create subsets ----
# female
female <- norm[,norm$gender == "F"]
pData(female)
dim(female)
# Features  Samples 
# 24354       16


#. Intensity boxplots of the log2 ----
oligo::boxplot(female, target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               las = 2,
               cex.axis=0.7)
f.log <- log2(Biobase::exprs(female))
oligo::boxplot(f.log, target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               las = 2,
               cex.axis=0.7,
               outline = FALSE)


# Perform a differential expression analysis ----
individual <- factor(Biobase::pData(female)$id)
timepoint <- factor(Biobase::pData(female)$timepoint)

design <- model.matrix(~ 0 + timepoint + individual)
colnames(design) <- c(levels(timepoint), levels(individual)[-1])
rownames(design) <- individual

fit <- lmFit(female, design)
contrast <- makeContrasts(post-pre, levels = design)
contr.fit <- eBayes(contrasts.fit(fit, contrast))

table <- topTable(contr.fit, coef = 1, number = Inf)

volcanoplot(contr.fit, main = "Pre - Post")

# workflow p35: differential expression analysis ----
write.csv(table, file = "res_GSE8479F.csv")


# Exploring Results
head(table)
nrow(subset(table, P.Value < 0.001))
nrow(subset(table, adj.P.Val < 0.05))

table(table$adj.P.Val < 0.05)

resSig <- subset(table, adj.P.Val < 0.05)
head(resSig[order(resSig$adj.P.Val), ])   # ascend
head(resSig[order(-resSig$adj.P.Val), ])  # descend

resOrdered <- table[order(table$adj.P.Val),]
head(resOrdered)



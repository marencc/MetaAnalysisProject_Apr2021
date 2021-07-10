# GSE1295
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
# Journal: https://journals.physiology.org/doi/full/10.1152/japplphysiol.00331.2004
# 


# Affymetrix Human Genome U133A (GPL96)



#### List of packages required for the workflow ####
# General Bioconductor packages
library(Biobase)
library(oligoClasses)

#Annotation packages
library(pd.hugene.1.1.st.v1)
# ??pd.hugene.1.1.st.v1  # This package is to be used in conjunction with the oligo package.
# ls("package:pd.hugene.1.1.st.v1")

library(hugene11sttranscriptcluster.db)
# ??hugene11sttranscriptcluster.db
# ls("package:hugene11sttranscriptcluster.db")

# library(hugene11stprobeset.db)
# hugene11stprobeset.db is based on exon probesets. For a more gene-centric view, use the transcriptcluster
# version of this package.



# Quality control and pre-processing packages
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
# library(rmarkdown)
# library(BiocStyle)
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
if (getOption('timeout') == 60) {
    options(timeout = 500)
}
getOption('timeout')
# [1] 500

# rawdata (skipped) ----
if (!file.exists("data")) {
    dir.create("data")
}
if (length(list.files("data")) < 1) {
    getGEOSuppFiles("GSE1295", makeDirectory = FALSE, baseDir = "data")
}
head(list.files("data"))
# 1st time
# "GSE1295_RAW.tar"
# from the 2nd time
# [1] "GSE1295_RAW.tar" "GSM19136.CEL.gz" "GSM19148.CEL.gz" "GSM19149.CEL.gz" "GSM19150.CEL.gz" "GSM19151.CEL.gz"

if (length(list.files("data")) == 1) {
    RAW.tar <- file.path("data", "GSE1295_RAW.tar")
    untar(RAW.tar, exdir = "data")
    head(list.files("data"))
    tail(list.files("data"))
}

#### Reading in raw data ####
celfiles <- list.files("data", pattern = "CEL.gz$", full = TRUE)
head(celfiles)
chipTypes <- sapply(celfiles, oligo:::getCelChipType, useAffyio=T)
chipTypes  ## Need to extract "HG-U133A"

dir.create("data/HG-U133A")
dir.create("data/HG_U95Av2")

celnames <- list.files("data", pattern = "CEL.gz$")
celnames

for (i in 1:length(celfiles)) {
    celname = celnames[i]
    celfile = celfiles[i]
    if (oligo:::getCelChipType(celfile, useAffyio=T) == "HG-U133A") {
        file.copy(from = celfile, to = paste("data/HG-U133A", celname, sep = "/"), overwrite = TRUE)
        file.remove(celfile)
        } else {    
        file.copy(from = celfile, to = paste("data/HG_U95Av2", celname, sep = "/"), overwrite = TRUE)
        file.remove(celfile)
            }
}

celfiles <- list.files("data/HG-U133A", full = TRUE)

rawData <- read.celfiles(celfiles, verbose = FALSE)
stopifnot(validObject(rawData))

rawData
# ExpressionFeatureSet (storageMode: lockedEnvironment)
# assayData: 506944 features, 24 samples 
# element names: exprs 
# protocolData
# rowNames: GSM19136.CEL.gz GSM19148.CEL.gz ... GSM20663.CEL.gz (24 total)
# varLabels: exprs dates
# varMetadata: labelDescription channel
# phenoData
# rowNames: GSM19136.CEL.gz GSM19148.CEL.gz ... GSM20663.CEL.gz (24 total)
# varLabels: index
# varMetadata: labelDescription channel
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation: pd.hg.u133a  

# storing rawData to gset for convenience ----
gset <- rawData

exprs(gset)[1:3, 1:3]  ## non-log-transformed data
#   GSM3302211_G084_A01_15_DSMT22.CEL.gz GSM3302212_G084_A02_34_06OSN_MEAT-MILK.CEL.gz GSM3302213_G084_A03_61_4007_PROM.CEL.gz
# 1                                103.3                                         100.5                                   157.3
# 2                               5291.8                                        4179.8                                  4762.4
# 3                                145.0                                         113.8                                   162.8


filename <- sampleNames(gset)
head(filename)
# [1] "GSM3302211_G084_A01_15_DSMT22.CEL.gz"          "GSM3302212_G084_A02_34_06OSN_MEAT-MILK.CEL.gz" "GSM3302213_G084_A03_61_4007_PROM.CEL.gz"      
# [4] "GSM3302214_G084_A04_101_8003_PRE_PROM.CEL.gz"  "GSM3302215_G084_A05_109_8011_PRE_PROM.CEL.gz"  "GSM3302216_G084_A06_117_8027_PRE_PROM.CEL.gz"
pData(gset)$filename <- filename

# phenodata (coursera) ----
geoMat <- getGEO("GSE1295")
geoMat
# $`GSE1295-GPL8300_series_matrix.txt.gz`
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 12625 features, 8 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM19150 GSM19165 ... GSM20662 (8 total)
# varLabels: title geo_accession ... data_row_count (24 total)
# varMetadata: labelDescription
# featureData
# featureNames: 1000_at 1001_at ... AFFX-YEL024w/RIP1_at (12625 total)
# fvarLabels: ID GB_ACC ... Gene Ontology Molecular Function (16 total)
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 15347626 
# Annotation: GPL8300 
# 
# $`GSE1295-GPL96_series_matrix.txt.gz`
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 22283 features, 24 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM19136 GSM19148 ... GSM20663 (24 total)
# varLabels: title geo_accession ... data_row_count (24 total)
# varMetadata: labelDescription
# featureData
# featureNames: 1007_s_at 1053_at ... AFFX-TrpnX-M_at (22283 total)
# fvarLabels: ID GB_ACC ... Gene Ontology Molecular Function (16 total)
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 15347626 
# Annotation: GPL96

pD.all <- pData(geoMat[[2]])
head(pD.all, 3)
tail(pD.all, 3)

pD <- pD.all[, c("title", # individual? & supplement
                 "description")] # subject id, gender

#. Renaming columns ----
# names(pD) <- c("individual", "id", "supplement", "gender", "group", "BMI")
head(pD, 3)
#                             individual     id              supplement gender            group   BMI
# GSM3302211          G084_A01_15_DSMT22 DSMT22 placebo / no supplement      M Young (baseline) 25.08
# GSM3302212 G084_A02_34_06OSN_MEAT_MILK  06OSN placebo / no supplement      M Young (baseline) 20.10
# GSM3302213       G084_A03_61_4007_PROM   4007 placebo / no supplement      M Frail (baseline) 23.34



# Tiding up phenodata ----
pD$gender <- as.factor(ifelse(str_detect(pD$title, "-M"), "male", "female")) # timepoint
pD$title <- str_replace(pD$title, "PSS-", "")
pD$title <- str_replace(pD$title, "-s2", "")
pD$id <- str_sub(pD$description, 7, 9)
pD

##.. supplement ----
levels(as.factor(pD$supplement))
# [1] "placebo / no supplement"
pD$supplement <- as.factor(str_sub(pD$supplement, 1, 7))
# Levels: placebo
tail(pD)  ## GSM3302469 supplement == "NA"


##.. id ----
pD$id <- paste("s", pD$id, sep = ".") # individual

##.. checking the final pData ----
head(pD, 3)
#                            individual     id supplement gender            group   BMI timepoint health
# GSM3302213      G084_A03_61_4007_PROM s.4007    placebo      M Frail (baseline) 23.34       pre  frail
# GSM3302215 G084_A05_109_8011_PRE_PROM s.8011    placebo      F Frail (baseline) 31.40       pre  frail
# GSM3302216 G084_A06_117_8027_PRE_PROM s.8027    placebo      F Frail (baseline) 32.60       pre  frail


# Subsetting samples ----
#. excluding "Young" subjects ----
# (they don't have "post" data)
levels(as.factor(pD$group))
# [1] "Frail (after training)"         "Frail (baseline)"               "Healthy older (after training)"
# [4] "Healthy older (baseline)"       "Young (baseline)" 

keep <- grep("Young", pD$group, invert = TRUE)
pD <- pD[keep,]
gset <- gset[,keep]
dim(pD)
# [1] 206   6
dim(gset)
# Features  Samples 
# 1178100      206

#. excluding "protein" induced subjects ----
keep <- grep("protein", pD$supplement, invert = TRUE)
pD <- pD[keep,]
gset <- gset[,keep]
dim(pD)
# [1] 128   6
dim(gset)
# Features  Samples 
# 1178100      128
dim(pD)


gset
# GeneFeatureSet (storageMode: lockedEnvironment)
# assayData: 1178100 features, 128 samples 
# element names: exprs 
# protocolData
# rowNames: GSM3302213_G084_A03_61_4007_PROM.CEL.gz
# GSM3302215_G084_A05_109_8011_PRE_PROM.CEL.gz ...
# GSM3302469_G087_H07_238_PAM_POST_STOUD.CEL.gz (128 total)
# varLabels: exprs dates
# varMetadata: labelDescription channel
# phenoData
# rowNames: GSM3302213_G084_A03_61_4007_PROM.CEL.gz
# GSM3302215_G084_A05_109_8011_PRE_PROM.CEL.gz ...
# GSM3302469_G087_H07_238_PAM_POST_STOUD.CEL.gz (128 total)
# varLabels: index filename
# varMetadata: labelDescription channel
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation: pd.hugene.1.1.st.v1

levels(as.factor(pD$group))
# [1] "Frail (after training)" "Frail (baseline)" "Healthy older (after training)" "Healthy older (baseline)" 
levels(as.factor(pD$supplement))
# [1] "placebo / no supplement"




#. Merging Target Data & gset ----
sampleNames(gset) <- str_sub(sampleNames(gset), 1, 10)
pD <- pD[sampleNames(gset),]
pData(gset) <- pD
gset
# GeneFeatureSet (storageMode: lockedEnvironment)
# assayData: 1178100 features, 128 samples 
# element names: exprs 
# protocolData
# rowNames: GSM3302213 GSM3302215 ... GSM3302469 (128 total)
# varLabels: exprs dates
# varMetadata: labelDescription channel
# phenoData
# rowNames: GSM3302213 GSM3302215 ... GSM3302469 (128 total)
# varLabels: individual id ... health (8 total)
# varMetadata: labelDescription channel
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation: pd.hugene.1.1.st.v1



#### new (coursera) ####
#### checking data intensities
oligo::boxplot(gset, target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               las = 2,
               cex.axis=0.7)


#### normalization ####
normData <- oligo::rma(gset, target = "core")
normData  ## rma returns "ExpressionSet"
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 33297 features, 128 samples 
# element names: exprs 
# protocolData
# rowNames: GSM3302213 GSM3302215 ... GSM3302469 (128 total)
# varLabels: exprs dates
# varMetadata: labelDescription channel
# phenoData
# rowNames: GSM3302213 GSM3302215 ... GSM3302469 (128 total)
# varLabels: individual id ... health (8 total)
# varMetadata: labelDescription channel
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation: pd.hugene.1.1.st.v1


# checking if data is normalized
oligo::boxplot(normData, target = "core",
               main = "Boxplot of log2-intensitites for the normalized data",
               las = 2,
               cex.axis=0.7)

exprs(normData[1:3, 1:3])  ## now exprs() returns log-transformed data
#         GSM3302213 GSM3302215 GSM3302216
# 7892501   3.759487   2.924124   4.010053
# 7892502   4.233000   3.765880   3.835612
# 7892503   1.695574   1.912988   1.922395


#### Quality control of the raw data ####
#. Create subsets ----
# male
male <- normData[,normData$gender == "M"]
pData(male)
dim(male)
# Features  Samples 
#    33297       90  


#. PCA analysis ----
myexp <- Biobase::exprs(male)
PCA <- prcomp(t(myexp), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     individual = head(pData(male))$id,
                     timepoint = pData(male)$timepoint,
                     health = pData(male)$health
)

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(. ~ health) +
    ggtitle("PCA_Norm timepoint(S) & individual(C)") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5))+
    coord_fixed(ratio = sd_ratio) +
    theme(legend.position = "none")






#### Filtering based on intensity ####
# filter out lowly expressed genes
# “soft” intensity based filtering here, since this is recommended by the limma
#. Histogram ----
medians <- rowMedians(Biobase::exprs(male))

hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

#set a manual threshold to exclude genes, but not too many
# man_threshold <- round(mean(medians))
man_threshold <- 4

hist_res <- hist(medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)

#. Exclude genes that below the threshold ----
# Transcripts that do not have intensities larger than the threshold in
# at least as many arrays as the smallest experimental group are excluded.
# get a list:
no_of_samples <- table(pData(male)$timepoint)
no_of_samples

samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(Biobase::exprs(male), 1,
                           function(x){
                               sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)


manfiltered <- subset(male, idx_man_threshold)

dim(exprs(male))  # 33297    90
dim(exprs(manfiltered))  # 21922    90


#### Differential Expression nalysis ####
individual <- factor(Biobase::pData(manfiltered)$id)
length(individual) ## 90
timepoint <- factor(Biobase::pData(manfiltered)$timepoint)
health <- factor(Biobase::pData(manfiltered)$health)

str(pData(manfiltered))
pData(manfiltered)$id
length(unique(pData(manfiltered)$id))
length(unique(pData(manfiltered)$id))


# normData (without supplement, health)
design <- model.matrix(~ 0 + timepoint + individual)
colnames(design) <- c(levels(timepoint), levels(individual)[-1])
rownames(design) <- individual


fit <- lmFit(manfiltered, design)
contrast <- makeContrasts(post-pre, levels = design)
contr.fit <- eBayes(contrasts.fit(fit, contrast))

table <- topTable(contr.fit, coef = 1, number = Inf)

volcanoplot(contr.fit, main = "Pre - Post")

# Write out the result ----
# write.csv(table, file = "res_GSE1295Mplacebo.csv")


# Exploring Results
head(table)
nrow(subset(table, P.Value < 0.001))  # 256
nrow(subset(table, adj.P.Val < 0.05))  # 151

resSig <- subset(table, adj.P.Val < 0.05)
head(resSig[order(resSig$adj.P.Val), ])
head(resSig[order(-resSig$adj.P.Val), ])

resOrdered <- table[order(table$adj.P.Val),]
head(resOrdered)





#### Annotation of the transcript clusters ####
# columns(hugene11sttranscriptcluster.db)
# keytypes(hugene11sttranscriptcluster.db)
# head(keys(hugene11sttranscriptcluster.db, keytype="GOALL"))
anno <- AnnotationDbi::select(hugene11sttranscriptcluster.db,
                              keys = (featureNames(manfiltered)),
                              columns = c("SYMBOL", "GENENAME", "ACCNUM"),
                              keytype = "PROBEID")
ls("package:hugene11sttranscriptcluster.db")
table(is.na(anno))

anno <- subset(anno, !is.na(SYMBOL))
dim(anno)

#. Removing multiple mappings ----
# group by PROBEID
anno_grouped <- group_by(anno, PROBEID)

# summarize
anno_summarized <-
    dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized)

# filter to extract cluster multiple gene that map to multiple gene symbols
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
head(anno_filtered)

probe_stats <- anno_filtered
nrow(probe_stats) # clusters that map to multiple gene symbols → remove

# remove above IDs(probe_stats)
ids_to_exlude <- (featureNames(manfiltered) %in% probe_stats$PROBEID)
table(ids_to_exlude)
final <- subset(manfiltered, !ids_to_exlude)
validObject(final)

# also exclude them from the feature data anno
head(anno)

# generate a column PROBEID in fData(final) and
# assign the row names of fData(final) to it:
fData(final)$PROBEID <- rownames(fData(final))

# left-join fData(final) with anno
fData(final) <- left_join(fData(final), anno)

rownames(fData(final)) <- fData(final)$PROBEID
validObject(final)

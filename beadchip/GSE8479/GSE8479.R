# GSE8479
# Resistance Exercise Reverses Aging in Human Skeletal Muscle
# n = 14 (m6 f8)
# older
# GPL2700
# Sentrix HumanRef-8 Expression BeadChip

library(Biobase)
library(oligoClasses)
library(oligo)
library(arrayQualityMetrics)
library(limma)
# library(gplots)
library(ggplot2)
# library(geneplotter)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
library(devtools)


#### loading the data ####
library(GEOquery)
if (!file.exists("data")) {
    dir.create("data")
    }

data <- getGEO("GSE8479", destdir = ".")
length(data)  # 1
gse <- data[[1]]

# exploring the data ----
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

dim(gse)
# Features  Samples 
# 24354       65

head(pData(gse))
exprs(gse)[1:5 , 1:2]

head(fData(gse))  # GB_ACC will be used for annotation
#                       ID         SequenceSource GB_ACC SPOT_ID
# GI_10047089-S GI_10047089-S         RefSeq NM_014332.1      NA
# GI_10047091-S GI_10047091-S         RefSeq NM_013259.1      NA
# GI_10047093-S GI_10047093-S         RefSeq NM_016299.1      NA
# GI_10047099-S GI_10047099-S         RefSeq NM_016303.1      NA
# GI_10047103-S GI_10047103-S         RefSeq NM_016305.1      NA
# GI_10047105-S GI_10047105-S         RefSeq NM_016352.1      NA


# checking if the loaded data is already normalized ----
head(pData(gse)$data_processing, 1)
# [1] "Raw data (see related publication for normalization procedure)"



#### Annotation ####
library(org.Hs.eg.db)  # database stored at: dbfile(org.Hs.eg.db)

# ref: https://support.bioconductor.org/p/119296/
refseq <- sapply(strsplit(fData(gse)[,3], "\\."), "[", 1)
head(refseq)
# [1] "NM_014332" "NM_013259" "NM_016299" "NM_016303" "NM_016305" "NM_016352"

sum(is.na(refseq))  # 4
refseq[is.na(refseq)] <- "NA"
# columns(org.Hs.eg.db)
anno <- lapply(c("ENTREZID","SYMBOL","GENENAME"), 
               function(x) mapIds(org.Hs.eg.db, refseq, x, "ACCNUM"))

annodf <- data.frame(PROBEID = fData(gse)[,1], 
                     ACCNUM = refseq,
                     ENTREZID = anno[[1]],
                     SYMBOL = anno[[2]],
                     GENENAME = anno[[3]],
                     row.names = fData(gse)[,1])

head(annodf, 3)
#                     PROBEID    ACCNUM ENTREZID SYMBOL                                      GENENAME
# GI_10047089-S GI_10047089-S NM_014332    23676   SMPX                 small muscle protein X-linked
# GI_10047091-S GI_10047091-S NM_013259    29114 TAGLN3                                  transgelin 3
# GI_10047093-S GI_10047093-S NM_016299    51182 HSPA14 heat shock protein family A (Hsp70) member 14

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


# filtering out probes with poor annotation (skipped) ----
# contents are logged in "memo.Rmd"



#### Tiding up the data ####
#. Subselect the columns of interest ----
dat <- gse
pD.all <- pData(dat)
head(pD.all, 3)
tail(pD.all, 3)
pD <- pD.all[, c("title", # individual
                 "Gender:ch1")] # gender

#. Rename ----
head(pD, 3)
#                 individual gender
# GSM210289 Subject code A30    F;F
# GSM210290 Subject code A44    M;M
# GSM210291 Subject code A41    M;M
names(pD) <- c("individual", "gender")
pD$individual <- gsub("Subject code ", "", pD$individual)
pD$id <- gsub("EB", "", pD$individual) # individual
pD$gender <- str_sub(pD$gender, 0, 1) # gender
pD$timepoint <- ifelse(str_detect(pD$individual, "EB"), "pre", "post") # timepoint


# extracting only older subjects have pre & post data ----
library(readxl)
idfile <- read_xls("./supplements/Table_S1.xls", sheet = 2)  # who did the biopsy after the training
                                                             # downloaded from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000465#s5
                                                             # @ Supporting Information
subjectids_after <- data.frame(idfile[1])

keep <- pD %>% 
    filter(individual %in% subjectids_after$Subject.Code. | str_detect(individual, "EB"))
keep$individual
# [1] "A2"    "A5"    "A17"   "A14"   "A25"   "A6"    "A20"   "A10"   "A11"   "A4"    "A8"    "A22"   "A12"  
# [14] "A19"   "A4EB"  "A2EB"  "A22EB" "A5EB"  "A11EB" "A17EB" "A20EB" "A10EB" "A25EB" "A6EB"  "A8EB"  "A12EB"
# [27] "A14EB" "A19EB"

nrow(keep)  # 28

kp <- pD$individual %in% keep$individual
dt <- dat[,kp]  # creating the dataset containing only elderly subjects mentioned above
pData(dt) <- keep
dim(dt)
# Features  Samples 
# 24354       28 

head(pData(dt), 3); tail(pData(dt), 3)
# > head(pData(dt), 3)
# individual gender  id timepoint
# GSM210315         A2      F  A2      post
# GSM210316         A5      F  A5      post
# GSM210317        A17      M A17      post
# > tail(pData(dt), 3)
# individual gender  id timepoint
# GSM210351      A12EB      F A12       pre
# GSM210352      A14EB      M A14       pre
# GSM210353      A19EB      M A19       pre



#### Normalization ####
library(beadarray)
eset <- dt
norm <- normaliseIllumina(eset)  # quantile normalization
# ?normaliseIllumina
# cf. original study uses normalizeBetweenArrays() of limma.
#     checked with norm_limma <- normalizeBetweenArrays(exprs(norm_limma), method = "quantile").
#     the boxplots (exprs(norm) & exprs(norm_limma)) are identical.

# comparing pre/post normalization
par(mfrow = c(1,2))
boxplot(log2(exprs(dt)), las = 2, cex.names = 0.5, ylab = expression(log [2](intensity)),
        outline = FALSE, main = "raw")
boxplot(log2(exprs(norm)), las = 2, cex.names = 0.5, ylab = expression(log [2](intensity)),
        outline = FALSE, main = "normalized")



#### Quality control ####
exprs(dt)[1:5, 1:5]
exprs(norm)[1:5, 1:5]

# log2
exp_raw <- log2(exprs(norm))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

#. PCA analysis ----
head(pData(norm), 3)
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     individual = pData(norm)$id,
                     gender = pData(norm)$gender,
                     timepoint = pData(norm)$timepoint
                     )
ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(gender ~ .)



#### Filtering based on intensity (skipped) ####
# filter out lowly expressed genes
# “soft” intensity based filtering here, since this is recommended by the limma
#. Histogram ----
medians <- rowMedians(log2(Biobase::exprs(norm)))
par(mfrow=c(1,1))
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")
# decided not to filter out the probes because there are not too many genes
#. below are skipped ----

# set a manual threshold to exclude genes, but not too many
# man_threshold <- round(mean(medians))
# man_threshold <- 4
# abline(v = man_threshold, col = "coral4", lwd = 2)
# 
# #. Exclude genes that below the threshold (skipped)----
# # Transcripts that do not have intensities larger than the threshold in
# # at least as many arrays as the smallest experimental group are excluded.
# # get a list:
# no_of_samples <- table(paste0(pDnorm$timepoint, "_", pDnorm$gender))
# no_of_samples
# 
# samples_cutoff <- min(no_of_samples)
# idx_man_threshold <- apply(Biobase::exprs(norm), 1,
#                            function(x){
#                                sum(x > man_threshold) >= samples_cutoff})
# table(idx_man_threshold)
# # FALSE  TRUE 
# # 11508 21789
# 
# manfiltered <- subset(norm, idx_man_threshold)
# dim(exprs(norm))  # 33297    70
# dim(exprs(manfiltered))  # 21789    70






##### Differential Expression nalysis ####
head(pData(norm), 3)
individual <- Biobase::pData(norm)$id
timepoint <- Biobase::pData(norm)$timepoint
gender <- Biobase::pData(norm)$gender

# normData (without supplement, health)
male <- individual[gender == "M"]
t_male <- timepoint[gender == "M"]
design_male <- model.matrix(~ 0 + t_male + male)
colnames(design_male)[1:2] <- c("post", "pre")
rownames(design_male) <- male

female <- individual[gender == "F"]
t_female <- timepoint[gender == "F"]
design_female <- model.matrix(~ 0 + t_female + female)
colnames(design_female)[1:2] <- c("post", "pre")
rownames(design_female) <- female

# inspect the design matrices
design_male[1:5, 1:5]; design_female[1:5, 1:5]
design_male[1:5, 1:5]
# > design_male[1:5, 1:5]
#       post pre maleA17 maleA19 maleA20
# A17    1   0       1       0       0
# A14    1   0       0       0       0
# A25    1   0       0       0       0
# A20    1   0       0       0       1
# A22    1   0       0       0       0
# > design_female[1:5, 1:5]
#       post pre femaleA11 femaleA12 femaleA2
# A2     1   0         0         0        1
# A5     1   0         0         0        0
# A6     1   0         0         0        0
# A10    1   0         0         0        0
# A11    1   0         1         0        0


# Contrasts and hypotheses tests: all genes ----
contrast_matrix_male <- makeContrasts(post-pre, levels = design_male)
fit_male <- eBayes(
    contrasts.fit(
        lmFit(
            norm[,gender == "M"],
            design = design_male),
        contrast_matrix_male))

contrast_matrix_female <- makeContrasts(post-pre, levels = design_female)
fit_female <- eBayes(
    contrasts.fit(
        lmFit(
            norm[,gender == "F"],
            design = design_female),
        contrast_matrix_female))


#### Extracting results: topTable() ####
library(RColorBrewer)
table_male <- topTable(fit_male, number = Inf)
head(table_male)
par(mfrow = c(1,2))
hist(table_male$P.Value, col = brewer.pal(3, name = "Set2")[1], xlab = "male p-values")
hist(table_male$adj.P.Val, col = brewer.pal(3, name = "Set2")[2], xlab = "male adj.p-values")

table_female <- topTable(fit_female, number = Inf)
head(table_female)
hist(table_female$P.Value, col = brewer.pal(3, name = "Set2")[1], xlab = "female p-values")
hist(table_female$adj.P.Val, col = brewer.pal(3, name = "Set2")[2], xlab = "female p-values")


#. Multiple testing FDR, and comparison with results from the original paper
nrow(subset(table_male, P.Value < 0.05))  # 5464 
tail(subset(table_male, P.Value < 0.05))
nrow(subset(table_male, adj.P.Val < 0.05))  # 161
nrow(subset(table_male, adj.P.Val < 0.01)) # 0
tail(subset(table_male, adj.P.Val < 0.01))

nrow(subset(table_female, P.Value < 0.05))  # 6416
tail(subset(table_female, P.Value < 0.05))
nrow(subset(table_female, adj.P.Val < 0.05))  # 1501
nrow(subset(table_female, adj.P.Val < 0.01))  # 92
tail(subset(table_female, adj.P.Val < 0.01))


#### writing out the result ####
write.csv(table_male, file = "res_GSE8479Mnew.csv")
write.csv(table_female, file = "res_GSE8479Fnew.csv")

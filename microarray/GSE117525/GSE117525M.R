# GSE117525
# Expression of protocadherin gamma in skeletal muscle tissue is associated with
# age and muscle weakness
# journal: https://onlinelibrary.wiley.com/doi/full/10.1002/jcsm.12099
# geo: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117525
# pubmed(abst): https://pubmed.ncbi.nlm.nih.gov/27239416/

## treatment_protocol_ch1
# GSM3302211 Training for both frail and healthy older subjects and consisted of
# progressive full‐body resistance‐type exercise training. However,
# the frail older group had training sessions twice per week,
# whereas the healthy older group trained three times per week.
# In addition, subjects took a protein or control drink for the duration of the study.
# The healthy older group received a 15g portion of milk protein or control supplement at breakfast.
# The frail older group received a similar drink containing 15g supplement drink (milk protein or control)
# at breakfast and lunch.

# Affymetrix Human Gene 1.1 ST Array (GPL20880)
# Annotation: pd.hugene.1.1.st.v1


#### List of packages required for the workflow ####
# General Bioconductor packages
library(Biobase)
library(oligoClasses)

#Annotation packages
library(pd.hugene.1.1.st.v1)
library(hugene11sttranscriptcluster.db)
# ??hugene11sttranscriptcluster.db
# ls("package:hugene11sttranscriptcluster.db")

# Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)
library(AnnotationDbi)

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
getGEOSuppFiles("GSE117525", makeDirectory = FALSE, baseDir = "data")
}
head(list.files("data"))
# 1st time
# "GSE117525_RAW.tar"
# from the 2nd time
# [1] "GPL20880_hugene11st_Hs_ENTREZG_desc.annot.txt.gz"
# [2] "GPL20880_hugene11st_Hs_ENTREZG_desc.txt.gz"      
# [3] "GPL20880_hugene11st_Hs_ENTREZG_mapping.txt.gz"   
# [4] "GPL20880_hugene11st_Hs_ENTREZG_probe_tab.txt.gz" 
# [5] "GPL20880_hugene11st_Hs_ENTREZG.cdf.gz"           
# [6] "GSE117525_RAW.tar"

if (length(list.files("data")) == 1) {
RAW.tar <- file.path("data", "GSE117525_RAW.tar")
untar(RAW.tar, exdir = "data")
head(list.files("data"))
tail(list.files("data"))
}

#### Reading in raw data ####
celfiles <- list.files("data", pattern = "CEL.gz$", full = TRUE)
head(celfiles)

rawData <- oligo::read.celfiles(celfiles)
stopifnot(validObject(rawData))

rawData
# GeneFeatureSet (storageMode: lockedEnvironment)
# assayData: 1178100 features, 259 samples 
# element names: exprs 
# protocolData
# rowNames: GSM3302211_G084_A01_15_DSMT22.CEL.gz
# GSM3302212_G084_A02_34_06OSN_MEAT-MILK.CEL.gz ...
# GSM3302469_G087_H07_238_PAM_POST_STOUD.CEL.gz (259 total)
# varLabels: exprs dates
# varMetadata: labelDescription channel
# phenoData
# rowNames: GSM3302211_G084_A01_15_DSMT22.CEL.gz
# GSM3302212_G084_A02_34_06OSN_MEAT-MILK.CEL.gz ...
# GSM3302469_G087_H07_238_PAM_POST_STOUD.CEL.gz (259 total)
# varLabels: index
# varMetadata: labelDescription channel
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation: pd.hugene.1.1.st.v1 

# storing rawData to gset for convenience ----
gset <- rawData

exprs(gset)[1:3, 1:3]  ## non-log-transformed data
#   GSM3302211_G084_A01_15_DSMT22.CEL.gz GSM3302212_G084_A02_34_06OSN_MEAT-MILK.CEL.gz GSM3302213_G084_A03_61_4007_PROM.CEL.gz
# 1                                103.3                                         100.5                                   157.3
# 2                               5291.8                                        4179.8                                  4762.4
# 3                                145.0                                         113.8                                   162.8


filenames <- sampleNames(gset)
head(filenames)
# [1] "GSM3302211_G084_A01_15_DSMT22.CEL.gz"          "GSM3302212_G084_A02_34_06OSN_MEAT-MILK.CEL.gz" "GSM3302213_G084_A03_61_4007_PROM.CEL.gz"      
# [4] "GSM3302214_G084_A04_101_8003_PRE_PROM.CEL.gz"  "GSM3302215_G084_A05_109_8011_PRE_PROM.CEL.gz"  "GSM3302216_G084_A06_117_8027_PRE_PROM.CEL.gz"
pData(gset)$filenames <- filenames

# phenodata (coursera) ----
# library(GEOquery)
geoMat <- getGEO("GSE117525")
geoMat[1]
# $GSE117525_series_matrix.txt.gz
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 19654 features, 259 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM3302211 GSM3302212 ... GSM3302469 (259 total)
# varLabels: title geo_accession ... weight (kg):ch1 (51 total)
# varMetadata: labelDescription
# featureData
# featureNames: 100009676_at 10000_at ... 9_at (19654 total)
# fvarLabels: ID Description SPOT_ID
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 27239416 
# Annotation: GPL20880

pD.all <- pData(geoMat[[1]])
head(pD.all, 3)
tail(pD.all, 3)
pD <- pD.all[, c("title", # individual? & supplement
                 "subjectid:ch1", # subject id
                 "supplement:ch1", # drink type including placebo
                 "Sex:ch1", # gender
                 "received training:ch1", # pre/post
                 "group:ch1" # yng/old
                 )] # just for curiosity

#. Renaming columns ----
names(pD)[2:6] <- c("subjectid", "supplement", "gender", "training", "agegroup")
head(pD, 3)
#                                  title subjectid              supplement gender training         agegroup age (yrs):ch1
# GSM3302211          G084_A01_15_DSMT22    DSMT22 placebo / no supplement      M       no Young (baseline)            21
# GSM3302212 G084_A02_34_06OSN_MEAT_MILK     06OSN placebo / no supplement      M       no Young (baseline)            22
# GSM3302213       G084_A03_61_4007_PROM      4007 placebo / no supplement      M       no Frail (baseline)            83
# 
#                    height (m):ch1 weight (kg):ch1 bmi (kg/m2):ch1
# GSM3302211           1.94           94.40           25.08
# GSM3302212           1.84           68.20           20.10
# GSM3302213           1.63           62.00           23.34

dim(pD)  # 259   6
length(unique(pD$subjectid))  # 187



# Tiding up phenodata ----
head(pD)
##. id ----
pD$subjectid <- str_replace(pD$subjectid, "pre", "")
pD$subjectid <- gsub(" ", "", pD$subjectid)
pD$subjectid <- paste("S", pD$subjectid, sep = ".") # individual
sum(table(pD$subjectid) == 1)  # 115: these are subjects that have no post data
                               # (number matches to the study)

# excluding the above 115 subjects ----
# checking if all subject ids are coupled
# table(pD$subjectid)
# below subjects are not coupled (match with GEO info)
pD %>%
    group_by(subjectid) %>%
    filter(n() == 1)
# # A tibble: 115 x 8
# Groups:   subjectid [115]
# title                       subjectid supplement gender training agegroup         timepoint health 
#   <chr>                       <chr>     <fct>      <chr>  <chr>    <chr>            <fct>     <fct>  
# 1 G084_A01_15_DSMT22          S.DSMT22  placebo    M      no       Young (baseline) NA        healthy
# 2 G084_A02_34_06OSN_MEAT_MILK S.06OSN   placebo    M      no       Young (baseline) NA        healthy
# 3 G084_A03_61_4007_PROM       S.4007    placebo    M      no       Frail (baseline) NA        frail  
# 4 G084_A11_303_4058_PROM      S.4058    placebo    M      no       Frail (baseline) NA        frail  
# 5 G084_B01_16_DSMT23          S.DSMT23  placebo    M      no       Young (baseline) NA        healthy
# 6 G084_B02_38_11ETK_MEAT_MILK S.11ETK   placebo    M      no       Young (baseline) NA        healthy
# 7 G084_B03_62_4010_PROM       S.4010    placebo    M      no       Frail (baseline) NA        frail  
# 8 G084_B11_304_4065_PROM      S.4065    placebo    M      no       Frail (baseline) NA        frail  
# 9 G084_C01_17_DSMT24          S.DSMT24  placebo    M      no       Young (baseline) NA        healthy
# 10 G084_C02_40_08ACN_MEAT_MILK S.08ACN  placebo    M      no       Young (baseline) NA        healthy
# # … with 105 more rows

table(table(pD$subjectid) == 2)
# FALSE  TRUE 
# 115    72

# keeping coupled subjects
prepost_subjects <- pD %>%
    group_by(subjectid) %>%
    filter(n() == 2)

keep <- pD$subjectid %in% prepost_subjects$subjectid
table(keep)

pD <- pD[keep,]
gset <- gset[,keep]
dim(pD)  # 144   6
dim(gset)
# Features  Samples 
# 1178100      144


#. timepoint (pre/post) ----
pD$timepoint <- as.factor(ifelse(str_detect(pD$title, "PRE"), "pre", 
                                 ifelse(str_detect(pD$title, "POST"), "post", "NA"))) # timepoint
pD$timepoint <- relevel(pD$timepoint, ref = "pre")  # setting "pre" as a reference level
levels(pD$timepoint)
# [1] "pre"  "post"


#. Frail/Healthy ----
unique(pD$agegroup)
# [1] "Frail (baseline)"               "Frail (after training)"         "Healthy older (baseline)"      
# [4] "Healthy older (after training)"
pD$health <- as.factor(ifelse(str_detect(pD$agegroup, "Frail"), "frail", "healthy"))
pD$health <- relevel(pD$health, ref = "healthy")
levels(pD$health)
# [1] "healthy" "frail"


#. supplement ----
# excluding "protein" induced subjects
keep <- grep("placebo", pD$supplement)
pD <- pD[keep,]
gset <- gset[,keep]
dim(pD)  # 70  8
dim(gset)
# Features  Samples 
# 1178100      70
levels(as.factor(pD$supplement))
# [1] "placebo / no supplement"

# renaming the label
pD$supplement <- ifelse(str_detect(pD$supplement, "placebo"), "placebo", "NA")
levels(as.factor(pD$supplement))
# [1] "placebo"


# checking the final pData ----
head(pD, 3)
#                                 title subjectid supplement gender training         agegroup timepoint  health
# GSM3302215 G084_A05_109_8011_PRE_PROM    S.8011    placebo      F      yes Frail (baseline)       pre  frail
# GSM3302216 G084_A06_117_8027_PRE_PROM    S.8027    placebo      F      yes Frail (baseline)       pre  frail
# GSM3302218 G084_A08_155_8082_PRE_PROM    S.8082    placebo      M      yes Frail (baseline)       pre  frail

str(pD)
# 'data.frame':	70 obs. of  8 variables:
# $ title     : chr  "G084_A05_109_8011_PRE_PROM" "G084_A06_117_8027_PRE_PROM" "G084_A08_155_8082_PRE_PROM" "G084_A09_165_8093_PRE_PROM" ...
# $ subjectid : chr  "S.8011" "S.8027" "S.8082" "S.8093" ...
# $ supplement: Factor w/ 1 level "placebo": 1 1 1 1 1 1 1 1 1 1 ...
# $ gender    : chr  "F" "F" "M" "F" ...
# $ training  : chr  "yes" "yes" "yes" "yes" ...
# $ agegroup  : chr  "Frail (baseline)" "Frail (baseline)" "Frail (baseline)" "Frail (baseline)" ...
# $ timepoint : Factor w/ 2 levels "pre","post": 1 1 1 1 2 2 2 2 1 1 ...
# $ health    : Factor w/ 2 levels "healthy","frail": 2 2 2 2 2 2 2 2 2 2 ...


#. Merging Target Data & gset ----
sampleNames(gset) <- str_sub(sampleNames(gset), 1, 10)
pD <- pD[sampleNames(gset),]
pData(gset) <- pD
gset
# GeneFeatureSet (storageMode: lockedEnvironment)
# assayData: 1178100 features, 70 samples 
# element names: exprs 
# protocolData
# rowNames: GSM3302215 GSM3302216 ... GSM3302455 (70 total)
# varLabels: exprs dates
# varMetadata: labelDescription channel
# phenoData
# rowNames: GSM3302215 GSM3302216 ... GSM3302455 (70 total)
# varLabels: title subjectid ... health (8 total)
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


#### Normalization ####
normData <- rma(gset, target = "core")
normData  ## rma returns "ExpressionSet"
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 33297 features, 70 samples 
# element names: exprs 
# protocolData
# rowNames: GSM3302215 GSM3302216 ... GSM3302455 (70 total)
# varLabels: exprs dates
# varMetadata: labelDescription channel
# phenoData
# rowNames: GSM3302215 GSM3302216 ... GSM3302455 (70 total)
# varLabels: title subjectid ... health (8 total)
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
#         GSM3302215 GSM3302216 GSM3302218
# 7892501   2.911407   4.108206   2.975292
# 7892502   3.693875   3.756861   3.724263
# 7892503   1.830768   1.809343   1.944873


#### Quality control ####
#. PCA analysis ----
myexp <- Biobase::exprs(normData)
PCA <- prcomp(t(myexp), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pDnorm <- pData(normData)
head(pDnorm, 3)
dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     individual = pDnorm$subjectid,
                     gender = pDnorm$gender,
                     timepoint = pDnorm$timepoint,
                     health = pDnorm$health
                     )

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(. ~ gender) +
    ggtitle("PCA_Norm timepoint & gender") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5))+
    coord_fixed(ratio = sd_ratio) +
    theme(legend.position = "none")



#### Filtering based on intensity ####
# filter out lowly expressed genes
# “soft” intensity based filtering here, since this is recommended by the limma
#. Histogram ----
medians <- rowMedians(Biobase::exprs(normData))
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

#set a manual threshold to exclude genes, but not too many
# man_threshold <- round(mean(medians))
man_threshold <- 4
abline(v = man_threshold, col = "coral4", lwd = 2)

#. Exclude genes that below the threshold ----
# Transcripts that do not have intensities larger than the threshold in
# at least as many arrays as the smallest experimental group are excluded.
# get a list:
no_of_samples <- table(paste0(pDnorm$timepoint, "_", pDnorm$gender))
no_of_samples

samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(Biobase::exprs(normData), 1,
                           function(x){
                               sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)
# FALSE  TRUE 
# 11508 21789

manfiltered <- subset(normData, idx_man_threshold)
dim(exprs(normData))  # 33297    70
dim(exprs(manfiltered))  # 21789    70



#### Annotation of the transcript clusters ####
# columns(hugene11sttranscriptcluster.db)
# keytypes(hugene11sttranscriptcluster.db)
# head(keys(hugene11sttranscriptcluster.db, keytype="GOALL"))
# ls("package:hugene11sttranscriptcluster.db")
anno <- AnnotationDbi::select(hugene11sttranscriptcluster.db,
                              keys = (featureNames(manfiltered)),
                              columns = c("SYMBOL", "GENENAME"),
                              keytype = "PROBEID")
table(is.na(anno))
# FALSE  TRUE 
# 69926 10750 


anno <- subset(anno, !is.na(SYMBOL))
sum(is.na(anno$SYMBOL)) # 0
nosymbols <- !(featureNames(manfiltered) %in% anno$PROBEID)
table(nosymbols)
# FALSE  TRUE 
# 16414  5375

withsymbols <- subset(manfiltered, !nosymbols)
dim(withsymbols)
# Features  Samples 
# 16414       70


#. Removing multiple mappings ----
# group by PROBEID
anno_grouped <- group_by(anno, PROBEID)
head(anno_grouped)

# summarize
anno_summarized <-
    dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized)

# filter to extract cluster multiple gene that map to multiple gene symbols
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
head(anno_filtered)

probe_stats <- anno_filtered
nrow(probe_stats) # 1643: clusters that map to multiple gene symbols → remove

# remove above IDs(probe_stats)
ids_to_exlude <- (featureNames(withsymbols) %in% probe_stats$PROBEID)
table(ids_to_exlude)
# FALSE  TRUE 
# 14771  1643

final <- subset(withsymbols, !ids_to_exlude)
validObject(final)
dim(final)  # 14771  70 


# also exclude them from the feature data anno
head(anno, 3)

# generate a column PROBEID in fData(final) and
# assign the row names of fData(final) to it:
fData(final)$PROBEID <- rownames(fData(final))

# left-join fData(final) with anno
fData(final) <- left_join(fData(final), anno)

rownames(fData(final)) <- fData(final)$PROBEID
validObject(final)




#### Differential Expression nalysis ####
head(pData(final), 3)
individual <- Biobase::pData(final)$subjectid
timepoint <- Biobase::pData(final)$timepoint
health <- Biobase::pData(final)$health
gender <- Biobase::pData(final)$gender

# normData (without supplement, health)
male <- individual[gender == "M"]
t_male <- timepoint[gender == "M"]
design_male <- model.matrix(~ 0 + t_male + male)
colnames(design_male)[1:2] <- c("pre", "post")
rownames(design_male) <- male

female <- individual[gender == "F"]
t_female <- timepoint[gender == "F"]
design_female <- model.matrix(~ 0 + t_female + female)
colnames(design_female)[1:2] <- c("pre", "post")
rownames(design_female) <- female

# inspect the design matrices
design_male[1:5, 1:5]
design_female[1:5, 1:5]


# Contrasts and hypotheses tests: all genes ----
contrast_matrix_male <- makeContrasts(post-pre, levels = design_male)
fit_male <- eBayes(
    contrasts.fit(
        lmFit(
            final[,gender == "M"],
            design = design_male),
        contrast_matrix_male))

contrast_matrix_female <- makeContrasts(post-pre, levels = design_female)
fit_female <- eBayes(
    contrasts.fit(
        lmFit(
            final[,gender == "F"],
            design = design_female),
        contrast_matrix_female))


#### Extracting results: topTable() ####
library(RColorBrewer)
table_male <- topTable(fit_male, number = Inf)
head(table_male)
hist(table_male$P.Value, col = brewer.pal(3, name = "Set2")[1], xlab = "male p-values")
hist(table_male$adj.P.Val, col = brewer.pal(3, name = "Set2")[2], xlab = "male adj.p-values")

table_female <- topTable(fit_female, number = Inf)
head(table_female)
hist(table_female$P.Value, col = brewer.pal(3, name = "Set2")[1], xlab = "female p-values")
hist(table_female$adj.P.Val, col = brewer.pal(3, name = "Set2")[2], xlab = "female p-values")


#. Multiple testing FDR, and comparison with results from the original paper
nrow(subset(table_male, P.Value < 0.05))  # 2120 (cf. before annotation 2698)
tail(subset(table_male, P.Value < 0.05))
nrow(subset(table_male, adj.P.Val < 0.01))  # 33
tail(subset(table_male, adj.P.Val < 0.01))

nrow(subset(table_female, P.Value < 0.05))  # 1652 (cf. before annotation 2254)
tail(subset(table_female, P.Value < 0.05))
nrow(subset(table_female, adj.P.Val < 0.01))  # 0
tail(subset(table_female, adj.P.Val < 0.01))


# Exploring Results
# resSig_male <- subset(table_male, P.Value < 0.05)
# head(resSig_male[order(resSig_male$P.Value), ])
# resSig_female <- subset(table_female, P.Value < 0.05)
# head(resSig_female[order(resSig_female$P.Value), ])



#### Write out the result ####

write.csv(table_male, file = "res_GSE1295M.csv")
write.csv(table_female, file = "res_GSE1295F.csv")

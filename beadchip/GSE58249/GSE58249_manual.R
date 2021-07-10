#### GSE58249 ####
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58249
# pubmed: https://pubmed.ncbi.nlm.nih.gov/25138607/
# ncbi full txt: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4200377/
# Illumina HumanHT-12 V4.0 expression beadchip
# BA Workflow: http://www.bioconductor.org/packages/release/data/experiment/vignettes/BeadArrayUseCases/inst/doc/BeadArrayUseCases.pdf
# getGEO(): normalized data
# resistance (n = 9) / aerobic (n = 8)
# 13 female, 5 male

#### Loading packages ####
library(Biobase)
library(oligoClasses)
library(oligo)
library(arrayQualityMetrics)
# args(arrayQualityMetrics)
# help(arrayQualityMetrics)
library(limma)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(stringr)
# library(matrixStats)
# library(genefilter)
library(openxlsx)
library(devtools)


#### Loading a data ####
library(GEOquery)
geo <- getGEO(GEO = "GSE58249")
length(geo)  # 1
gse <- geo[[1]]
gse
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 47323 features, 34 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM1404698 GSM1404699 ... GSM1404731 (34 total)
# varLabels: title geo_accession ... tissue type:ch1 (40 total)
# varMetadata: labelDescription
# featureData
# featureNames: ILMN_1343291 ILMN_1343295 ... ILMN_3311190 (47323 total)
# fvarLabels: ID Species ... GB_ACC (30 total)
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 25138607 
# Annotation: GPL10558 

dim(gse)
# Features  Samples 
# 47323       34 
exprs(gse)[1:5, 1:5]
#              GSM1404698 GSM1404699 GSM1404700 GSM1404701 GSM1404702
# ILMN_1343291         NA         NA   10.62050   10.57280   10.47720
# ILMN_1343295         NA         NA   11.79760   12.10130   11.98290
# ILMN_1651199         NA         NA    6.25077    6.12807    6.14016
# ILMN_1651209         NA         NA    6.19532    6.24520    6.25682
# ILMN_1651210         NA         NA    6.20404    6.34881    6.16962


# checking if the loaded data is already normalized ----
head(pData(gse)$data_processing, 1)
# [1] "Data was processed in Genome Studio, quantile normalized and
#      removal of probes not detacted in at least 10% of the samples"

oligo::boxplot(log2(exprs(gse)), target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               las = 2,
               outline = FALSE,
               cex.axis=0.7)


#### Data cleaning ####
#. Subselect the columns of interest ----
pD.all <- pData(gse)
head(pD.all, 3)
tail(pD.all, 3)
pD <- pD.all[, c("title",  # individual
                 "male or female:ch1",  # gender
                 "exercise type (aerobic or resistance):ch1",  # ex type
                 "pre exercise/ post exercise:ch1",  # timepoint
                 "description",
                 "description.1"
                 )]

#. Renaming ----
names(pD)[1:4] <- c("individual", "gender", "exercise.type", "timepoint")
head(pD, 3)
#               individual gender exercise.type timepoint           description                                                   description.1
# GSM1404698  Pre-2-PRT-M      M    RESISTANCE       PRE  Pre-2 (4 replicates) This sample not included in final analysis (no normalized data)
# GSM1404699 Post-2-PRT-M      M    RESISTANCE      POST Post-2 (3 replicates) This sample not included in final analysis (no normalized data)
# GSM1404700  Pre-4-PRT-F      F    RESISTANCE       PRE

pD$individual <- str_sub(pD$individual, -8, -7)
pD$individual <- str_replace_all(pD$individual, "-", "")
pD$individual <- paste0("S", pD$individual)

head(pD)
tail(pD)
dim(pD)
levels(factor(pD$exercise.type))  # "AEROBIC"    "RESISTANCE"
levels(factor(pD$timepoint))  # "POST" "PRE" 
pD$timepoint <- relevel(factor(tolower(pD$timepoint)), ref = "pre")
levels(factor(pD$timepoint))  # "pre"  "post"
sum(table(unique(pD$individual)))  # 17

# note: original study didn't include below subjects.
# this analysis also doesn't included them based on the ration of NAs later 
# (ref: gse <- gse[!bad.gene,!bad.sample])
pD %>% filter(grepl("not included", pD$description) | grepl("not included", pD$description.1))
# individual gender exercise.type timepoint                                                     description
# GSM1404698   Pre-2-PRT-M      M    RESISTANCE       PRE                                            Pre-2 (4 replicates)
# GSM1404699  Post-2-PRT-M      M    RESISTANCE      POST                                           Post-2 (3 replicates)
# GSM1404720  Pre-19-AER-F      F       AEROBIC       PRE This sample not included in final analysis (no normalized data)
# GSM1404721 Post-19-AER-F      F       AEROBIC      POST This sample not included in final analysis (no normalized data)
# 
# description.1
# GSM1404698 This sample not included in final analysis (no normalized data)
# GSM1404699 This sample not included in final analysis (no normalized data)
# GSM1404720                                                                
# GSM1404721



#### Quality control of the raw data ####
# filtering out probes with poor annotation ----
# ref: http://www.bioconductor.org/packages/release/data/experiment/vignettes/BeadArrayUseCases/inst/doc/BeadArrayUseCases.pdf
library(illuminaHumanv4.db)
ids <- featureNames(gse)
head(ids)
# [1] "ILMN_1343291" "ILMN_1343295" "ILMN_1651199" "ILMN_1651209" "ILMN_1651210" "ILMN_1651221"

probeQuality <- unlist(mget(ids, illuminaHumanv4PROBEQUALITY, ifnotfound = NA))
head(probeQuality)
# ILMN_1343291 ILMN_1343295 ILMN_1651199 ILMN_1651209 ILMN_1651210 ILMN_1651221 
# "Perfect"    "Perfect"        "Bad"    "Perfect"        "Bad"        "Bad"

table(probeQuality)
# probeQuality
#   Bad   Good     Good***    Good****    No match     Perfect  Perfect*** Perfect**** 
# 12379   1159         137         507         468       27367        2911        2395

remove <- probeQuality == "No match" | probeQuality == "Bad" | is.na(probeQuality)
gse <- gse[!remove,]

dim(gse)
# Features  Samples 
# 34476       34

# Take a look
exprs(gse)[1:5, 1:5]
exprs(gse)[,1:2]

# if is.na() === TRUE ----
bad.sample <- colMeans(is.na(exprs(gse))) > 0.8
bad.gene <- rowMeans(is.na(exprs(gse))) > 0.5
gse <- gse[!bad.gene,!bad.sample]
dim(gse)
# Features  Samples 
# 34476       30 

table(bad.sample)
keep <- rownames(pD) %in% sampleNames(gse)
pD <- pD[keep,]
head(pD, 3)
#             individual gender exercise.type timepoint description description.1
# GSM1404700         S4      F    RESISTANCE       pre                          
# GSM1404701         S4      F    RESISTANCE      post                          
# GSM1404702         S6      F       AEROBIC       pre
dim(pD)  # 30  6


#### Normalization (skipped) ####
# this data is already normalized.
head(pData(gse)$data_processing, 1)
# [1] "Data was processed in Genome Studio, quantile normalized and
#      removal of probes not detacted in at least 10% of the samples"


# PCA Analysis ----
# log2
exp_raw <- log2(exprs(gse))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

# PCA plot
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     individual = pD$individual,
                     gender = pD$gender,
                     timepoint = pD$timepoint,
                     exercise.type = pD$exercise.type
                     )

# gender × exercise.type
ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(gender ~ exercise.type)



#### Differencial Expression Analysis ####
head(pD, 3)
individual <- pD$individual
gender <- pD$gender
exercise.type <- pD$exercise.type
timepoint <- pD$timepoint

#### lmFit() ####
# male ----
#. male * resistance ----
i_MR <- individual[gender == "M" & exercise.type == "RESISTANCE"]
t_MR <- timepoint[gender == "M" & exercise.type == "RESISTANCE"]
design_MR <- model.matrix(~ 0 + t_MR + i_MR)
colnames(design_MR)[1:2] <- c("pre", "post")
rownames(design_MR) <- i_MR
gse_male_resistance <- gse[, gender == "M" & exercise.type == "RESISTANCE"]
fit_MR <- lmFit(gse_male_resistance, design_MR)
contrast_MR <- makeContrasts(post-pre, levels = design_MR)
contr.fit.MR <- eBayes(contrasts.fit(fit_MR, contrast_MR))

#. male * aerobic ----
i_MA <- individual[gender == "M" & exercise.type == "AEROBIC"]
t_MA <- timepoint[gender == "M" & exercise.type == "AEROBIC"]
design_MA <- model.matrix(~ 0 + t_MA + i_MA)
colnames(design_MA)[1:2] <- c("pre", "post")
rownames(design_MA) <- i_MA
gse_male_aerobic <- gse[, gender == "M" & exercise.type == "AEROBIC"]
fit_MA <- lmFit(gse_male_aerobic, design_MA)
contrast_MA <- makeContrasts(post-pre, levels = design_MA)
contr.fit.MA <- eBayes(contrasts.fit(fit_MA, contrast_MA))

# female ----
#. female * resistance ----
i_FR <- individual[gender == "F" & exercise.type == "RESISTANCE"]
t_FR <- timepoint[gender == "F" & exercise.type == "RESISTANCE"]
design_FR <- model.matrix(~ 0 + t_FR + i_FR)
colnames(design_FR)[1:2] <- c("pre", "post")
rownames(design_FR) <- i_FR
gse_female_resistance <- gse[, gender == "F" & exercise.type == "RESISTANCE"]
fit_FR <- lmFit(gse_female_resistance, design_FR)
contrast_FR <- makeContrasts(post-pre, levels = design_FR)
contr.fit.FR <- eBayes(contrasts.fit(fit_FR, contrast_FR))

#. female * aerobic ----
i_FA <- individual[gender == "F" & exercise.type == "AEROBIC"]
t_FA <- timepoint[gender == "F" & exercise.type == "AEROBIC"]
design_FA <- model.matrix(~ 0 + t_FA + i_FA)
colnames(design_FA)[1:2] <- c("pre", "post")
rownames(design_FA) <- i_FA
gse_female_aerobic <- gse[, gender == "F" & exercise.type == "AEROBIC"]
fit_FA <- lmFit(gse_female_aerobic, design_FA)
contrast_FA <- makeContrasts(post-pre, levels = design_FA)
contr.fit.FA <- eBayes(contrasts.fit(fit_FA, contrast_FA))


#### annotation ####
# male * resistance ----
anno_ids_mr <- rownames(gse_male_resistance)
entrezid <- mget(anno_ids_mr, illuminaHumanv4ENTREZID, ifnotfound = NA)
symbols <- mget(anno_ids_mr, illuminaHumanv4SYMBOL, ifnotfound = NA)
genename <- mget(anno_ids_mr, illuminaHumanv4GENENAME, ifnotfound = NA)
annotationDf <- data.frame(Symbol = as.character(symbols),
                           Name = as.character(genename),
                           EntrezID = as.numeric(entrezid))
contr.fit.MR$genes <- annotationDf


# male * aerobic ----
anno_ids_ma <- rownames(gse_male_aerobic)
entrezid <- mget(anno_ids_ma, illuminaHumanv4ENTREZID, ifnotfound = NA)
symbols <- mget(anno_ids_ma, illuminaHumanv4SYMBOL, ifnotfound = NA)
genename <- mget(anno_ids_ma, illuminaHumanv4GENENAME, ifnotfound = NA)
annotationDf <- data.frame(Symbol = as.character(symbols),
                           Name = as.character(genename),
                           EntrezID = as.numeric(entrezid))
contr.fit.MA$genes <- annotationDf


# female * resistance ----
anno_ids_fr <- rownames(gse_female_resistance)
entrezid <- mget(anno_ids_fr, illuminaHumanv4ENTREZID, ifnotfound = NA)
symbols <- mget(anno_ids_fr, illuminaHumanv4SYMBOL, ifnotfound = NA)
genename <- mget(anno_ids_fr, illuminaHumanv4GENENAME, ifnotfound = NA)
annotationDf <- data.frame(Symbol = as.character(symbols),
                           Name = as.character(genename),
                           EntrezID = as.numeric(entrezid))
contr.fit.FR$genes <- annotationDf


# female * aerobic ----
anno_ids_fa <- rownames(gse_female_aerobic)
entrezid <- mget(anno_ids_fa, illuminaHumanv4ENTREZID, ifnotfound = NA)
symbols <- mget(anno_ids_fa, illuminaHumanv4SYMBOL, ifnotfound = NA)
genename <- mget(anno_ids_fa, illuminaHumanv4GENENAME, ifnotfound = NA)
annotationDf <- data.frame(Symbol = as.character(symbols),
                           Name = as.character(genename),
                           EntrezID = as.numeric(entrezid))
contr.fit.FA$genes <- annotationDf


#### Results ####
# toptable() ----
table_mr <- topTable(contr.fit.MR, coef = 1, number = Inf)
table_ma <- topTable(contr.fit.MA, coef = 1, number = Inf)
table_fr <- topTable(contr.fit.FR, coef = 1, number = Inf)
table_fa <- topTable(contr.fit.FA, coef = 1, number = Inf)


# writing out the results ----
write.csv(table_mr, file = "res_GSE58249MR.csv")
write.csv(table_ma, file = "res_GSE58249MA.csv")
write.csv(table_fr, file = "res_GSE58249FR.csv")
write.csv(table_fa, file = "res_GSE58249FA.csv")



# histgram ----
library(RColorBrewer)
par(mfrow = c(2,2))

#. male ----
hist(table_mr$P.Value, col = brewer.pal(3, name = "Set2")[1], 
     xlab = "mr Pval", main = "mr PVal")
hist(table_mr$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
     xlab = "mr adjPval", main = "mr adjPval")
hist(table_ma$P.Value, col = brewer.pal(3, name = "Set2")[1],
     xlab = "ma Pval", main = "ma Pval")
hist(table_ma$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
     xlab = "ma adjPval", main = "ma adjPval")

#. female ----
hist(table_fr$P.Value, col = brewer.pal(3, name = "Set2")[1], 
     xlab = "fr Pval", main = "fr PVal")
hist(table_fr$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
     xlab = "fr adjPval", main = "fr adjPval")
hist(table_fa$P.Value, col = brewer.pal(3, name = "Set2")[1],
     xlab = "fa Pval", main = "fa Pval")
hist(table_fa$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
     xlab = "fa adjPval", main = "fa adjPval")


# tables ----
#. male * resistance ---- 
nrow(subset(table_mr, P.Value < 0.05))  # 961
tail(subset(table_mr, P.Value < 0.05))
nrow(subset(table_mr, adj.P.Val < 0.05))  # 0
tail(subset(table_mr, adj.P.Val < 0.05))

#. male * aerobic ----
nrow(subset(table_ma, P.Value < 0.05))  # 2225
tail(subset(table_ma, P.Value < 0.05))
nrow(subset(table_ma, adj.P.Val < 0.05))  # 59
tail(subset(table_ma, adj.P.Val < 0.05))

#. female * resistance ----
nrow(subset(table_fr, P.Value < 0.05))  # 2683 (cf. before annotation 2254)
tail(subset(table_fr, P.Value < 0.05))
nrow(subset(table_fr, adj.P.Val < 0.05))  # 53
tail(subset(table_fr, adj.P.Val < 0.05))

#. female * aerobic ----
nrow(subset(table_fa, P.Value < 0.05))  # 2295
tail(subset(table_fa, P.Value < 0.05))
nrow(subset(table_fa, adj.P.Val < 0.05))  # 56
tail(subset(table_fa, adj.P.Val < 0.01))

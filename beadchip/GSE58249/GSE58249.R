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
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
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


#### Data cleaning ####
# featuredata ----
head(fData(gse), 3)
fData(gse) <- fData(gse)[,c("Probe_Id",
                            "Accession",
                            "Source",
                            "Symbol",
                            "Entrez_Gene_ID",
                            "Synonyms"
                            )]
names(fData(gse)) <- c("PROBEID", "ACCNUM", "SOURCE", "SYMBOL","ENTREZID", "SYNONYMS")
head(fData(gse), 3)


# phenodata ----
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
pD$group <- paste0(pD$gender , str_sub(pD$exercise.type, 1,1))

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
#             individual gender exercise.type timepoint description description.1 group
# GSM1404700         S4      F    RESISTANCE       pre                              FR
# GSM1404701         S4      F    RESISTANCE      post                              FR
# GSM1404702         S6      F       AEROBIC       pre                              FA
dim(pD)  # 30  7


#### Normalization (skipped) ####
# this data is already normalized.
head(pData(gse)$data_processing, 1)
# [1] "Data was processed in Genome Studio, quantile normalized and
#      removal of probes not detacted in at least 10% of the samples"

oligo::boxplot(log2(exprs(gse)), target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               las = 2,
               outline = FALSE,
               cex.axis=0.7)


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
                     exercise.type = pD$exercise.type,
                     group = pD$group
)

# gender × exercise.type
ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(gender ~ exercise.type)


####. More elaborate quality control plots (skipped) ####
# codes are moved to memo.Rmd



#### Differencial Expression Analysis ####
head(pD, 3)
individual <- pD$individual
timepoint <- pD$timepoint
group <- pD$group
# gender <- pD$gender
# exercise.type <- pD$exercise.type
# timepoint <- pD$timepoint

#### lmFit() ####
for (i in unique(group)) {
    print(i)
    subjects <- individual[group == i]
    pre_post <- timepoint[group == i]
    data <- gse[, group == i]
    
    design <- model.matrix(~ 0 + pre_post + subjects)
        colnames(design)[1:2] <- c("pre", "post")
        rownames(design) <- subjects
    fit <- lmFit(data, design)
    contrast <- makeContrasts(post-pre, levels = design)
    contr.fit <- eBayes(contrasts.fit(fit, contrast))
    
    #### Results ####
    table <- topTable(contr.fit, coef = 1, number = Inf)
    write.csv(table, file = paste0("res_GSE58249", i, ".csv"))
}


#### Visualizing results ####
library(RColorBrewer)
result_files <- list.files(pattern = ".csv")
par(mfrow = c(2,2))

for (i in 1:length(result_files)) {
    file <- result_files[i]
    results <- read.csv(file)
    group <- str_sub(file, -6, -5) 
    
    ## histgram ##
    hist(results$P.Value, col = brewer.pal(3, name = "Set2")[1], 
         main = paste(group, "Pval"), xlab  = NULL)
    hist(results$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
         main = paste(group, "adj.Pval"), xlab = NULL)

    ## some numbers ##
    cat(group, "\n")
    cat("p < 0.05:", nrow(subset(results, P.Value < 0.05)), "\n")
    cat("adj.P < 0.05:", nrow(subset(results, adj.P.Val < 0.05)),"\n")
    cat("adj.P < 0.01:", nrow(subset(results, adj.P.Val < 0.01)),"\n\n")
    }

# to explore more (ex)
# tail(subset(results, adj.P.Val < 0.05))




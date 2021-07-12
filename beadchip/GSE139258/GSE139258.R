#### GSE139258 ####
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
geo <- getGEO(GEO = "GSE139258")
length(geo)  # 1
gse <- geo[[1]]
gse
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 47231 features, 24 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM4134936 GSM4134937 ... GSM4134959 (24 total)
# varLabels: title geo_accession ... tissue:ch1 (45 total)
# varMetadata: labelDescription
# featureData
# featureNames: ILMN_1343291 ILMN_1343295 ... ILMN_3311190 (47231 total)
# fvarLabels: ID Species ... GB_ACC (30 total)
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 32562350 
# Annotation: GPL10558


dim(gse)
# Features  Samples 
# 47231       24 
exprs(gse)[1:5, 1:5]  # non log-transformed data
#              GSM4134936 GSM4134937 GSM4134938 GSM4134939 GSM4134940
# ILMN_1343291 22362.8200 23690.2000 19798.9500 40214.4500 25279.4700
# ILMN_1343295 58602.3800 65495.7300 60236.0700 57477.4600 57218.5900
# ILMN_1651199   130.5734   128.6970   127.2887   124.2555   127.0545
# ILMN_1651209   147.5510   129.2871   175.1209   172.1879   166.6653
# ILMN_1651210   121.1727   141.3416   127.7274   133.1908   140.4606



#### Data cleaning ####
# featuredata ----
head(fData(gse), 3)
names(fData(gse))
fData(gse) <- fData(gse)[,c("Probe_Id",
                            "GB_ACC",
                            "Symbol",
                            "Entrez_Gene_ID"
)]
names(fData(gse)) <- c("PROBEID", "ACCNUM", "SYMBOL","ENTREZID")
head(fData(gse), 3)
#                   PROBEID      ACCNUM    SYMBOL ENTREZID
# ILMN_1343291 ILMN_1343291 NM_001402.5    EEF1A1     1915
# ILMN_1343295 ILMN_1343295 NM_002046.3     GAPDH     2597
# ILMN_1651199 ILMN_1651199 XM_944551.1 LOC643334   643334

### memo: GENENAME may be needed in the later analysis


# phenodata ----
#. Subselect the columns of interest ----
pD.all <- pData(gse)
head(pD.all, 3)
tail(pD.all, 3)
pD <- pD.all[, c("title", # individual
                 "Sex:ch1", # gender
                 "groups:ch1", # group type
                 "group:ch1", # group type
                 "time:ch1"
                 )]
names(pD) <- c("individual", "gender", "groupShort", "groupLong", "timepoint")
# cf. LOS == Lean/Overweight Sedentary

pD
# (hide) ----
#                                              dividual  gender groupShort     group timepoint
# GSM4134936     LOS_PRE Male LEAN 4 days 10001322038_A    Male        LOS      LEAN      PRE
# GSM4134937   LOS_POST Male LEAN 23 days 10001322038_B    Male        LOS      LEAN     POST
# GSM4134938   LOS_PRE Female LEAN 4 days 10001322038_C  Female        LOS      LEAN      PRE
# GSM4134939 LOS_POST Female LEAN 23 days 10001322038_D  Female        LOS      LEAN     POST
# GSM4134940   LOS_PRE Female LEAN 4 days 10001322038_E  Female        LOS      LEAN      PRE
# GSM4134941 LOS_POST Female LEAN 23 days 10001322038_F  Female        LOS      LEAN     POST
# GSM4134942   LOS_PRE Female LEAN 4 days 10001322038_G  Female        LOS      LEAN      PRE
# GSM4134943 LOS_POST Female LEAN 23 days 10001322038_H  Female        LOS      LEAN     POST
# GSM4134944   LOS_PRE Female LEAN 4 days 10001322038_I  Female        LOS      LEAN      PRE
# GSM4134945 LOS_POST Female LEAN 23 days 10001322038_J  Female        LOS      LEAN     POST
# GSM4134946     LOS_PRE Male LEAN 4 days 10001322038_K    Male        LOS      LEAN      PRE
# GSM4134947   LOS_POST Male LEAN 23 days 10001322038_L    Male        LOS      LEAN     POST
# GSM4134948   LA_PRE Male ATHLETE 4 days 10001322039_A    Male         LA   ATHLETE      PRE
# GSM4134949   LA_PRE Male ATHLETE 4 days 10001322039_B    Male         LA   ATHLETE      PRE
# GSM4134950   LA_PRE Male ATHLETE 4 days 10001322039_C    Male         LA   ATHLETE      PRE
# GSM4134951   LA_PRE Male ATHLETE 4 days 10001322039_D    Male         LA   ATHLETE      PRE
# GSM4134952   LA_PRE Male ATHLETE 4 days 10001322039_E    Male         LA   ATHLETE      PRE
# GSM4134953   LA_PRE Male ATHLETE 4 days 10001322039_F    Male         LA   ATHLETE      PRE
# GSM4134954    OBS_PRE Male obese 4 days 10001322039_G    Male        OBS     obese      PRE
# GSM4134955  OBS_PRE Female obese 4 days 10001322039_H  Female        OBS     obese      PRE
# GSM4134956  OBS_PRE Female obese 4 days 10001322039_I  Female        OBS     obese      PRE
# GSM4134957  OBS_PRE Female obese 4 days 10001322039_J  Female        OBS     obese      PRE
# GSM4134958  OBS_PRE Female obese 4 days 10001322039_K  Female        OBS     obese      PRE
# GSM4134959  OBS_PRE Female obese 4 days 10001322039_L  Female        OBS     obese      PRE

# excluding LA & OBS (no data for post exercise) ----
keep <- pD$groupShort == "LOS"
gse <- gse[,keep]
pD <- pD[keep,]


pD$gender <- str_sub(pD$gender, 1, 1)
levels(factor(pD$timepoint))  # "POST" "PRE" 
pD$timepoint <- relevel(factor(tolower(pD$timepoint)), ref = "pre")
levels(factor(pD$timepoint))  # "pre"  "post"


# skipped (subjects cannot be identified)
# pData(dat)$individual <- factor(
#     c("A", "A", "B", "B", "C", "C", "D", "D", "E", "E", "F", "F"))
# male: A,F
# female: B,C,D,E

# pD$id <- paste(pD$individual, pD$timepoint, sep = ".")



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
# 34476       12

exprs(gse)[1:5, 1:5]


# if is.na() === TRUE ----
bad.sample <- colMeans(is.na(exprs(gse))) > 0.8
bad.gene <- rowMeans(is.na(exprs(gse))) > 0.5
gse <- gse[!bad.gene,!bad.sample]
dim(gse)
# Features  Samples 
# 34476       12 

table(bad.sample)
# bad.sample
# FALSE 
# 12



#### Normalization (skipped) ####
# this data is already normalized.
head(pData(gse)$data_processing, 1)
# [1] "The data were normalised using cubic spline normalization with IlluminaGUI in R"

oligo::boxplot(log2(exprs(gse)), target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               las = 2,
               outline = FALSE,
               cex.axis=0.7)



#### PCA Analysis ####
# log2
exp_raw <- log2(exprs(gse))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

# head(pD, 3)
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     # individual = pD$individual,  # cannot identify
                     gender = pD$gender,
                     timepoint = pD$timepoint
                     )

# gender × exercise.type
ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(. ~ gender)


####. More elaborate quality control plots (skipped) ####
# codes are moved to memo.Rmd



#### Differencial Expression Analysis ####
head(pD, 3)
# individual <- pD$individual
timepoint <- pD$timepoint
gender <- pD$gender

#### lmFit() ####
exprs(gse) <- log2(exprs(gse))
for (i in unique(gender)) {
    print(i)
    # subjects <- individual[group == i]
    pre_post <- timepoint[gender == i]
    data <- gse[, gender == i]
    
    design <- model.matrix(~ 0 + pre_post)  # removed: + subjects
    colnames(design)[1:2] <- c("pre", "post")
    # rownames(design) <- subjects
    fit <- lmFit(data, design)
    contrast <- makeContrasts(post - pre, levels = design)
    contr.fit <- eBayes(contrasts.fit(fit, contrast))
    
    #### Results ####
    table <- topTable(contr.fit, coef = 1, number = Inf)
    write.csv(table, file = paste0("res_GSE139258re", i, ".csv"))
    print(head(table))
}


#### Visualizing results ####
library(RColorBrewer)
result_files <- list.files(pattern = ".csv")
result_files
par(mfrow = c(2,2))

for (i in 1:length(result_files)) {
    file <- result_files[i]
    results <- read.csv(file)
    group <- str_sub(file, -5, -5)
    
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




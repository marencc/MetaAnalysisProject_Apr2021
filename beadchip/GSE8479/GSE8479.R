# GSE8479
# Resistance Exercise Reverses Aging in Human Skeletal Muscle
# n = 14 (m6 f8)
# older
# GPL2700
# Sentrix HumanRef-8 Expression BeadChip

#### Loading packages ####
library(Biobase)
library(oligoClasses)
library(oligo)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(genefilter)
library(devtools)


#### loading the data ####
library(GEOquery)
if (!file.exists("data")) {
    dir.create("data")
    }

geo <- getGEO("GSE8479", destdir = "data")
length(geo)  # 1
gse <- geo[[1]]

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
exprs(gse)[1:5 , 1:2] # non log-transformed data

head(fData(gse))  # GB_ACC will be used for annotation
#                          ID   SequenceSource    GB_ACC SPOT_ID
# GI_10047089-S GI_10047089-S         RefSeq NM_014332.1      NA
# GI_10047091-S GI_10047091-S         RefSeq NM_013259.1      NA
# GI_10047093-S GI_10047093-S         RefSeq NM_016299.1      NA
# GI_10047099-S GI_10047099-S         RefSeq NM_016303.1      NA
# GI_10047103-S GI_10047103-S         RefSeq NM_016305.1      NA
# GI_10047105-S GI_10047105-S         RefSeq NM_016352.1      NA



#### Annotation ####
library(org.Hs.eg.db)  # database stored at: dbfile(org.Hs.eg.db)

# ref: https://support.bioconductor.org/p/119296/
GB_ACC <- sapply(strsplit(fData(gse)[,3], "\\."), "[", 1)
head(GB_ACC)
# [1] "NM_014332" "NM_013259" "NM_016299" "NM_016303" "NM_016305" "NM_016352"

sum(is.na(GB_ACC))  # 4
GB_ACC[is.na(GB_ACC)] <- "NA"
# columns(org.Hs.eg.db)
anno <- lapply(c("ENTREZID","SYMBOL","GENENAME"), 
               function(x) mapIds(org.Hs.eg.db, GB_ACC, x, "ACCNUM"))

annodf <- data.frame(PROBEID = fData(gse)[,1], 
                     ACCNUM = GB_ACC,
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
ids <- featureNames(gse)
head(featureNames(gse))
# [1] "GI_10047089-S" "GI_10047091-S" "GI_10047093-S" "GI_10047099-S" "GI_10047103-S"  "GI_10047105-S"

library(illuminaHumanv1.db)  # memo: CHIPNAME: HumanWG6v1
probe_quality <- unlist(mget(ids, illuminaHumanv1PROBEQUALITY, ifnotfound = NA))
length(ids)  # [1] 24354
length(probe_quality)  # [1] 24361

# checking the cause of the length difference
ids_df <- data.frame(ids = ids, original_ids = ids)
pq_df <- data.frame(ids = names(probe_quality), pq_ids = names(probe_quality))
pq_ids_df <-full_join(ids_df, pq_df)
pq_ids_df[is.na(pq_ids_df$pq_ids),]
#                 ids  original_ids pq_ids
# 1492  GI_14141192-S GI_14141192-S   <NA>
# 7608  GI_25453469-S GI_25453469-S   <NA>
# 12538 GI_34304116-S GI_34304116-S   <NA>
# 21010  GI_4507728-S  GI_4507728-S   <NA>
# 21013  GI_4507744-S  GI_4507744-S   <NA>
# 22035  GI_5016088-S  GI_5016088-S   <NA>
# 23365  GI_7669491-S  GI_7669491-S   <NA>

# nrow(pq_ids_df[is.na(pq_ids_df$pq_ids),])
# [1] 7   

#### memo:
# the length difference probably comes from above 7 names
# length(ids)  # [1] 24354
# length(probe_quality)  # [1] 24361 

# somehow the above 7 names are not counted as NA for "probe_quality" when check the following
probe_quality[is.na(probe_quality)]
# lysA pheA thrB trpF 
# NA   NA   NA   NA

table(is.na(probe_quality))
# FALSE  TRUE 
# 24357     4

# removing all 7 names from original ids
remove <- ids %in% pq_ids_df[is.na(pq_ids_df$pq_ids),]$ids
table(remove)
# FALSE  TRUE 
# 24347     7
ids2 <- ids[!remove]
probe_quality2 <- unlist(mget(ids2, illuminaHumanv1PROBEQUALITY, ifnotfound = NA))
length(ids2)  # [1] 24347
length(probe_quality2)  # [1] 24347


# skipped because this length difference may be causing error for "gse[!rem2,]" below
table(probe_quality2)
# probe_quality2
# Bad        Good     Good***    Good****    No match     Perfect  Perfect*** Perfect**** 
#     4040         767          24         159          88       17536         657        1072

# cf.
# probe_quality
#       Bad        Good     Good***    Good****    No match     Perfect  Perfect***  Perfect**** 
#     4040         767          24         160          88       17548         657          1073


probes_remove <- probe_quality2 == "No match" | probe_quality2 == "Bad"  | is.na(probe_quality2)
table(probes_remove)
# probes_remove
# FALSE  TRUE 
# 20215  4132  sum == 24347

# cf. when using "probe_quality"
# FALSE  TRUE 
# 20229  4132  sum == 24361

gse_removed <- gse[!probes_remove,]


# if is.na() === TRUE: remove NAs
bad.sample <- colMeans(is.na(exprs(gse_removed))) > 0.8
bad.gene <- rowMeans(is.na(exprs(gse_removed))) > 0.5
gse_removed <- gse_removed[!bad.gene,!bad.sample]

# Take a look to make a subset
gse_removed
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 20221 features, 65 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM210289 GSM210290 ... GSM210353 (65 total)
# varLabels: title geo_accession ... Sample Group:ch1 (44 total)
# varMetadata: labelDescription
# featureData
# featureNames: GI_10047089-S GI_10047091-S ... thrB (20221 total)
# fvarLabels: PROBEID ACCNUM ... GENENAME (5 total)
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 17520024 
# Annotation: GPL2700
gse <- gse_removed
# ( annotation: OK)



#### Tiding up the data ####
#. Subselect the columns of interest ----
pD.all <- pData(gse)
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

nrow(keep)  # 28 --> nrow(pD) should be 28

kp <- pD$individual %in% keep$individual
gse <- gse[,kp]  # creating the dataset containing only elderly subjects mentioned above
pD <- pD[kp,]
dim(pD)
dim(gse)  # 28  3
# Features  Samples 
# 20221       28 

head(pD, 3); tail(pD, 3)
#           individual gender  id
# GSM210315         A2      F  A2
# GSM210316         A5      F  A5
# GSM210317        A17      M A17
#           individual gender  id
# GSM210351      A12EB      F A12
# GSM210352      A14EB      M A14
# GSM210353      A19EB      M A19

pD$timepoint <- factor(ifelse(str_detect(pD$individual, "EB"), "pre", "post"))
pD$timepoint <- relevel(pD$timepoint, ref = "pre")
levels(pD$timepoint)
# [1] "pre"  "post"

head(pD, 3)
# individual gender  id timepoint
# GSM210315         A2      F  A2      post
# GSM210316         A5      F  A5      post
# GSM210317        A17      M A17      post
nrow(pD) # 28: OK



#### Normalization ####
head(pData(gse)$data_processing, 1)
# [1] "Raw data (see related publication for normalization procedure)"
     # "gse" is not normalized yet

library(beadarray)
norm <- normaliseIllumina(gse)  # quantile normalization
# ?normaliseIllumina
# cf. original study uses normalizeBetweenArrays() of limma.
#     checked with norm_limma <- normalizeBetweenArrays(exprs(norm_limma), method = "quantile").
#     the boxplots (exprs(norm) & exprs(norm_limma)) are identical.

# comparing pre/post normalization
par(mfrow = c(1,2))
boxplot(log2(exprs(gse)), las = 2, cex.axis = 0.7, ylab = expression(log [2](intensity)),
        outline = FALSE, main = "raw")
boxplot(log2(exprs(norm)), las = 2,  cex.axis=0.7, ylab = "",
        outline = FALSE, main = "normalized")



#### Quality control ####
exprs(gse)[1:5, 1:5]
# GSM210315  GSM210316  GSM210317  GSM210318 GSM210319
# GI_10047089-S 5726.23500 7315.88400 6586.29300 8897.02600 7697.3590
# GI_10047091-S   98.85793   98.09677   79.31383   98.20589  108.5538
# GI_10047093-S  483.68300  678.18030  671.47300  596.94870  575.3476
# GI_10047099-S  126.38240  140.13030  119.27290  205.39560  124.6043
# GI_10047103-S  762.17880 1231.18200  625.13700  947.06090  763.3690

exprs(norm)[1:5, 1:5]  # non log-transformed data
# GSM210315  GSM210316  GSM210317  GSM210318 GSM210319
# GI_10047089-S 4749.95064 4813.36339 5272.25471 6981.80814 5693.1232
# GI_10047091-S   91.22922   89.01135   76.44373   84.46068  101.6320
# GI_10047093-S  412.42558  415.60395  548.23660  458.58666  471.0787
# GI_10047099-S  114.17500  111.91306  110.69317  163.25491  115.0227
# GI_10047103-S  634.04039  735.52104  510.47403  732.95137  615.4030



# log2
exp_raw <- log2(exprs(norm))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

#. PCA analysis ----
head(pD, 3)
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     individual = pD$id,
                     gender = pD$gender,
                     timepoint = pD$timepoint
                     )
ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(gender ~ .)



#### Filtering based on intensity (skipped) ####
# filter out lowly expressed genes
# “soft” intensity based filtering here, since this is recommended by the limma
#. Histogram ----
medians <- rowMedians(log2(exprs(norm)))
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
head(pD, 3)
individual <- pD$id
timepoint <- pD$timepoint
gender <- pD$gender


#### lmFit() ####
# make sure to use norm not gse for data
exprs(norm) <- log2(exprs(norm))
for (i in unique(gender)) {
    print(i)
    subjects <- individual[gender == i]
    pre_post <- timepoint[gender == i]
    data <- norm[, gender == i]  # using norm not gse: OK 
    
    design <- model.matrix(~ 0 + pre_post + subjects)
    colnames(design)[1:2] <- c("pre", "post")
    rownames(design) <- subjects
    fit <- lmFit(data, design)
    contrast <- makeContrasts(post-pre, levels = design)
    contr.fit <- eBayes(contrasts.fit(fit, contrast))
    
    #### Results ####
    table <- topTable(contr.fit, coef = 1, number = Inf)
    write.csv(table, file = paste0("res_GSE8479", i, ".csv"))
    print(head(table))
}


#### Visualizing results ####
library(RColorBrewer)
result_files <- list.files(pattern = ".csv")
result_files
# [1] "res_GSE58249normLogF.csv" "res_GSE58249normLogM.csv"
par(mfrow = c(2,2))

for (i in 1:length(result_files)) {
    file <- result_files[i]
    results <- read.csv(file)
    group <- str_sub(file, -5, -5) 
    
    ## histgram ##
    hist(results$P.Value, col = brewer.pal(3, name = "Set2")[1], 
         main = paste(file, "Pval"), xlab  = NULL)
    hist(results$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
         main = paste(file, "adj.Pval"), xlab = NULL)
    
    ## some numbers ##
    cat(file, "\n")
    cat("p < 0.05:", nrow(subset(results, P.Value < 0.05)), "\n")
    cat("p < 0.01:", nrow(subset(results, P.Value < 0.01)), "\n")
    cat("adj.P < 0.05:", nrow(subset(results, adj.P.Val < 0.05)),"\n")
    cat("adj.P < 0.01:", nrow(subset(results, adj.P.Val < 0.01)),"\n\n")
}

# to explore more (ex)
# tail(subset(results, adj.P.Val < 0.05))




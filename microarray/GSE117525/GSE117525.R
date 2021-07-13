# GSE117525
#### Loading packages ####
library(Biobase)
library(oligoClasses)
library(oligo)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(devtools)



#### Loading the data from GEO ####
library(GEOquery)
geo <- getGEO("GSE117525")
length(geo) # 1
gse <- geo[[1]]
gse
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


exprs(gse)[1:3, 1:3]  ## exprs(gse) is already normalized & log-transformed
# GSM3302211 GSM3302212 GSM3302213
# 100009676_at   5.441495   5.146010   5.385128
# 10000_at       6.619549   6.528196   6.275187
# 10001_at       7.193538   7.593082   7.610006



#### Data cleaning ####
# featuredata ----
head(fData(gse), 3)  # (Description & ID can be omitted)
                     # According to the platform page (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL20880)
                     # "SPOT_ID" is "ENTREZ_GENE_ID (ENTREZID)"
#                        ID                                   Description   SPOT_ID
# 100009676_at 100009676_at                        ZBTB11 antisense RNA 1 100009676
# 10000_at         10000_at v-akt murine thymoma viral oncogene homolog 3     10000
# 10001_at         10001_at                    mediator complex subunit 6     10001


# phenodata ----
pD.all <- pData(gse)
head(pD.all, 3)
tail(pD.all, 3)
pD <- pD.all[, c("title", # individual? & supplement
                 "subjectid:ch1", # subject id
                 "supplement:ch1", # drink type including placebo
                 "Sex:ch1", # gender
                 "received training:ch1", # pre/post
                 "group:ch1" # yng/old
                 )] 
names(pD)[2:6] <- c("subjectid", "supplement", "gender", "training", "agegroup")
head(pD, 3)
#                                  title subjectid              supplement gender training         agegroup
# GSM3302211          G084_A01_15_DSMT22    DSMT22 placebo / no supplement      M       no Young (baseline)
# GSM3302212 G084_A02_34_06OSN_MEAT_MILK     06OSN placebo / no supplement      M       no Young (baseline)
# GSM3302213       G084_A03_61_4007_PROM      4007 placebo / no supplement      M       no Frail (baseline)
dim(pD)  # 259   6
length(unique(pD$subjectid))  # 187



# Tiding up phenodata ----
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
# FALSE  TRUE 
# 115   144

pD <- pD[keep,]
gse <- gse[,keep]
dim(pD)  # 144   6
dim(gse)
# Features  Samples 
# 19654      144


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


#. group ----
pD$group <- paste0(pD$gender, str_sub(toupper(pD$health), 1, 1))


#. supplement ----
# excluding "protein" induced subjects
keep <- grep("placebo", pD$supplement)
pD <- pD[keep,]
gse <- gse[,keep]
dim(pD)  # 70  9
dim(gse)
# Features  Samples 
# 19654      70
levels(as.factor(pD$supplement))
# [1] "placebo / no supplement"

# renaming the label
pD$supplement <- ifelse(str_detect(pD$supplement, "placebo"), "placebo", "NA")
levels(as.factor(pD$supplement))
# [1] "placebo"


# checking the final pData ----
head(pD, 3)
#                                 title subjectid supplement gender training         agegroup timepoint health group
# GSM3302215 G084_A05_109_8011_PRE_PROM    S.8011    placebo      F      yes Frail (baseline)       pre  frail    FF
# GSM3302216 G084_A06_117_8027_PRE_PROM    S.8027    placebo      F      yes Frail (baseline)       pre  frail    FF
# GSM3302218 G084_A08_155_8082_PRE_PROM    S.8082    placebo      M      yes Frail (baseline)       pre  frail    MF

str(pD)
# 'data.frame':	70 obs. of  9 variables:
# $ title     : chr  "G084_A05_109_8011_PRE_PROM" "G084_A06_117_8027_PRE_PROM" "G084_A08_155_8082_PRE_PROM" "G084_A09_165_8093_PRE_PROM" ...
# $ subjectid : chr  "S.8011" "S.8027" "S.8082" "S.8093" ...
# $ supplement: chr  "placebo" "placebo" "placebo" "placebo" ...
# $ gender    : chr  "F" "F" "M" "F" ...
# $ training  : chr  "yes" "yes" "yes" "yes" ...
# $ agegroup  : chr  "Frail (baseline)" "Frail (baseline)" "Frail (baseline)" "Frail (baseline)" ...
# $ timepoint : Factor w/ 2 levels "pre","post": 1 1 1 1 2 2 2 2 1 1 ...
# $ health    : Factor w/ 2 levels "healthy","frail": 2 2 2 2 2 2 2 2 2 2 ...
# $ group     : chr  "FF" "FF" "MF" "FF" ...




#### Normalization ####
head(pData(gse)$data_processing, 1)
# [1] "Expression estimates were calculated applying the robust microarray analysis (RMA) algorithm."

oligo::boxplot(gse, target = "core",
               main = "Boxplot of log2-intensitites for the normalized data",
               las = 2,
               outline = FALSE,
               cex.axis=0.7)

exprs(gse[1:3, 1:3])
#              GSM3302215 GSM3302216 GSM3302218
# 100009676_at   4.943804   5.450881   5.379767
# 10000_at       6.339496   6.636859   6.446695
# 10001_at       7.311173   7.571778   7.314539




#### Quality control ####
#. PCA analysis ----
myexp <- exprs(gse)
PCA <- prcomp(t(myexp), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     individual = pD$subjectid,
                     gender = pD$gender,
                     timepoint = pD$timepoint,
                     health = pD$health,
                     group = pD$group
                     )

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(health ~ gender) +
    ggtitle("PCA_Norm timepoint & gender & health") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%"))

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(. ~ gender) +
    ggtitle("PCA_Norm timepoint & gender") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%"))



#### Filtering based on intensity ####
# filter out lowly expressed genes
# “soft” intensity based filtering here, since this is recommended by the limma
#. Histogram ----
medians <- rowMedians(exprs(gse))
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
no_of_samples <- table(paste0(pD$gender, "_", pD$timepoint))
no_of_samples
# F_post  F_pre M_post  M_pre 
# 14     13     21     22

samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(exprs(gse), 1,
                           function(x){
                               sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)
# FALSE  TRUE 
# 5747 13907

manfiltered <- subset(gse, idx_man_threshold)
dim(exprs(gse))  # 19654    70
dim(exprs(manfiltered))  # 13907    70



#### Annotation of the transcript clusters ####
library(hugene11sttranscriptcluster.db)
# columns(hugene11sttranscriptcluster.db)
# keytypes(hugene11sttranscriptcluster.db)
# head(keys(hugene11sttranscriptcluster.db, keytype="ENTREZID"))
# ls("package:hugene11sttranscriptcluster.db")
fData(gse)[3]
gse_entrezid <- as.character(fData(manfiltered)$SPOT_ID)  # fData(gse)$SPOT_ID == ENTREZID
anno <- AnnotationDbi::select(hugene11sttranscriptcluster.db,
                              keys = (gse_entrezid),
                              columns = c("SYMBOL", "GENENAME"),
                              keytype = "ENTREZID")
head(anno)
#    ENTREZID     SYMBOL                           GENENAME
# 1 100009676 ZBTB11-AS1             ZBTB11 antisense RNA 1
# 2     10000       AKT3      AKT serine/threonine kinase 3
# 3     10001       MED6         mediator complex subunit 6
# 4 100033413 SNORD116-1 small nucleolar RNA, C/D box 116-1
# 5 100033414 SNORD116-2 small nucleolar RNA, C/D box 116-2
# 6 100033416 SNORD116-4 small nucleolar RNA, C/D box 116-4


table(is.na(anno))
# FALSE  TRUE 
# 41601   120  


anno <- subset(anno, !is.na(SYMBOL))
sum(is.na(anno$SYMBOL)) # 0
nosymbols <- !(gse_entrezid %in% anno$ENTREZID)
table(nosymbols)
# FALSE  TRUE 
# 13847    60

withsymbols <- subset(manfiltered, !nosymbols)
dim(withsymbols)
# Features  Samples 
# 13847       70


#. Removing multiple mappings ----
# group by PROBEID
anno_grouped <- group_by(anno, ENTREZID)
head(anno_grouped)

# summarize
anno_summarized <-
    dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized)

# filter to extract cluster multiple gene that map to multiple gene symbols
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
head(anno_filtered)
# A tibble: 0 x 2
# … with 2 variables: ENTREZID <chr>, no_of_matches <int>

probe_stats <- anno_filtered
nrow(probe_stats) # 0: clusters that map to multiple gene symbols → remove
# remove above IDs(probe_stats)
ids_to_exlude <- (featureNames(withsymbols) %in% probe_stats$ENTREZID)
table(ids_to_exlude)
# FALSE
# 13847

final <- subset(withsymbols, !ids_to_exlude)
validObject(final)
dim(final)  # 13847  70 


# also exclude them from the feature data anno
head(anno, 3)

# generate a column PROBEID in fData(final) and
# assign the row names of fData(final) to it:
fData(final)$PROBEID <- featureNames(final)
fData(final)$ENTREZID <- anno[match(fData(final)$SPOT_ID, anno$ENTREZID), 1]  # adding anno$ENTREZID
head(fData(final))

# left-join fData(final) with anno
fData(final) <- left_join(fData(final), anno, by = "ENTREZID")
head(fData(final))

fData(final) <- fData(final)[, c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME")]
rownames(fData(final)) <- fData(final)$PROBEID
head(fData(final))

validObject(final)
dim(final)
# Features  Samples 
# 13847       70 



#### Differential Expression nalysis ####
head(pD, 3)
individual <- pD$subjectid
timepoint <- pD$timepoint
gender <- pD$gender
# health <- pD$health
# group <- pD$group

#### lmFit() ####
exprs(final)[1:3, 1:3]  # checking if log-transformed
# exprs(final) <- log2(exprs(final))
for (i in unique(gender)) {
    print(i)
    subjects <- individual[gender == i]
    pre_post <- timepoint[gender == i]
    data <- final[, gender == i]
    
    design <- model.matrix(~ 0 + pre_post + subjects)
    colnames(design)[1:2] <- c("pre", "post")
    rownames(design) <- subjects
    fit <- lmFit(data, design)
    contrast <- makeContrasts(post - pre, levels = design)
    contr.fit <- eBayes(contrasts.fit(fit, contrast))
    
    #### Results ####
    table <- topTable(contr.fit, coef = 1, number = Inf)
    write.csv(table, file = paste0("res_GSE117525re", i, ".csv"))
    print(head(table))
}


#### Visualizing results ####
library(RColorBrewer)
result_files <- list.files(pattern = "re..csv")
result_files
# [1] "res_GSE117525reF.csv" "res_GSE117525reM.csv"

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




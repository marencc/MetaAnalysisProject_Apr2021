#### GSE28422: Resistance ####
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


#### Loading the data from GEO ####
library(GEOquery)
geo <- getGEO("GSE28422")
length(geo) # 1
gse <- geo[[1]]
gse
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 54675 features, 110 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM702357 GSM702358 ... GSM702466 (110 total)
# varLabels: title geo_accession ... training state:ch1 (39 total)
# varMetadata: labelDescription
# featureData
# featureNames: 1007_s_at 1053_at ... AFFX-TrpnX-M_at (54675 total)
# fvarLabels: ID GB_ACC ... Gene Ontology Molecular Function (16 total)
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 22302958 
# Annotation: GPL570


exprs(gse)[1:3, 1:3]
#           GSM702357 GSM702358 GSM702359
# 1007_s_at   1161.83  1071.500  1957.930
# 1053_at     1747.50  1478.970  3048.850
# 117_at       391.97   470.198   222.002



#### Data cleaning ####
# featuredata ----
head(fData(gse), 3)
names(fData(gse))
fData(gse) <- fData(gse)[,c("ID",  # ex. 1007_s_at
                            "GB_ACC",
                            "Gene Symbol",
                            "ENTREZ_GENE_ID"
                            )]
names(fData(gse)) <- c("PROBEID", "ACCNUM", "SYMBOL","ENTREZID")
head(fData(gse), 3)
#             PROBEID ACCNUM           SYMBOL          ENTREZID
# 1007_s_at 1007_s_at U48705 DDR1 /// MIR4640 780 /// 100616237
# 1053_at     1053_at M87338             RFC2              5982
# 117_at       117_at X51757            HSPA6              3310


# phenodata ----
pD.all <- pData(gse)
head(pD.all, 3)
tail(pD.all, 3)

pD <- pD.all[, c("title", # identification?
                 "age:ch1",
                 "gender:ch1", # gender
                 "time point:ch1", # timepoint
                 "training state:ch1",
                 "description" # id?
                 )]
head(pD, 3)
unique(pD$title)

#. Renaming ----
names(pD)[2:5] <- c("age", "gender", "timepoint", "training")
pD$description <- gsub("EA08017_.*_AGEX_OF_|EA08017_.*_AGEX_YF_|EA08017_.*_AGEX_OM_|EA08017_.*_AGEX_YM_",
                       "", pD$description)
pD$description <- gsub("_H133\\+.CHP", "", pD$description)
pD$id <- paste0("S", gsub("T._", "", pD$description))

pD$biopsy <- str_sub(pD$title, 1, 2)  # ex. T1, T2, T3, T4
# this analysis needs "T1" & "T3"
keep <- pD$biopsy %in% c("T1", "T3")
pD <- pD[keep,]
gse <- gse[,keep]
dim(gse)
# Features  Samples 
# 54675       54

pD$timepoint <- factor(ifelse(str_detect(pD$biopsy, "T1"), "pre", "post"))
pD$timepoint <- relevel(pD$timepoint, ref = "pre")
levels(pD$timepoint)

# excluding subjects that do not have both pre & post data
complete_subjects <- pD %>% group_by(id) %>% 
    summarise(complete_subjects = n_distinct(biopsy)) %>% 
    filter(complete_subjects == 2)
keep <- pD$id %in% complete_subjects$id
pD <- pD[keep,]
gse <- gse[,keep]
dim(pD)  # 52  8
dim(gse)
# Features  Samples 
# 54675       52

pD$group <- paste0(str_sub(pD$gender, 1, 1), str_sub(pD$age, 1, 1))
head(pD, 3)
#                               title   age gender timepoint  training description id biopsy group
# GSM702357 T1_Pre_Male_Young (81373) Young   Male       pre Untrained        T1_4 S4     T1    MY
# GSM702358 T1_Pre_Male_Young (81379) Young   Male       pre Untrained        T1_6 S6     T1    MY
# GSM702359 T1_Pre_Male_Young (81390) Young   Male       pre Untrained        T1_9 S9     T1    MY



#### Normalization (skipped) ####
head(pData(gse)$data_processing, 1)
# [1] "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using 
# Affymetrix default analysis settings and global scaling as normalization method.
# The trimmed mean target intensity of each array was arbitrarily set to 1500."

#### checking data intensities
oligo::boxplot(log2(exprs(gse)), target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               outline = FALSE,
               las = 2,
               cex.axis=0.7)
# memo: based on the boxplot, the data seems to be already normalized.



#### Quality control of the raw data ####
# Take a look
exprs(gse)[1:5, 1:5]
#            GSM702357 GSM702358 GSM702359 GSM702360  GSM702361
# 1007_s_at 1161.83000 1071.5000 1957.9300 1631.5800 1257.48000
# 1053_at   1747.50000 1478.9700 3048.8500 2239.7700 1964.19000
# 117_at     391.97000  470.1980  222.0020   94.4564   75.25550
# 121_at     129.94200  162.2530  303.0210  163.4050  194.74900
# 1255_g_at    6.23465   30.2465   11.2966   19.5840    8.09821



# log2
exp_raw <- log2(exprs(gse))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

#. PCA plot of the log ----
head(pD, 3)
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     individual = pD$id,
                     timepoint = pD$biopsy,
                     group = pD$group,
                     gender = pD$gender
                     )

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(. ~ gender) +
    ggtitle("PCA: timepoint & gender") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%"))

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(group ~ .) +
    ggtitle("PCA: timepoint & gender & age") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%"))




#### Filtering based on intensity ####
# filter out lowly expressed genes
#. Histogram ----
medians <- rowMedians(log2(exprs(gse)))
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

#set a manual threshold to exclude genes, but not too many
# man_threshold <- round(mean(medians))  # 8
man_threshold <- 5
abline(v = man_threshold, col = "coral4", lwd = 2)

#. Excluding genes that below the threshold ----
no_of_samples <- table(pD$group)
no_of_samples
# FO FY MO MY 
# 12 12 12 16

samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(log2(exprs(gse)), 1,
                           function(x){
                               sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)
# idx_man_threshold
# FALSE  TRUE 
# 5109 49566

manfiltered <- subset(gse, idx_man_threshold)

dim(exprs(gse))  # 54675    52
dim(exprs(manfiltered))  # 49566    52



#### Annotation of the transcript clusters ####
library(hgu133plus2.db)
anno <- AnnotationDbi::select(hgu133plus2.db,
                              keys = (featureNames(manfiltered)),
                              columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                              keytype = "PROBEID")
table(is.na(anno))
# FALSE   TRUE 
# 185169  26523


anno <- subset(anno, !is.na(SYMBOL))
sum(is.na(anno$SYMBOL)) # 0
nosymbols <- !(featureNames(manfiltered) %in% anno$PROBEID)
table(nosymbols)
# FALSE  TRUE 
# 40725  8841

withsymbols <- subset(manfiltered, !nosymbols)
dim(withsymbols)
# Features  Samples 
# 40725       52


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
nrow(probe_stats) # 2058: clusters that map to multiple gene symbols → remove

# remove above IDs(probe_stats)
ids_to_exlude <- (featureNames(withsymbols) %in% probe_stats$PROBEID)
table(ids_to_exlude)
# FALSE  TRUE 
# 38667  2058 

final <- subset(withsymbols, !ids_to_exlude)
validObject(final)
dim(final)  # 38667  52


# also exclude them from the feature data anno
# generate a column PROBEID in fData(final) and
# assign the row names of fData(final) to it:
fData(final)$PROBEID <- rownames(fData(final))

# left-join fData(final) with anno
fData(final) <- left_join(fData(final), anno)

rownames(fData(final)) <- fData(final)$PROBEID
validObject(final)
dim(final)
# Features  Samples 
# 38667       52



#### Differential Expression Analysis ####
head(pD, 3)
individual <- pD$id
timepoint <- pD$timepoint
group <- pD$group

#### lmFit() ####
exprs(final) <- log2(exprs(final))
for (i in unique(group)) {
    print(i)
    subjects <- individual[group == i]
    pre_post <- timepoint[group == i]
    data <- final[, group == i]
    
    design <- model.matrix(~ 0 + pre_post + subjects)
    colnames(design)[1:2] <- c("pre", "post")
    rownames(design) <- subjects
    fit <- lmFit(data, design)
    contrast <- makeContrasts(post - pre, levels = design)
    contr.fit <- eBayes(contrasts.fit(fit, contrast))
    
    #### Results ####
    table <- topTable(contr.fit, coef = 1, number = Inf)
    write.csv(table, file = paste0("res_GSE28422re", i, ".csv"))
    print(head(table))
}


#### Visualizing results ####
library(RColorBrewer)
result_files <- list.files(pattern = ".csv")
result_files
# [1] "res_GSE28422FO.csv" "res_GSE28422FY.csv" "res_GSE28422MO.csv" "res_GSE28422MY.csv"
par(mfrow = c(2,2))

for (i in 1:length(result_files)) {
    file <- result_files[i]
    results <- read.csv(file)
    group <- str_sub(file, -6, -5) 
    
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
    cat("adj.P < 0.05:", nrow(subset(results, adj.P.Val < 0.01)),"\n\n")
}

# to explore more (ex)
# tail(subset(results, adj.P.Val < 0.05))


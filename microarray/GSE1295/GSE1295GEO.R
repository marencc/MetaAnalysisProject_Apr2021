# GSE1295
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1295
# Journal: https://journals.physiology.org/doi/full/10.1152/japplphysiol.00331.2004
# Affymetrix Human Genome U133A (GPL96) --> hgu133a.db


#### Loading packages ####
library(Biobase)
library(oligoClasses)
library(hgu133a.db)
# library(pd.hg.u133a)
library(oligo)
# library(arrayQualityMetrics)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(devtools)


#### Importing raw data from GEO ####
library(GEOquery)
geo <- getGEO("GSE1295")
length(geo)
names(geo)
# [1] "GSE1295-GPL8300_series_matrix.txt.gz" "GSE1295-GPL96_series_matrix.txt.gz" 
    ### For this analysis, we need a data from platform GPL96
gse <- geo[[2]]
gse
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


exprs(gse)[1:3, 1:3]
#            GSM19136  GSM19148  GSM19149
# 1007_s_at 249.13837 249.63313 261.45288
# 1053_at    78.84882  18.68847  13.88993
# 117_at     63.93585  45.46224  86.41307


head(fData(gse), 3)
names(fData(gse))

#### Data cleaning ####
# featuredata ----
head(fData(gse), 3)
fData(gse) <- fData(gse)[,c("ID",  # ex. 1007_s_at
                            "Gene Symbol",
                            "ENTREZ_GENE_ID"
                            )]
names(fData(gse)) <- c("PROBEID", "SYMBOL","ENTREZID")
head(fData(gse), 3)


# phenodata ----
pD.all <- pData(gse)
head(pD.all, 3)
tail(pD.all, 3)

pD <- pD.all[, c("title", # identification?
                 "description" # identification, gender, timepoint
                 )]

pD$gender <- as.factor(ifelse(str_detect(pD$title, "-M"), "M", "F")) # gender
pD$title <- str_replace(pD$title, "PSS-", "")
pD$title <- str_replace(pD$title, "aUA-s2", "")
pD$timepoint <- ifelse(str_detect(pD$title, "0h"), "pre", (ifelse(str_detect(pD$title, "24h"), "post", "NA")))
pD$group <- paste(pD$gender, pD$timepoint, sep = "_")
pD$id <- paste0("S_", pD$gender, str_sub(pD$title, -1, -1))
head(pD, 3)
#           title                                               description gender timepoint group   id
# GSM19136  F0h-1                                                       N/A      F       pre F_pre S_F1
# GSM19148 M96h-5 LC ID 208 Male Age 40 to 65  25 mg Needle biopsy 96 hours      M        NA  M_NA S_M5
# GSM19149 M96h-4 LC ID 254 Male Age 40 to 65  25 mg Needle biopsy 96 hours      M        NA  M_NA S_M4


# subsetting ----
#. keep only pre/post exercise ----
keep <- grep("pre|post", pD$timepoint)
pD <- pD[keep,]
gse <- gse[,keep]
pD
#           title                                                        description gender timepoint  group   id
# GSM19136  F0h-1                                                                N/A      F       pre  F_pre S_F1
# GSM19158  F0h-2          LC ID 082 Female Age 40 to 65  25 mg Needle biopsy 0 hour      F       pre  F_pre S_F2
# GSM19159  F0h-3          LC ID 093 Female Age 40 to 65  25 mg Needle biopsy 0 hour      F       pre  F_pre S_F3
# GSM19160 F24h-3         LC ID 093 Female Age 40 to 65  25 mg Needle biopsy 24 hour      F      post F_post S_F3
# GSM19161 F24h-1                                                                N/A      F      post F_post S_F1
# GSM19162 F24h-2         LC ID 082 Female Age 40 to 65  25 mg Needle biopsy 24 hour      F      post F_post S_F2
# GSM19170  M0h-4           LC ID 254 Male Age 40 to 65  25 mg Needle biopsy 0 hours      M       pre  M_pre S_M4
# GSM19171  M0h-5           LC ID 208 Male Age 40 to 65  25 mg Needle biopsy 0 hours      M       pre  M_pre S_M5
# GSM19172  M0h-1            LC ID 079 Male Age 40 to 65  25 mg Needle biopsy 0 hour      M       pre  M_pre S_M1
# GSM19173 M24h-1           LC ID 079 Male Age 40 to 65  25 mg Needle biopsy 24 hour      M      post M_post S_M1
# GSM19174 M24h-5          LC ID 208 Male Age 40 to 65  25 mg Needle biopsy 24 hours      M      post M_post S_M5
# GSM20659 M24h-4 LC ID 254 Male needle biopsy within 24 hours of last exercise bout      M      post M_post S_M4

dim(pD)  # 12  6
dim(gse)
# Features  Samples 
# 22283       12


#### normalization ####
# memo: no information about normalization in original phenodata

#### checking data intensities
oligo::boxplot(log2(exprs(gse)), target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               # outline = FALSE,
               las = 2,
               cex.axis=0.7)
# memo: based on this boxplot, the data seems to be already normalized
# cf. (from paper, section "Gene expression data analysis")
#       "Background and scaled noise were similarly averaged for all chips before analysis."
#       https://journals.physiology.org/doi/full/10.1152/japplphysiol.00331.2004


#### Quality control of the raw data ####
#. PCA analysis ----
normData <- gse
myexp <- log2(exprs(normData))
PCA <- prcomp(t(myexp), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

head(pD, 3)
dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     individual = pD$id,
                     gender = pD$gender,
                     timepoint = pD$timepoint,
                     group = pD$group
                     )

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
medians <- rowMedians(log2(exprs(gse)))
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

#set a manual threshold to exclude genes, but not too many
# man_threshold <- round(mean(medians))  # 6
man_threshold <- 4
abline(v = man_threshold, col = "coral4", lwd = 2)

#. Excluding genes that below the threshold ----
# Transcripts that do not have intensities larger than the threshold in
# at least as many arrays as the smallest experimental group are excluded.
# get a list:
no_of_samples <- table(paste(pD$gender, pD$timepoint, sep = "_"))
no_of_samples

samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(exprs(gse), 1,
                           function(x){
                               sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)


manfiltered <- subset(gse, idx_man_threshold)

dim(exprs(gse))  # 22283    12
dim(exprs(manfiltered))  # 22233    12



#### Annotation of the transcript clusters ####
library(hgu133a.db)
anno <- AnnotationDbi::select(hgu133a.db,
                              keys = (featureNames(manfiltered)),
                              columns = c("SYMBOL", "GENENAME"),
                              keytype = "PROBEID")
table(is.na(anno))
# FALSE  TRUE 
# 70836  2406


anno <- subset(anno, !is.na(SYMBOL))
sum(is.na(anno$SYMBOL)) # 0
nosymbols <- !(featureNames(manfiltered) %in% anno$PROBEID)
table(nosymbols)
# FALSE  TRUE 
# 21030  1203

withsymbols <- subset(manfiltered, !nosymbols)
dim(withsymbols)
# Features  Samples 
# 21030       12


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
# 19808  1222

final <- subset(withsymbols, !ids_to_exlude)
validObject(final)
dim(final)  # 19808    12 


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
# 19808       12



#### Differential Expression Analysis ####
head(pD, 3)
individual <- pD$id
timepoint <- pD$timepoint
gender <- pD$gender

#### lmFit() ####
for (i in unique(gender)) {
    print(i)
    subjects <- individual[gender == i]
    pre_post <- timepoint[gender == i]
    data <- final[, gender == i]
    
    design <- model.matrix(~ 0 + pre_post + subjects)
    colnames(design)[1:2] <- c("pre", "post")
    rownames(design) <- subjects
    fit <- lmFit(data, design)
    contrast <- makeContrasts(post-pre, levels = design)
    contr.fit <- eBayes(contrasts.fit(fit, contrast))
    
    #### Results ####
    table <- topTable(contr.fit, coef = 1, number = Inf)
    write.csv(table, file = paste0("res_GSE1295", i, ".csv"))
}


#### Visualizing results ####
library(RColorBrewer)
result_files <- list.files(pattern = ".csv")
result_files
# [1] "res_GSE1295F.csv" "res_GSE1295M.csv"
par(mfrow = c(2,2))

for (i in 1:length(result_files)) {
    file <- result_files[i]
    results <- read.csv(file)
    group <- str_sub(file, -5, -5) 
    
    ## histgram ##
    hist(results$P.Value, col = brewer.pal(3, name = "Set2")[1], 
         main = paste(group, "Pval"), xlab  = NULL)
    hist(results$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
         main = paste(group, "Pval"), xlab = NULL)
    
    ## some numbers ##
    cat(group, "\n")
    cat("p < 0.05:", nrow(subset(results, P.Value < 0.05)), "\n")
    cat("adj.P < 0.05:", nrow(subset(results, adj.P.Val < 0.05)),"\n\n")
}

# to explore more (ex)
# tail(subset(results, adj.P.Val < 0.05))


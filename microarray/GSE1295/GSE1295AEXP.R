#### GSE1295: Aerobics, Overweight Men & Women  ####
# 3 male, 3 female
# metabolic syndrom
# biopsy: 9M later (24h, 96h: 4D, 336h: 14D)  ---> compare: base & 24h (9M later)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1295
# https://journals.physiology.org/doi/full/10.1152/japplphysiol.00331.2004
# E-GEOD-1295

#### List of packages required for the workflow ####
# General Bioconductor packages
library(Biobase)
library(oligoClasses)

# Annotation and data import packages
library(ArrayExpress)
library(pd.hg.u133a)
library(hgu133a.db)

# Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

# Analysis and statistics packages
library(limma)
# library(topGO)
# library(ReactomePA)
library(clusterProfiler)

# Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)

# Helpers
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
library(devtools)


#### Download the raw data from ArrayExpress ####
raw_data_dir <- tempdir()

anno_AE <- getAE("E-GEOD-1295", path = raw_data_dir, type = "raw")


#### Import of annotation data and microarray expression data as “ExpressionSet” ####
sdrf_location <- file.path(raw_data_dir, "E-GEOD-1295.sdrf.txt")
SDRF <- read.delim(sdrf_location)
?read.delim
pData(SDRF)
sampleNames(SDRF)

rownames(SDRF) <- SDRF$Array.Data.File # sample names is stored here
SDRF <- AnnotatedDataFrame(SDRF) # to create an ExpressionSet later on

#### Reading in raw data ####
celfiles <- list.files(raw_data_dir, pattern = "CEL$", full = TRUE)
head(celfiles)
chipTypes <- sapply(celfiles, oligo:::getCelChipType, useAffyio = TRUE)
chipTypes  ## Need to extract "HG-U133A"

dir.create(paste(raw_data_dir, "HG-U133A", sep = "/"))
u133adir <- paste(raw_data_dir, "HG-U133A", sep = "/")

celnames <- list.files(raw_data_dir, pattern = "CEL$")
celnames

for (i in 1:length(celfiles)) {
    celname = celnames[i]
    celfile = celfiles[i]
    if (oligo:::getCelChipType(celfile, useAffyio = TRUE) == "HG-U133A") {
        file.copy(from = celfile, to = paste(u133adir, celname, sep = "/"), overwrite = TRUE)
        file.remove(celfile)
    } else {    
        file.remove(celfile)
    }
}

celfiles <- list.files(u133adir, full = TRUE)
celnames <- list.files(u133adir)
pData(SDRF) <- pData(SDRF)[celnames,]  ## Now, SDRF has rows that this analysis need

#### Create Expression Set ####
# Specified SDRF file we created earlier as phenoData.
# Thus, we had to make sure to import the CEL files in the order corresponds to the SDRF table.
# To enforce this, we used the column Array.Data.File of the SDRF table as the filenames argument.
raw_data <- oligo::read.celfiles(filenames = file.path(u133adir, SDRF$Array.Data.File),
                                 verbose = FALSE, phenoData = SDRF)
stopifnot(validObject(raw_data))

# Take a look to make a subset
head(Biobase::pData(raw_data))

#. Subselect the columns of interest ----
pD.all <- pData(raw_data)
pD <- pD.all[, c("Hybridization.Name", # identification?
                 "Characteristics..Sex.", # gender
                 "Factor.Value..sampling.time.point.", # timepoint
                 "Unit..TimeUnit.", # unit (day)
                 "Description" # identification?
                 )]

names(pD) <- c("names", "gender", "day", "unit", "description")

pD$gender <- as.factor(ifelse(str_detect(pD$title, "-M"), "male", "female")) # timepoint
pD$names <- str_replace(pD$names, "PSS-", "")
pD$names <- str_replace(pD$names, "-s2", "")
pD$timepoint <- paste0("d", pD$day)
pD$id <- paste0(str_sub(pD$gender, 1, 1), str_sub(pD$names, -4, -1))
# sum(table(unique(paste0(str_sub(pD$gender, 1, 1), str_sub(pD$names, -4, -1)))))  == 6 (male3, female3)
pD


#. keep only d0(base) & d1(9M later) ----
keep <- grep("d0|d1$", pD$timepoint)
pD <- pD[keep,]
raw_data <- raw_data[,keep]
pD
#                  names gender day unit                                                        description timepoint    id
# GSM19136.CEL  F0h-1aUA female   0  day                                                                           d0 f1aUA
# GSM19158.CEL  F0h-2aUA female   0  day          LC ID 082 Female Age 40 to 65  25 mg Needle biopsy 0 hour        d0 f2aUA
# GSM19159.CEL  F0h-3aUA female   0  day          LC ID 093 Female Age 40 to 65  25 mg Needle biopsy 0 hour        d0 f3aUA
# GSM19160.CEL F24h-3aUA female   1  day         LC ID 093 Female Age 40 to 65  25 mg Needle biopsy 24 hour        d1 f3aUA
# GSM19161.CEL F24h-1aUA female   1  day                                                                           d1 f1aUA
# GSM19162.CEL F24h-2aUA female   1  day         LC ID 082 Female Age 40 to 65  25 mg Needle biopsy 24 hour        d1 f2aUA
# GSM19170.CEL  M0h-4aUA   male   0  day           LC ID 254 Male Age 40 to 65  25 mg Needle biopsy 0 hours        d0 m4aUA
# GSM19171.CEL  M0h-5aUA   male   0  day           LC ID 208 Male Age 40 to 65  25 mg Needle biopsy 0 hours        d0 m5aUA
# GSM19172.CEL  M0h-1aUA   male   0  day                                                                           d0 m1aUA
# GSM19173.CEL M24h-1aUA   male   1  day           LC ID 079 Male Age 40 to 65  25 mg Needle biopsy 24 hour        d1 m1aUA
# GSM19174.CEL M24h-5aUA   male   1  day          LC ID 208 Male Age 40 to 65  25 mg Needle biopsy 24 hours        d1 m5aUA
# GSM20659.CEL M24h-4aUA   male   1  day LC ID 254 Male needle biopsy within 24 hours of last exercise bout        d1 m4aUA

dim(pD)  # 12  7
dim(raw_data)
# Features  Samples 
# 506944       12


#### Quality control of the raw data ####
# Take a look
Biobase::exprs(raw_data)[1:5, 1:5]
#   GSM19136.CEL GSM19158.CEL GSM19159.CEL GSM19160.CEL GSM19161.CEL
# 1        167.0        258.0        448.0        545.0        178.8
# 2       9706.3      12328.0       4292.8      13747.3       9071.0
# 3        193.8        251.3        499.5        579.0        243.5
# 4      10016.5      12817.3       4484.0      13825.8       9388.0
# 5        110.8        147.8        156.5        168.5        134.0



# log2
exp_raw <- log2(Biobase::exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

# PCA plot of the log
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     timepoint = pD$timepoint,
                     gender = pD$gender,
                     id = pD$id)

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(shape = id, colour = timepoint)) +
    ggtitle("PCA plot of the log-transformed raw expression data") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    # theme(plot.title = element_text(hjust = 0.5))+
    # coord_fixed(ratio = sd_ratio) +
    # scale_shape_manual(values = c(4,15)) +
    scale_color_manual(values = c("darkorange2", "dodgerblue4"))

# Intensity boxplots of the log2
oligo::boxplot(raw_data, target = "core",
               main = "Boxplot of log2-intensitites for the raw data")

# more elaborate quality control plots
arrayQualityMetrics(expressionset = raw_data,
                    outdir = tempdir(),
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = c("Factor.Value.disease.", "Factor.Value.phenotype."))
dev.off()

#### Background adjustment, calibration, summarization and annotation ####
head(ls("package:hgu133a.db"))
# [1] "hgu133a"          "hgu133a_dbconn"   "hgu133a_dbfile"   "hgu133a_dbInfo"   "hgu133a_dbschema"
# [6] "hgu133a.db"

#### Normalization: rma() ####
# normalize = FALSE を入れないver (targetをベースにnormalizeされる)
eset_norm <- oligo::rma(raw_data)  ## removed "target = 'core' " due to error: unused argument (target = "core")
eset_norm
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 22283 features, 12 samples 
# element names: exprs 
# protocolData
# rowNames: GSM19136.CEL GSM19158.CEL ... GSM20659.CEL (12 total)
# varLabels: exprs dates
# varMetadata: labelDescription channel
# phenoData
# rowNames: GSM19136.CEL GSM19158.CEL ... GSM20659.CEL (12 total)
# varLabels: Source.Name Characteristics..OrganismPart. ...
# Comment..Derived.ArrayExpress.FTP.file. (25 total)
# varMetadata: labelDescription channel
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation: pd.hg.u133a


#. PCA analysis ----
exp <- Biobase::exprs(eset_norm)
PCA <- prcomp(t(exp), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     timepoint = pD$timepoint,
                     gender = pD$gender,
                     id = pD$id)

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(shape = id, colour = timepoint)) +
    ggtitle("PCA plot of the calibrated, summarized data") +
    facet_grid(. ~ gender) + 
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    scale_color_manual(values = c("darkorange2", "dodgerblue4"))


#. Heatmap clustering analysis (skipped) ----
phenotype_names <- ifelse(str_detect(pData(eset_norm)$Factor.Value.phenotype.,
                                     "non"), "non_infl.", "infl.")

disease_names <- ifelse(str_detect(pData(eset_norm)$Factor.Value.disease.,
                                   "Crohn"), "CD", "UC")

annotation_for_heatmap <-
    data.frame(Phenotype = phenotype_names, Disease = disease_names)

row.names(annotation_for_heatmap) <- row.names(pData(eset_norm))

# change the default distance calculation method(Euclidean) to "manhattan"
dists <- as.matrix(dist(t(exp), method = "manhattan"))

rownames(dists) <- row.names(pData(eset_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

ann_colors <- list(
    Phenotype = c(non_infl. = "chartreuse4", infl. = "burlywood3"),
    Disease = c(CD = "blue4", UC = "cadetblue2")
)

# heatmap
pheatmap(dists, col = (hmcol),
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE,
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE),
                           max(dists, na.rm = TRUE)),
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the calibrated samples")

# note: in order to stay as close as possible to the original paper,
# we continue with the complete set of samples.

#### Filtering based on intensity ####
# filter out lowly expressed genes
# “soft” intensity based filtering here, since this is recommended by the limma
#. Histogram ----
medians <- rowMedians(Biobase::exprs(eset_norm))

hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

#set a manual threshold to exclude genes, but not too many
man_threshold <- 6.5 

hist_res <- hist(medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)

#. Exclude genes that below the threshold ----
# Transcripts that do not have intensities larger than the threshold in
# at least as many arrays as the smallest experimental group are excluded.
# get a list:
no_of_samples <- table(paste0(pD$timepoint, "_", pD$gender))
no_of_samples
# d0_female   d0_male d1_female   d1_male 
# 3         3         3         3

samples_cutoff <- min(no_of_samples)  # 3


idx_man_threshold <- apply(Biobase::exprs(eset_norm), 1,
                           function(x){
                               sum(x > man_threshold) >= samples_cutoff
                               }
                           )
table(idx_man_threshold)

manfiltered <- subset(eset_norm, idx_man_threshold)
manfiltered
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 14678 features, 12 samples 
# element names: exprs 
# protocolData
# rowNames: GSM19136.CEL GSM19158.CEL ... GSM20659.CEL (12 total)
# varLabels: exprs dates
# varMetadata: labelDescription channel
# phenoData
# rowNames: GSM19136.CEL GSM19158.CEL ... GSM20659.CEL (12 total)
# varLabels: Source.Name Characteristics..OrganismPart. ...
# Comment..Derived.ArrayExpress.FTP.file. (25 total)
# varMetadata: labelDescription channel
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation: pd.hg.u133a



#### Annotation of the transcript clusters ####
# columns(hgu133a.db)
# keytypes(hgu133a.db)
# head(keys(hgu133a.db, keytype="PROBEID"))
anno <- AnnotationDbi::select(hgu133a.db,
                               keys = (featureNames(manfiltered)),
                               columns = c("SYMBOL", "GENENAME"),
                               keytype = "PROBEID")

anno <- subset(anno, !is.na(SYMBOL))


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

#### Linear models ####
#. limma: Linear models for microarrays ----
# make a subset
individual <- pD$id
timepoint <- pD$timepoint
gender <- pD$gender

# make a group by individual and disease type
male <- individual[gender == "male"]
t_male <- timepoint[gender == "male"]
design_male <- model.matrix(~ 0 + t_male + male)
colnames(design_male)[1:2] <- unique(timepoint[gender == "male"])
rownames(design_male) <- male

female <- individual[gender == "female"]
t_female <- timepoint[gender == "female"]
design_female <- model.matrix(~ 0 + t_female + female)
colnames(design_female)[1:2] <- unique(timepoint[gender == "female"])
rownames(design_female) <- female

# inspect the design matrices
head(design_male)
head(design_female)


#. Differential expression analysis of the CRAT gene. ----
# Contrasts and hypotheses tests: all genes
contrast_matrix_male <- makeContrasts(d1-d0, levels = design_male)
fit_male <- eBayes(
                contrasts.fit(
                    lmFit(
                        final[,gender == "male"],
                        design = design_male),
                        contrast_matrix_male))

contrast_matrix_female <- makeContrasts(d1-d0, levels = design_female)
fit_female <- eBayes(
                contrasts.fit(
                    lmFit(
                        final[,gender == "female"],
                        design = design_female),
                    contrast_matrix_female))


#### Extracting results: topTable() ####
table_male <- topTable(fit_male, number = Inf)
head(table_male)
hist(table_male$P.Value, col = brewer.pal(3, name = "Set2")[1], xlab = "male p-values")
hist(table_male$adj.P.Val, col = brewer.pal(3, name = "Set2")[2], xlab = "male adj.p-values")

table_female <- topTable(fit_female, number = Inf)
head(table_female)
hist(table_female$P.Value, col = brewer.pal(3, name = "Set2")[1], xlab = "female p-values")
hist(table_female$adj.P.Val, col = brewer.pal(3, name = "Set2")[2], xlab = "female p-values")


#. Multiple testing FDR, and comparison with results from the original paper
nrow(subset(table_male, P.Value < 0.05))  # 853
tail(subset(table_male, P.Value < 0.05))
# original: Using the same filtering criteria for the male skeletal muscle samples,
#           we found significantly fewer expression changes
#           (182 genes with paired t-test <0.05 in men compared with 492 in women).


nrow(subset(table_female, P.Value < 0.05))  # 1308
tail(subset(table_female, P.Value < 0.05))
# original: 429 genes were differentially expressed in skeletal muscle
#           of women after 9 mo of aerobic exercise training on the basis of P value
#           with no adjustment for multiple testing.


#### Write out the result ####
write.csv(table_male, file = "res_GSE1295M.csv")
write.csv(table_female, file = "res_GSE1295F.csv")


# Exploring Results
head(table_male)
nrow(subset(table_male, P.Value < 0.05))  # 853
nrow(subset(table_male, adj.P.Val < 0.01))  # 0

nrow(subset(table_female, P.Value < 0.05))  # 1308
nrow(subset(table_female, adj.P.Val < 0.01))  # 0

resSig_male <- subset(table_male, P.Value < 0.05)
head(resSig_male[order(resSig_male$P.Value), ])

resSig_female <- subset(table_female, P.Value < 0.05)
head(resSig_female[order(resSig_female$P.Value), ])

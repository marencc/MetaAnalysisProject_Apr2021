#### Endurance: 10 male, 12wks 60-80min cycling/D ####
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116801
# GSE116801

#### List of packages required for the workflow ####
# General Bioconductor packages
library(Biobase)
library(oligoClasses)

# For importing data
library(ArrayExpress)

# Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

# Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)

# Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)

# Formatting/documentation packages
library(rmarkdown)
library(BiocStyle)
library(dplyr)
library(tidyr)

# Helpers:
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
library(devtools)


#### Download the raw data from GEO ####
library(GEOquery)
getGEOSuppFiles("GSE116801", makeDirectory = F, baseDir = raw_data_dir)
RAW.tar <- file.path(raw_data_dir, "GSE116801_RAW.tar")
CEL.dir <- paste(raw_data_dir, "CEL", sep="/")
untar(RAW.tar, exdir = CEL.dir)
head(list.files(CEL.dir, "CEL"))

celFiles <- list.files(CEL.dir, pattern = "CEL.gz$", full = TRUE)
sapply(celFiles, gunzip, overwrite = TRUE)
head(list.files(CEL.dir, "CEL"))

#. For PhenoData ----
geoMat <- getGEO("GSE116801")
pData <- pData(geoMat[[1]]) # all
head(pData)
dim(pData)

colnames(raw_data) <- gsub("_Jos.*", "", colnames(raw_data))
all(rownames(pData)==colnames(exprs(raw_data)))

metadata <- data.frame(labelDescription=colnames(pData),
                       row.names=colnames(pData))

pData <- AnnotatedDataFrame(data=pData, varMetadata=metadata)


#### Create Expression Set ####
raw_data <- oligo::read.celfiles(filenames = file.path("CEL",list.files( "CEL")),
                                  verbose = FALSE, phenoData = pData)
stopifnot(validObject(raw_data))

# Take a look to make a subset
# make sure if rownames(row_data) corresponds to geo_accession
raw_data
varLabels(raw_data)

#. Annotation database based on the chip used. ----
library(pd.hg.u133.plus.2)
library(hgu133plus2.db)

#. Subselect the columns of interest ----
head(Biobase::pData(raw_data))
tail(Biobase::pData(raw_data))
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[, c("title", # individual
                                                         "title", # timepoint
                                                         "geo_accession", # just in case
                                                         "characteristics_ch1")] # to ensure if all is male

#. Rename ----
names(pData(raw_data)) <- c("individual", "timepoint", "geo", "gender")
pData(raw_data)$individual <- gsub("Exercise_Subject|Sedentary_Subject","",pData(raw_data)$individual)
pData(raw_data)$individual <- (paste0("S", pData(raw_data)$individual))
pData(raw_data)$timepoint <- str_sub(pData(raw_data)$timepoint, 1, 2)
pData(raw_data)$id <- paste(pData(raw_data)$individual, pData(raw_data)$timepoint, sep=".")
head(pData(raw_data))
tail(pData(raw_data))
dim(pData(raw_data))

#### Quality control of the raw data ####
# Take a look
Biobase::exprs(raw_data)[1:5, 1:5]

# log2
exp_raw <- log2(Biobase::exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

#. PCA plot of the log ----
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     individual = pData(raw_data)$individual,
                     timepoint = pData(raw_data)$timepoint,
                     id = pData(raw_data)$id
)

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(shape = timepoint, colour = individual)) +
    ggtitle("PCA timepoint(S) & individual(C)") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5))+
    coord_fixed(ratio = sd_ratio) +
    theme(legend.position = "none")
# outlierが3つくらいありそう


#. Intensity boxplots of the log2 ----
oligo::boxplot(raw_data, target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               las = 2,
               cex.axis=0.7)

# GSM3262055, 059, 063 (多分上記ggplotで思った3つと同じ)
attention <- rownames(pData(raw_data)) %in% c("GSM3262055", "GSM3262059", "GSM3262063")
attention <- pData(raw_data)[attention, ]
attention


#. more elaborate quality control plots ----
arrayQualityMetrics(expressionset = raw_data,
                    outdir = "QCreport",
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = "timepoint")
dev.off()

# arrayQualityMetricsによるoutlier
pData(raw_data)[c(9,13,17),]

#### Background adjustment, calibration, summarization and annotation ####
head(ls("package:hgu133plus2.db"))

#. oligo::rma() ----
# Do the Bcground adjstmt, calbratn, sumrztn
eset <- oligo::rma(raw_data, normalize = FALSE)

#. RLE: Relative Log Expression data quality analysis ----
row_medians_assayData <- 
    Biobase::rowMedians(as.matrix(Biobase::exprs(eset)))

RLE_data <- sweep(Biobase::exprs(eset), 1, row_medians_assayData)
RLE_data <- as.data.frame(RLE_data)

RLE_data_gathered <-
    tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

# Warning message:
# Removed 22521 rows containing non-finite values (stat_boxplot).
ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
    geom_boxplot(outlier.shape = NA) + 
    ylim(c(-2, 2)) + 
    theme(axis.text.x = element_text(colour = "aquamarine4",
                                     angle = 60, size = 6.5, hjust = 1 ,
                                     face = "bold"))
# 055, 059, 063 outlier (同上)
pData(raw_data)[c(9,13,17),]


#### Auto-normalize with RMA ####
# normalize = FALSE を入れないver (targetをベースにnormalizeされる)
eset_norm <- oligo::rma(raw_data)

#. PCA analysis ----
myexp <- Biobase::exprs(eset_norm)
PCA <- prcomp(t(myexp), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     individual = pData(raw_data)$individual,
                     timepoint = pData(raw_data)$timepoint,
                     id = pData(raw_data)$id
)

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(shape = timepoint, colour = individual)) +
    ggtitle("PCA_Norm timepoint(S) & individual(C)") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5))+
    coord_fixed(ratio = sd_ratio) +
    theme(legend.position = "none")


#. Heatmap clustering analysis ----
i.sensitivity_names <- ifelse(str_detect(pData(eset_norm)$i.sensitivity,
                                         "non"), "non_infl.", "infl.")

timepoint_names <- ifelse(str_detect(pData(eset_norm)$timepoint,
                                     "Crohn"), "CD", "UC")

annotation_for_heatmap <-
    data.frame(i.sensitivity = i.sensitivity_names, timepoint = timepoint_names)

row.names(annotation_for_heatmap) <- row.names(pData(eset_norm))

# change the default distance calculation method(Euclidean) to "manhattan"
dists <- as.matrix(dist(t(myexp), method = "manhattan"))

rownames(dists) <- row.names(pData(eset_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

ann_colors <- list(
    i.sensitivity = c(non_infl. = "chartreuse4", infl. = "burlywood3"),
    timepoint = c(CD = "blue4", UC = "cadetblue2")
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
man_threshold <- round(mean(medians))

hist_res <- hist(medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)

#. Exclude genes that below the threshold ----
# Transcripts that do not have intensities larger than the threshold in
# at least as many arrays as the smallest experimental group are excluded.
# get a list:
no_of_samples <- table(pData(eset_norm)$timepoint)
no_of_samples
no_of_samples <- no_of_samples*0.8  # 10だと多すぎるかと思ったので8割にした

samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(Biobase::exprs(eset_norm), 1,
                           function(x){
                               sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)


manfiltered <- subset(eset_norm, idx_man_threshold)

# 変化確認
head(exprs(eset_norm))
head(exprs(manfiltered))

dim(exprs(eset_norm))
dim(exprs(manfiltered))


#### Annotation of the transcript clusters ####
anno <- AnnotationDbi::select(hgu133plus2.db,
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
#. . design
individual <- factor(Biobase::pData(final)$individual)
timepoint <- factor(Biobase::pData(final)$timepoint)

# make a group by individual and timepoint type
design <- model.matrix(~ 0 + timepoint + individual)
colnames(design) <- c(levels(timepoint),levels(individual)[-1])
rownames(design) <- individual

# inspect the design matrices
head(design)


#. . Contrasts matrix ####
contrastMatrix <- makeContrasts(Ex-Se, levels = design)

fit <- eBayes(contrasts.fit(lmFit(final, design = design),
                            contrastMatrix))


#### Extracting results: topTable() ####
#. Multiple testing FDR
table <- topTable(fit, number = Inf)
head(table)
nrow(subset(table, P.Value < 0.001)) # 17
nrow(subset(table, adj.P.Val < 0.05)) # 0
tail(subset(table, P.Value < 0.001))
tail(subset(table, adj.P.Val < 0.05))

write.csv(table, "Res_DEG_GSE116801.csv")


hist(table$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "inflamed vs non-inflamed - Crohn’s timepoint", xlab = "p-values")
hist(table$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
     main = "inflamed vs non-inflamed - Crohn’s timepoint", xlab = "adj.p-values")


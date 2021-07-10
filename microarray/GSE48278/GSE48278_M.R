#### GSE48278: Various exercise groups, both sex ####
# A: Mild aerobic exercise - low amount/mod intensity n=10
# B: Inactive control n=10
# C: High aerobic exercise - high amt/vig intensity n=10
# D: Moderate aerobic exercise - low amt/vig intensity n=10
# E: Mod aerobic + resistance training n=8
# F: Resistance training only n=9

# https://link.springer.com/article/10.1007/s00125-014-3343-4
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48278
# Affymetrix Human Genome U133 Plus 2.0 Array
# E-GEOD-48278
# GSE48278

# library("maEndToEnd")

#### List of packages required for the workflow ####
# Some Helper/Styling packages have been commented here,
# as they are not strictly neccesary to execute the workflow.
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


#### Download the raw data from ArrayExpress ####
raw_data_dir <- tempdir()

if (!dir.exists(raw_data_dir)) {
    dir.create(raw_data_dir)
}

anno_AE <- getAE("E-GEOD-48278", path = raw_data_dir, type = "raw")


#### Import of annotation data and microarray expression data as “ExpressionSet” ####
sdrf_location <- file.path(raw_data_dir, "E-GEOD-48278.sdrf.txt")
SDRF <- read.delim(sdrf_location)

rownames(SDRF) <- SDRF$Array.Data.File # sample names is stored here
SDRF <- AnnotatedDataFrame(SDRF) # to create an ExpressionSet later on


#### Create Expression Set ####
# Specified SDRF file we created earlier as phenoData.
# Thus, we had to make sure to import the CEL files in the order corresponds to the SDRF table.
# To enforce this, we used the column Array.Data.File of the SDRF table as the filenames argument.
raw_data <- oligo::read.celfiles(filenames = file.path(raw_data_dir,
                                                       SDRF$Array.Data.File),
                                 verbose = FALSE, phenoData = SDRF)
stopifnot(validObject(raw_data))

# Take a look to make a subset
raw_data
varLabels(raw_data)

#. Annotation database based on the chip used. ----
library(pd.hg.u133.plus.2)
library(hgu133plus2.db)

#. Subselect the columns of interest ----
head(Biobase::pData(raw_data))
tail(Biobase::pData(raw_data))
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[, c("Characteristics..subject.id.", # individual
                                                         "Characteristics..Sex.",        # gender
                                                         "Comment..Sample_title.",       # timepoint
                                                         "Characteristics..exercise.group.")] # group

#. Rename ----
names(pData(raw_data)) <- c("individual", "gender", "timepoint", "group")
pData(raw_data)$individual <- (paste0("S", pData(raw_data)$individual))
pData(raw_data)$timepoint <- ifelse(str_detect(pData(raw_data)$timepoint, "POST"), "POST", "PRE")
pData(raw_data)$grp.short <- str_sub(Biobase::pData(raw_data)$group, 1, 1)
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
                     gender = pData(raw_data)$gender,
                     timepoint = pData(raw_data)$timepoint,
                     group = pData(raw_data)$grp.short
)

# all gender × group
ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(gender ~ group)
# female × group C may have sth?


#. Create subsets ----
# female × exercised
male <- raw_data[,raw_data$gender == "male"]
keep <- c("A","C","D","E","F")
filter <- colnames(male)[male$grp.short %in% keep]
length(filter)

male <- male[,filter]
male
pData(male)

#. Intensity boxplots of the log2 ----
oligo::boxplot(male, target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               las = 2,
               cex.axis=0.7)

#. more elaborate quality control plots ----
arrayQualityMetrics(expressionset = male,
                    outdir = "QCreport",
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = c("timepoint", "grp.short", "gender"))
dev.off()

# Outlier
pData(male)[c(7,8,15,16,19,20,41,42,43,44),]

# exclude outliers
outlier_samples <- c(7,8,15,16,19,20,41,42,43,44)
male <- male[,-outlier_samples]


#### Background adjustment, calibration, summarization and annotation ####
head(ls("package:hgu133plus2.db"))

#. oligo::rma() ----
# Do the Bcground adjstmt, calbratn, sumrztn
eset <- oligo::rma(male, normalize = FALSE)

#. RLE: Relative Log Expression data quality analysis ----
row_medians_assayData <- 
    Biobase::rowMedians(as.matrix(Biobase::exprs(eset)))

RLE_data <- sweep(Biobase::exprs(eset), 1, row_medians_assayData)
RLE_data <- as.data.frame(RLE_data)

RLE_data_gathered <-
    tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

# Warning message:
# Removed 4246 rows containing non-finite values (stat_boxplot).
ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
    geom_boxplot(outlier.shape = NA) + 
    ylim(c(-2, 2)) + 
    theme(axis.text.x = element_text(colour = "aquamarine4",
                                     angle = 60, size = 6.5, hjust = 1 ,
                                     face = "bold"))


#### Auto-normalize with RMA ####
# normalize = FALSE を入れないver (targetをベースにnormalizeされる)
eset_norm <- oligo::rma(male)

#. PCA analysis ----
myexp <- Biobase::exprs(eset_norm)
PCA <- prcomp(t(myexp), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     individual = pData(male)$individual,
                     gender = pData(male)$gender,
                     timepoint = pData(male)$timepoint,
                     group = pData(male)$grp.short
)

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid( ~ group)

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint))

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
man_threshold <- 3.7 # why 4?

hist_res <- hist(medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)

#. Exclude genes that below the threshold ----
# Transcripts that do not have intensities larger than the threshold in
# at least as many arrays as the smallest experimental group are excluded.
# get a list:
no_of_samples <-
    table(paste0(pData(eset_norm)$grp.short, "_",
                 pData(eset_norm)$timepoint))
no_of_samples

samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(Biobase::exprs(eset_norm), 1,
                           function(x){
                               sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)

manfiltered <- subset(eset_norm, idx_man_threshold)

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
# make a subset
individual <- factor(Biobase::pData(final)$individual)
timepoint <- factor(Biobase::pData(final)$timepoint)
group <- factor(Biobase::pData(final)$grp.short)

# make a group by individual and timepoint type
i_A <- individual[group == "A"]
design_A <- model.matrix(~ 0 + timepoint[group == "A"] + i_A)
colnames(design_A) <- c(levels(timepoint),levels(i_A)[-1])
rownames(design_A) <- i_A
head(design_A)

i_C <- individual[group == "C"]
design_C <- model.matrix(~ 0 + timepoint[group == "C"] + i_C)
colnames(design_C) <- c(levels(timepoint),levels(i_C)[-1])
rownames(design_C) <- i_C
head(design_C)

i_D <- individual[group == "D"]
design_D <- model.matrix(~ 0 + timepoint[group == "D"] + i_D)
colnames(design_D) <- c(levels(timepoint),levels(i_D)[-1])
rownames(design_D) <- i_D
head(design_D)

i_E <- individual[group == "E"]
design_E <- model.matrix(~ 0 + timepoint[group == "E"] + i_E)
colnames(design_E) <- c(levels(timepoint),levels(i_E)[-1])
rownames(design_E) <- i_E
head(design_E)

i_F <- individual[group == "F"]
design_F <- model.matrix(~ 0 + timepoint[group == "F"] + i_F)
colnames(design_F) <- c(levels(timepoint),levels(i_F)[-1])
rownames(design_F) <- i_F
head(design_F)

# inspect the design matrices
head(design_A)
head(design_C)
head(design_D)
head(design_E)
head(design_F)


#### Contrasts and hypotheses tests: all genes ####
contrast_matrix_A <- makeContrasts(POST-PRE, levels = design_A)
fit_A <- eBayes(contrasts.fit(lmFit(final[, group == "A"],
                                    design = design_A),
                              contrast_matrix_A))

contrast_matrix_C <- makeContrasts(POST-PRE, levels = design_C)
fit_C <- eBayes(contrasts.fit(lmFit(final[, group == "C"],
                                    design = design_C),
                              contrast_matrix_C))

contrast_matrix_D <- makeContrasts(POST-PRE, levels = design_D)
fit_D <- eBayes(contrasts.fit(lmFit(final[, group == "D"],
                                    design = design_D),
                              contrast_matrix_D))

contrast_matrix_E <- makeContrasts(POST-PRE, levels = design_E)
fit_E <- eBayes(contrasts.fit(lmFit(final[, group == "E"],
                                    design = design_E),
                              contrast_matrix_E))

contrast_matrix_F <- makeContrasts(POST-PRE, levels = design_F)
fit_F <- eBayes(contrasts.fit(lmFit(final[, group == "F"],
                                    design = design_F),
                              contrast_matrix_F))


#### Extracting results: topTable() ####
#. Multiple testing FDR, and comparison with results from the original paper
# male
table_A <- topTable(fit_A, number = Inf)
head(table_A)
nrow(subset(table_A, P.Value < 0.001))
nrow(subset(table_A, adj.P.Val < 0.05))
tail(subset(table_A, P.Value < 0.001))
tail(subset(table_A, adj.P.Val < 0.05))

write.csv(table_A, "Res_DEG_GSE48278_MA.csv")


table_C <- topTable(fit_C, number = Inf)
head(table_C)
nrow(subset(table_C, P.Value < 0.001))
nrow(subset(table_C, adj.P.Val < 0.05))
tail(subset(table_C, P.Value < 0.001))
tail(subset(table_C, adj.P.Val < 0.05))

write.csv(table_C, "Res_DEG_GSE48278_MC.csv")


table_D <- topTable(fit_D, number = Inf)
head(table_D)
nrow(subset(table_D, P.Value < 0.001))
nrow(subset(table_D, adj.P.Val < 0.05))
tail(subset(table_D, P.Value < 0.001))
tail(subset(table_D, adj.P.Val < 0.05))

write.csv(table_D, "Res_DEG_GSE48278_MD.csv")


table_E <- topTable(fit_E, number = Inf)
head(table_E)
nrow(subset(table_E, P.Value < 0.001))
nrow(subset(table_E, adj.P.Val < 0.05))
tail(subset(table_E, P.Value < 0.001))
tail(subset(table_E, adj.P.Val < 0.05))

write.csv(table_E, "Res_DEG_GSE48278_ME.csv")


table_F <- topTable(fit_F, number = Inf)
head(table_F)
nrow(subset(table_F, P.Value < 0.001))
nrow(subset(table_F, adj.P.Val < 0.05))
tail(subset(table_F, P.Value < 0.001))
tail(subset(table_F, adj.P.Val < 0.05))

write.csv(table_F, "Res_DEG_GSE48278_MF.csv")

### 2021.05.13 ここまで


hist(table_A$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "inflamed vs non-inflamed - Crohn’s timepoint", xlab = "p-values")
hist(table_A$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
     main = "inflamed vs non-inflamed - Crohn’s timepoint", xlab = "adj.p-values")


#### Visualization of DE analysis results ####
#. Volcano plot ----
volcano_names <- ifelse(abs(fit_M$coefficients)>=1,
                        fit_M$genes$SYMBOL, NA)

volcanoplot(fit_M, coef = 1L, style = "p-value", highlight = 100,
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)


#### Gene ontology based enrichment analysis ####
# choose FDR cutoff at 10%
DE_genes_M <- subset(table_M, adj.P.Val < 0.1)$PROBEID

# Matching the background set of genes
back_genes_idx <- genefilter::genefinder(final,
                                         as.character(DE_genes_M),
                                         method = "manhattan", scale = "none")

# extract the PROBEIDs: correspond to the indices
# matrix with the DE-genes as column names and the PROBEIDs of the corresponding background genes in the cells
back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)

# create a vector back_genes containing all background gene PROBEIDs
back_genes <- featureNames(final)[back_genes_idx]
back_genes <- setdiff(back_genes, DE_genes_M)
intersect(back_genes, DE_genes_M)
length(back_genes)

# Multidensity plot: x=mean expression
multidensity(list(
    all = table_M[,"AveExpr"] ,
    fore = table_M[DE_genes_M , "AveExpr"],
    back = table_M[rownames(table_M) %in% back_genes, "AveExpr"]),
    col = c("#e46981", "#ae7ee2", "#a7ad4a"),
    xlab = "mean expression",
    main = "DE genes for CD-background-matching")

#. Running topGO ----
# GO has 3 types: Cellular component (CC), biological processes (BP), and molecular function (MF)
# for illustrative purposes we limit ourselves to the BP category here
# create a named vector all_genes with all genes to be analyzed, i.e. DE-genes and background genes:
gene_IDs <- rownames(table_M)
in_universe <- gene_IDs %in% c(DE_genes_M, back_genes)
in_selection <- gene_IDs %in% DE_genes_M

all_genes <- in_selection[in_universe]
all_genes <- factor(as.integer(in_selection[in_universe]))
names(all_genes) <- gene_IDs[in_universe]

top_GO_data <- new("topGOdata", ontology = "BP", allGenes = all_genes,
                   nodeSize = 10, annot = annFUN.db, affyLib = "hgu133plus2.db")

result_top_GO_elim <-
    runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")
result_top_GO_classic <-
    runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")

res_top_GO <- GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                       Fisher.classic = result_top_GO_classic,
                       orderBy = "Fisher.elim" , topNodes = 100)

genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                           chip = "hgu133plus2.db", geneCutOff = 1000)

res_top_GO$sig_genes <- sapply(genes_top_GO, function(x){
    str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"],";"),
          collapse = "")
})

head(res_top_GO[,1:8], 20)

#. Visualization of the GO-analysis results ----
showSigOfNodes(top_GO_data, score(result_top_GO_elim), firstSigNodes = 3,
               useInfo = 'def')

#. pathway: reactome ----
entrez_ids <- mapIds(hgu133plus2.db,
                     keys = rownames(table_M),
                     keytype = "PROBEID",
                     column = "ENTREZID")

reactome_enrich <- enrichPathway(gene = entrez_ids[DE_genes_M],
                                 universe = entrez_ids[c(DE_genes_M,
                                                         back_genes)],
                                 organism = "human",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.9,
                                 readable = TRUE)

reactome_enrich@result$Description <- paste0(str_sub(
    reactome_enrich@result$Description, 1, 20),
    "...")

head(summary(reactome_enrich))[1:6]

#. barplot() ----
barplot(reactome_enrich)

#. emapplot() ----
emapplot(reactome_enrich, showCategory = 10)



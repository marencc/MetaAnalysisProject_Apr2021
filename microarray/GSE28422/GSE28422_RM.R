#### GSE28422: Resistance ####
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28422
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3365403/
# 28 young/old adults
# yng: n=16, 24±1y 8m & 8f
# old: n=12, 84±1y 6m & 6f
# male: n=14, male: n=14
# 12 wks of progressive resistance training (3D/wk = 36 sessions)
# biopsies (4): pre, post of 1st & 36th training
# E-GEOD-28422
# GSE28422

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

anno_AE <- getAE("E-GEOD-28422", path = raw_data_dir, type = "raw")


#### Import of annotation data and microarray expression data as “ExpressionSet” ####
sdrf_location <- file.path(raw_data_dir, "E-GEOD-28422.sdrf.txt")
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

# Take a look the data
raw_data
varLabels(raw_data)

#. Annotation database based on the chip used. ----
library(pd.hg.u133.plus.2)
library(hgu133plus2.db)

#. Subselect the columns of interest ----
head(Biobase::pData(raw_data))
tail(Biobase::pData(raw_data))
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[, c("Comment..Sample_description.", # individual
                                                         "Characteristics.gender.")] # gender

#. Rename ----
names(pData(raw_data)) <- c("individual", "gender")
pData(raw_data)$individual <- gsub("EA08017_.*_AGEX_OF_|EA08017_.*_AGEX_YF_|EA08017_.*_AGEX_OM_|EA08017_.*_AGEX_YM_","",pData(raw_data)$individual)
pData(raw_data)$individual <- str_sub(Biobase::pData(raw_data)$individual, end = -11)
pData(raw_data)$id <- paste0("S", str_sub(Biobase::pData(raw_data)$individual, 4, 5)) # individual (ID)

pData(raw_data)$timepoint <- str_sub(Biobase::pData(raw_data)$individual, 1, 2)

head(pData(raw_data))
tail(pData(raw_data))
dim(pData(raw_data))

#. Create subsets ----
# male
male <- raw_data[,raw_data$gender == "male"]
keep <- c("T1", "T3")
filter <- colnames(male)[male$timepoint %in% keep]
length(filter)

male <- male[,filter]
male
pData(male)

#### Quality control of the raw data ####
# Take a look
Biobase::exprs(male)[1:5, 1:5]

# log2
exp_raw <- log2(Biobase::exprs(male))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

#. PCA plot of the log ----
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     individual = pData(male)$id,
                     timepoint_id = pData(male)$individual,
                     timepoint = pData(male)$timepoint
)

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(shape = timepoint, colour = individual)) +
    ggtitle("PCA timepoint(S) & individual(C)") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5))+
    coord_fixed(ratio = sd_ratio)

#. Intensity boxplots of the log2 ----
oligo::boxplot(male, target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               las = 2,
               cex.axis=0.7)

#. more elaborate quality control plots ----
arrayQualityMetrics(expressionset = male,
                    outdir = "QCreport",
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = c("timepoint", "id"))
dev.off()

# outlierっぽい
pData(male)[c(22),]

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
                     individual = pData(male)$id,
                     id.time = pData(male)$individual,
                     timepoint = pData(male)$timepoint
)


ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(shape = timepoint, colour = individual)) +
    ggtitle("PCA_Norm timepoint(S) & individual(C)") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5))+
    coord_fixed(ratio = sd_ratio)

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
man_threshold <- 3.4 # how to decide?

hist_res <- hist(medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 3.4)

#. Exclude genes that below the threshold ----
# Transcripts that do not have intensities larger than the threshold in
# at least as many arrays as the smallest experimental group are excluded.
# get a list:
no_of_samples <-
    table(pData(eset_norm)$timepoint)
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
individual <- factor(Biobase::pData(final)$id)
id.time <- factor(Biobase::pData(final)$individual)
timepoint <- factor(Biobase::pData(final)$timepoint)
table(pData(final))

# make a design
design <- model.matrix(~ 0 + timepoint + individual)
colnames(design) <- c(levels(timepoint),levels(individual)[-1])
rownames(design) <- individual
table(design)

# inspect the design matrices
head(design)

#### Contrasts and hypotheses tests: all genes ####
contrast_matrix <- makeContrasts(T1-T3, levels = design)

fit_T1T3 <- eBayes(contrasts.fit(lmFit(final,
                                       design = design),
                                 contrast_matrix))

#### Extracting results: topTable() ####
#. Multiple testing FDR, and comparison with results from the original paper
# male
table_T1T3 <- topTable(fit_T1T3, number = Inf)
head(table_T1T3)
nrow(subset(table_T1T3, P.Value < 0.001))
nrow(subset(table_T1T3, adj.P.Val < 0.05))
tail(subset(table_T1T3, P.Value < 0.001))
tail(subset(table_T1T3, adj.P.Val < 0.05))

write.csv(table_T1T3, "Res_GSE28422_RM_T1T3.csv")


hist(table_T1T3$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "inflamed vs non-inflamed - Crohn’s timepoint", xlab = "p-values")
hist(table_T1T3$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
     main = "inflamed vs non-inflamed - Crohn’s timepoint", xlab = "adj.p-values")

#### Gene ontology based enrichment analysis ####
# choose FDR cutoff at 10%
DE_genes <- subset(table_T1T3, adj.P.Val < 0.1)$PROBEID

# Matching the background set of genes
back_genes_idx <- genefilter::genefinder(final,
                                         as.character(DE_genes),
                                         method = "manhattan", scale = "none")

# extract the PROBEIDs: correspond to the indices
# matrix with the DE-genes as column names and the PROBEIDs of the corresponding background genes in the cells
back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)

# create a vector back_genes containing all background gene PROBEIDs
back_genes <- featureNames(final)[back_genes_idx]
back_genes <- setdiff(back_genes, DE_genes)
intersect(back_genes, DE_genes)
length(back_genes)

# Multidensity plot: x=mean expression
multidensity(list(
    all = table_T1T3[,"AveExpr"] ,
    fore = table_T1T3[DE_genes , "AveExpr"],
    back = table_T1T3[rownames(table_T1T3) %in% back_genes, "AveExpr"]),
    col = c("#e46981", "#ae7ee2", "#a7ad4a"),
    xlab = "mean expression",
    main = "DE genes for CD-background-matching")

#. Running topGO ----
# GO has 3 types: Cellular component (CC), biological processes (BP), and molecular function (MF)
# for illustrative purposes we limit ourselves to the BP category here
# create a named vector all_genes with all genes to be analyzed, i.e. DE-genes and background genes:
gene_IDs <- rownames(table_T1T3)
in_universe <- gene_IDs %in% c(DE_genes, back_genes)
in_selection <- gene_IDs %in% DE_genes

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
                     keys = rownames(table_T1T3),
                     keytype = "PROBEID",
                     column = "ENTREZID")

reactome_enrich <- enrichPathway(gene = entrez_ids[DE_genes],
                                 universe = entrez_ids[c(DE_genes,
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
# ErroError in has_pairsim(x) : 
# Term similarity matrix not available.
# Please use pairwise_termsim function to deal with the results of enrichment analysis.

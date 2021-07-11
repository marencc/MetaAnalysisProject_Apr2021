#### GSE48278 ####
#### Loading packages ####
library(Biobase)
library(oligoClasses)
library(oligo)
# library(arrayQualityMetrics)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(devtools)


#### Loading the data from GEO ####
library(GEOquery)
geo <- getGEO("GSE48278")
length(geo) # 1
gse <- geo[[1]]
gse
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 54675 features, 114 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM1174042 GSM1174043 ... GSM1174155 (114 total)
# varLabels: title geo_accession ... tissue:ch1 (43 total)
# varMetadata: labelDescription
# featureData
# featureNames: 1007_s_at 1053_at ... AFFX-TrpnX-M_at (54675 total)
# fvarLabels: ID GB_ACC ... Gene Ontology Molecular Function (16 total)
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 25091629 
# Annotation: GPL570


exprs(gse)[1:3, 1:3]
#           GSM1174042 GSM1174043 GSM1174044
# 1007_s_at   681.0750    770.095   613.2700
# 1053_at     138.6410    133.144   114.8640
# 117_at       22.9964     30.508    22.1818



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

pD <- pD.all[, c("title",  # individual & timepoint
                 "Sex:ch1",  # gender
                 "exercise group:ch1"  # group
                 )] 
names(pD)[2:3] <- c("gender", "ex.group")
head(pD, 3)
unique(pD$title)
unique(pD$`exercise group:ch1`)
# [1] "A (Mild aerobic exercise)"         "B (Inactive control)"              "C (High aerobic exercise)"        
# [4] "D (Moderate aerobic exercise)"     "E (Aerobic + Resistance exercise)" "F (Resistance exercise)" 


#. Renaming ----
pD$id <- str_sub(pD$title, 1, 4)  # length(unique(pD$id))  # 57
pD$ex.group.short <- str_sub(pD$ex.group, 1, 1)
pD$group <- paste0(pD$gender, pD$ex.group.short)
pD$timepoint <- factor(ifelse(str_detect(pD$title, "PRE"), "pre", "post"))
pD$timepoint <- relevel(pD$timepoint, ref = "pre")
levels(pD$timepoint)
# [1] "pre"  "post"


#. excluding group B (Inactive control) ----
keep <- grep("B", pD$ex.group.short, invert = TRUE)
pD <- pD[keep,]
gse <- gse[,keep]
dim(pD)  # 94  7
dim(gse)
# Features  Samples 
# 54675       94

head(pD)
# title gender                  ex.group   id ex.group.short timepoint group
# GSM1174042 S137_A_POST      F A (Mild aerobic exercise) S137              A      post    FA
# GSM1174043  S137_A_PRE      F A (Mild aerobic exercise) S137              A       pre    FA
# GSM1174044 S139_A_POST      M A (Mild aerobic exercise) S139              A      post    MA
# GSM1174045  S139_A_PRE      M A (Mild aerobic exercise) S139              A       pre    MA
# GSM1174050 S146_C_POST      M C (High aerobic exercise) S146              C      post    MC
# GSM1174051  S146_C_PRE      M C (High aerobic exercise) S146              C       pre    MC
tail(pD)
# title gender                          ex.group   id ex.group.short timepoint group
# GSM1174150 S379_E_POST      F E (Aerobic + Resistance exercise) S379              E      post    FE
# GSM1174151  S379_E_Pre      F E (Aerobic + Resistance exercise) S379              E      post    FE
# GSM1174152 S387_E_POST      M E (Aerobic + Resistance exercise) S387              E      post    ME
# GSM1174153  S387_E_PRE      M E (Aerobic + Resistance exercise) S387              E       pre    ME
# GSM1174154 S401_F_POST      F           F (Resistance exercise) S401              F      post    FF
# GSM1174155  S401_F_PRE      F           F (Resistance exercise) S401              F       pre    FF




#### Normalization ####
head(pData(gse)$data_processing, 1)
# [1] "CHP files were produced in Expression Console using the PLIER algorithm.
# CHP files were imported into Partek Genomics Suite for statistical analysis.
# RM ANCOVA was used to test group and time factors, using sex and age as covariates.
# Regression analyses were used in genotype/gene expression or genotype/phenotype association testing."


#### checking data intensities
for (i in unique(pD$group)) {
    check <- gse[pD$group == i,]
    oligo::boxplot(log2(exprs(check)), target = "core",
                   main = paste(i, "Boxplot of log2-intensitites for the raw data"),
                   outline = FALSE,
                   las = 2,
                   cex.axis=0.7)    
    }
oligo::boxplot(log2(exprs(gse)), target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               outline = FALSE,
               las = 2,
               cex.axis=0.7)
# memo: based on the boxplot, the data doesn't seem to be normalized.

library(affyPLM)
norm <- normalize.ExpressionSet.quantiles(gse)
exprs(norm)[1:5, 1:5]
oligo::boxplot(log2(exprs(norm)), target = "core",
               main = "Boxplot of log2-intensitites for the raw data",
               outline = FALSE,
               las = 2,
               cex.axis=0.7)



#### Quality control of the raw data ####
# Take a look
exprs(gse)[1:5, 1:5]
#           GSM1174042 GSM1174043 GSM1174044 GSM1174045 GSM1174050
# 1007_s_at   681.0750   770.0950   613.2700   790.6360  317.61700
# 1053_at     138.6410   133.1440   114.8640   108.7390    7.12879
# 117_at       22.9964    30.5080    22.1818    47.2096   10.08110
# 121_at      300.0870   211.5810   194.9230   273.4350  129.47100
# 1255_g_at    61.3044    73.7449    86.1957    53.0646   49.00420


exprs(norm)[1:5, 1:5]
# GSM1174042 GSM1174043 GSM1174044 GSM1174045 GSM1174050
# 1007_s_at  609.49438  757.40280  625.78678  682.21046  589.19390
# 1053_at    126.79606  120.76283  113.29286   95.04914   12.62705
# 117_at      21.66635   28.80891   21.76640   43.73073   16.04390
# 121_at     272.18590  195.15600  195.21189  232.72597  226.12387
# 1255_g_at   56.81878   66.84979   85.48805   49.04166   75.02016


# log2
#. PCA plot of the log ----
# head(pD, 3)
data_types <- c(gse, norm)
for (i in 1:length(data_types)) {
    if (i == 1) { print("gse")
        data <- gse
        name <- "gse"
    } else { print("norm")
        data <- norm
        name <- "norm"
        }
    exp_raw <- log2(exprs(data))
    PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
    percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
    sd_ratio <- sqrt(percentVar[2] / percentVar[1])
    
    dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                         individual = pD$id,
                         timepoint = pD$timepoint,
                         group = pD$group,
                         gender = pD$gender
                         )
    
    print(ggplot(dataGG, aes(PC1, PC2)) +
        geom_point(aes(colour = timepoint)) +
        facet_wrap( ~ group, ncol = 5) +
        ggtitle(paste(name, "PCA: timepoint & gender & exercise type")) +
        xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
        ylab(paste0("PC2, VarExp: ", percentVar[2], "%")))
    
    }
# memo: based on the PCA plots, using "norm" data would work well



#### Filtering based on intensity ####
# filter out lowly expressed genes
#. Histogram ----
medians <- rowMedians(log2(exprs(norm)))
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

#set a manual threshold to exclude genes, but not too many
# man_threshold <- round(mean(medians))  # 6
man_threshold <- 4
abline(v = man_threshold, col = "coral4", lwd = 2)

#. Excluding genes that below the threshold ----
no_of_samples <- table(pD$group)
no_of_samples
# FA FC FD FE FF MA MC MD ME MF 
# 10 10 10  8 10 10 10 10  8  8

samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(exprs(norm), 1,
                           function(x){ sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)
# FALSE  TRUE 
#    29 54646 

manfiltered <- subset(gse, idx_man_threshold)

dim(exprs(gse))  # 54675    94
dim(exprs(manfiltered))  # 54646    94



#### Annotation of the transcript clusters ####
library(hgu133plus2.db)
anno <- AnnotationDbi::select(hgu133plus2.db,
                              keys = (featureNames(manfiltered)),
                              columns = c("SYMBOL", "GENENAME"),
                              keytype = "PROBEID")
table(is.na(anno))
# FALSE   TRUE 
# 154179  20700

anno <- subset(anno, !is.na(SYMBOL))
sum(is.na(anno$SYMBOL)) # 0
nosymbols <- !(featureNames(manfiltered) %in% anno$PROBEID)
table(nosymbols)
# FALSE  TRUE 
# 44296 10350

withsymbols <- subset(manfiltered, !nosymbols)
dim(withsymbols)
# Features  Samples 
# 44296       94


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
nrow(probe_stats) # 2227: clusters that map to multiple gene symbols → remove

# remove above IDs(probe_stats)
ids_to_exlude <- (featureNames(withsymbols) %in% probe_stats$PROBEID)
table(ids_to_exlude)
# FALSE  TRUE 
# 42069  2227 

final <- subset(withsymbols, !ids_to_exlude)
validObject(final)
dim(final)  # 42069  94


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
# 42069       94



#### Differential Expression Analysis ####
head(pD, 3)
individual <- pD$id
timepoint <- pD$timepoint
group <- pD$group

#### lmFit() ####
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
    write.csv(table, file = paste0("res_GSE48278", i, ".csv"))
}


#### Visualizing results ####
library(RColorBrewer)
result_files <- list.files(pattern = ".csv")
result_files
# [1] "res_GSE48278FA.csv" "res_GSE48278FC.csv" "res_GSE48278FD.csv"
# [4] "res_GSE48278FE.csv" "res_GSE48278FF.csv" "res_GSE48278MA.csv"
# [7] "res_GSE48278MC.csv" "res_GSE48278MD.csv" "res_GSE48278ME.csv"
# [10] "res_GSE48278MF.csv"

par(mfrow = c(2,2))

for (i in 1:length(result_files)) {
    file <- result_files[i]
    results <- read.csv(file)
    group <- str_sub(file, -6, -5) 
    
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


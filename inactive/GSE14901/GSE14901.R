#### GSE14901: Aerobics, Overweight Men & Women  ####
# E-GEOD-14901

#### Loading packages ####
library(Biobase)
library(oligoClasses)
library(ArrayExpress)
library(pd.hg.u133a)
library(hgu133a.db)
library(oligo)
# library(arrayQualityMetrics)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(devtools)


#### Importing raw data from GEO ####
library(GEOquery)
geo <- getGEO("GSE14901")
length(geo)  # 1
names(geo)  # "GSE14901_series_matrix.txt.gz"

gse <- geo[[1]]
gse
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 54675 features, 72 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM373146 GSM373147 ... GSM373217 (72 total)
# varLabels: title geo_accession ... time:ch1 (38 total)
# varMetadata: labelDescription
# featureData
# featureNames: 1007_s_at 1053_at ... AFFX-TrpnX-M_at (54675 total)
# fvarLabels: ID GB_ACC ... Gene Ontology Molecular Function (16 total)
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 19654872 
# Annotation: GPL570 

head(pData(gse)$data_processing, 1)
# [1] "Global changes in gene expression were examined using the HG U133 plus gene micro-array 
# (Affymetrix, Santa Clara, CA) which contains 54,675 human probe sequences. 
# The background correction, normalization and derivation of expression measure were based on 
# the Affymetrix signal (MAS 5.0 algorithm)."

exprs(gse)[1:3, 1:3]
#           GSM373146 GSM373147 GSM373148
# 1007_s_at    5350.0    7323.1    4831.2
# 1053_at      1044.1    1049.2    1007.7
# 117_at        959.6     973.0    1334.1



#### Data cleaning ####
# featuredata ----
# head(fData(gse), 3)
# names(fData(gse))
# fData(gse) <- fData(gse)[,c("ID",  # ex. 1007_s_at
#                             "GB_ACC",
#                             "Gene Symbol",
#                             "ENTREZ_GENE_ID"
# )]
# names(fData(gse)) <- c("PROBEID", "ACCNUM", "SYMBOL","ENTREZID")
# head(fData(gse), 3)
#             PROBEID ACCNUM           SYMBOL          ENTREZID
# 1007_s_at 1007_s_at U48705 DDR1 /// MIR4640 780 /// 100616237
# 1053_at     1053_at M87338             RFC2              5982
# 117_at       117_at X51757            HSPA6              3310


# phenodata ----
pD.all <- pData(gse)
head(pD.all, 3)
tail(pD.all, 3)

pD <- pD.all[, c("title", # id
                 "gender:ch1", # gender
                 "time:ch1" # time
                 )]
names(pD)[2:3] <- c("gender", "day")

pD$gender <- as.factor(ifelse(str_detect(pD$gender, "female"), "F", "M")) # gender
pD$day <- str_replace(pD$day, " days", "")
# pD$day <- paste0("D", str_replace(pD$day, " days", ""))
pD$timepoint <- ifelse(str_detect(pD$day, "0"), "pre", 
                       ifelse(str_detect(pD$day, "14"), "post", "mid"))
pD$id <- str_sub(pD$title, -3, -1)
pD$id_time <- paste(pD$id, pD$timepoint, sep = ".")
head(pD, 3)
#                      title gender day timepoint  id  id_time
# GSM373146 Male_PRECAST_D01      M   0       pre D01  D01.pre
# GSM373147 Male_CAST02D_D01      M   2       mid D01  D01.mid
# GSM373148 Male_CAST14D_D01      M  14      post D01 D01.post

# subsetting ----
#. keep only pre/post exercise ----
dim(pD)  # 72  6
keep <- grep("pre|post", pD$timepoint)
pD <- pD[keep,]
gse <- gse[,keep]
dim(pD)  # 48  6
dim(gse)
# Features  Samples 
# 54675       48


#### normalization ####
# memo: no information about normalization in original phenodata

#### checking data intensities
oligo::boxplot(log2(exprs(gse)), target = "core",
               main = "Boxplot of log2-intensitites",
               outline = FALSE,
               las = 2,
               cex.axis=0.7)
# memo: the data is already normalized
head(pData(gse)$data_processing, 1)
# [1] "Global changes in gene expression were examined using the HG U133 plus gene micro-array 
# (Affymetrix, Santa Clara, CA) which contains 54,675 human probe sequences. 
# The background correction, normalization and derivation of expression measure were based on 
# the Affymetrix signal (MAS 5.0 algorithm)."


#### Quality control of the raw data ####
#. PCA analysis ----
myexp <- log2(exprs(gse))
PCA <- prcomp(t(myexp), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

head(pD, 3)
dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     individual = pD$id,
                     gender = pD$gender,
                     timepoint = pD$timepoint,
                     id_time = pD$id_time
)

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    facet_grid(. ~ gender) +
    ggtitle("PCA_Norm timepoint & gender") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%"))

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = gender)) +
    ggtitle("PCA_Norm gender") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%"))

ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = timepoint)) +
    ggtitle("PCA_Norm timepoint") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%"))



#### Filtering based on intensity ####
#. Histogram ----
medians <- rowMedians(log2(exprs(gse)))
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

#set a manual threshold to exclude genes, but not too many
# man_threshold <- round(mean(medians))  # 9
man_threshold <- 7
abline(v = man_threshold, col = "coral4", lwd = 2)

#. Excluding genes that below the threshold ----
# Transcripts that do not have intensities larger than the threshold in
# at least as many arrays as the smallest experimental group are excluded.
# get a list:
no_of_samples <- table(paste(pD$gender, pD$timepoint, sep = "_"))
no_of_samples
# F_post  F_pre M_post  M_pre 
# 12     12     12     12 

samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(log2(exprs(gse)), 1,
                           function(x){
                               sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)
# FALSE  TRUE 
# 3665 51010

manfiltered <- subset(gse, idx_man_threshold)

dim(exprs(gse))  # 54675    48
dim(exprs(manfiltered))  # 51010    4



#### Annotation of the transcript clusters ####
library(hgu133a.db)
anno <- AnnotationDbi::select(hgu133a.db,
                              keys = (featureNames(manfiltered)),
                              columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                              keytype = "PROBEID")
table(is.na(anno))
# FALSE   TRUE 
# 119726  92574


anno <- subset(anno, !is.na(SYMBOL))
sum(is.na(anno$SYMBOL)) # 0
nosymbols <- !(featureNames(manfiltered) %in% anno$PROBEID)
table(nosymbols)
# FALSE  TRUE 
# 20152 30858

withsymbols <- subset(manfiltered, !nosymbols)
dim(withsymbols)
# Features  Samples 
# 20152       48


#. Removing multiple mappings ----
# group by PROBEID
anno_grouped <- group_by(anno, PROBEID)
head(anno_grouped)
# A tibble: 6 × 4
# Groups:   PROBEID [5]
#   PROBEID   SYMBOL  GENENAME                                     ENTREZID 
#     <chr>     <chr>   <chr>                                        <chr>    
# 1 1007_s_at DDR1    discoidin domain receptor tyrosine kinase 1  780      
# 2 1007_s_at MIR4640 microRNA 4640                                100616237
# 3 1053_at   RFC2    replication factor C subunit 2               5982     
# 4 117_at    HSPA6   heat shock protein family A (Hsp70) member 6 3310     
# 5 121_at    PAX8    paired box 8                                 7849     
# 6 1255_g_at GUCA1A  guanylate cyclase activator 1A               2978


# summarize
anno_summarized <-
    dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized)
# A tibble: 6 × 2
#   PROBEID   no_of_matches
#   <chr>             <int>
# 1 1007_s_at             2
# 2 1053_at               1
# 3 117_at                1
# 4 121_at                1
# 5 1255_g_at             1
# 6 1294_at               2


# filter to extract cluster multiple gene that map to multiple gene symbols
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
head(anno_filtered)

probe_stats <- anno_filtered
nrow(probe_stats) # 1168: clusters that map to multiple gene symbols → remove

# remove above IDs(probe_stats)
ids_to_exlude <- (featureNames(withsymbols) %in% probe_stats$PROBEID)
table(ids_to_exlude)
# FALSE  TRUE 
# 18984  1168

final <- subset(withsymbols, !ids_to_exlude)
validObject(final)
dim(final)  # 18984    48 


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
# 18984       48



#### Differential Expression Analysis ####
head(pD, 3)
individual <- pD$id
timepoint <- pD$timepoint
gender <- pD$gender

#### lmFit() ####
exprs(final) <- log2(exprs(final))

# check fData to avoid the result table to be too large
head(fData(final))
ncol(fData(final))
names(fData(final))
fData(final) <- fData(final)[17:20] # keeping necessary columns
head(fData(final))

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
    write.csv(table, file = paste0("res_GSE14901", i, ".csv"))
    print(head(table))
}


#### Visualizing results ####
library(RColorBrewer)
result_files <- list.files(pattern = ".csv")
result_files
# [1] "res_GSE14901F.csv" "res_GSE14901M.csv"
par(mfrow = c(2,2))

#. histgram (p $ adjp vals) ----
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
    cat("adj.P < 0.1:", nrow(subset(results, adj.P.Val < 0.05)),"\n")
    cat("adj.P < 0.05:", nrow(subset(results, adj.P.Val < 0.05)),"\n")
    cat("adj.P < 0.01:", nrow(subset(results, adj.P.Val < 0.01)),"\n\n")
}

## Results  
# 7/29/2021  
# res_GSE14901F.csv 
# p < 0.05: 3183 
# p < 0.01: 1301 
# adj.P < 0.1: 295 
# adj.P < 0.05: 295 
# adj.P < 0.01: 4 
# 
# res_GSE14901M.csv 
# p < 0.05: 2795 
# p < 0.01: 1128 
# adj.P < 0.1: 237 
# adj.P < 0.05: 237 
# adj.P < 0.01: 2


# to explore more (ex)
# tail(subset(results, adj.P.Val < 0.05))



#. volcano plot ----
library(EnhancedVolcano)
library(patchwork)

for (i in 1:length(result_files)) {
    file <- result_files[i]
    results <- read.csv(file)
    name <- str_sub(file, 5, 13)
    group <- tolower(str_sub(file, -5, -5))
    
    p1 <-
    EnhancedVolcano(results,
                    lab = results$SYMBOL,
                    # selectLab = 'Peptide9A',
                    x = 'logFC',
                    y = 'adj.P.Val',
                    xlim = c(floor(min(results$logFC)), ceiling(max(abs(results$logFC)))),
                    ylim = c(0, ceiling(max(-log10(results$adj.P.Val)))),
                    xlab = bquote(~Log[2]~ 'fold change'),
                    ylab = bquote(~-Log[10]~adjusted~italic(P)),
                    pCutoff = 0.01,
                    FCcutoff = 0.5,
                    pointSize = 2.0, 
                    labSize = 3.0,
                    col = c('black', 'red', 'orange', 'blue'),
                    legendPosition = '',
                    legendLabels = c('NS','Log2 FC',
                                     'adj.p.val','adj.p.val & Log2 FC')) +
        
        labs(title = paste(name, "adj.p < 0.01"),
             titleLabSize = 10,
             subtitle = "pre vs post",
             caption = "")
    
    
    p2 <- 
    EnhancedVolcano(results,
                    lab = results$SYMBOL,
                    # selectLab = 'Peptide9A',
                    x = 'logFC',
                    y = 'adj.P.Val',
                    xlim = c(floor(min(results$logFC)), ceiling(max(abs(results$logFC)))),
                    ylim = c(0, ceiling(max(-log10(results$adj.P.Val)))),
                    xlab = bquote(~Log[2]~ 'fold change'),
                    ylab = bquote(~-Log[10]~adjusted~italic(P)),
                    pCutoff = 0.05,
                    FCcutoff = 0.5,
                    pointSize = 2.0, 
                    labSize = 3.0,
                    col = c('black', 'red', 'orange', 'blue'),
                    legendPosition = '',
                    legendLabels = c('NS','Log2 FC',
                                     'adj.p.val','adj.p.val & Log2 FC')) +
        
        labs(title = paste(name, "adj.p < 0.05"),
             titleLabSize = 10,
             subtitle = "pre vs post",
             caption = "")

        
    p3 <- 
    EnhancedVolcano(results,
                    lab = results$SYMBOL,
                    # selectLab = 'Peptide9A',
                    x = 'logFC',
                    y = 'adj.P.Val',
                    xlim = c(floor(min(results$logFC)), ceiling(max(abs(results$logFC)))),
                    ylim = c(0, ceiling(max(-log10(results$adj.P.Val)))),
                    xlab = bquote(~Log[2]~ 'fold change'),
                    ylab = bquote(~-Log[10]~adjusted~italic(P)),
                    pCutoff = 0.1,
                    FCcutoff = 0.5,
                    pointSize = 2.0, 
                    labSize = 3.0,
                    col = c('black', 'red', 'orange', 'blue'),
                    legendPosition = '',
                    legendLabels = c('NS','Log2 FC',
                                     'adj.p.val','adj.p.val & Log2 FC')) +
        
        labs(title = paste(name, "adj.p < 0.1"),
             titleLabSize = 10,
             subtitle = "pre vs post",
             caption = "")
    
    
    p4 <- 
    EnhancedVolcano(results,
                    lab = results$SYMBOL,
                    # selectLab = 'Peptide9A',
                    x = 'logFC',
                    y = 'P.Value',
                    xlim = c(floor(min(results$logFC)), ceiling(max(abs(results$logFC)))),
                    ylim = c(0, ceiling(max(-log10(results$P.Value)))),
                    xlab = bquote(~Log[2]~ 'fold change'),
                    ylab = bquote(~-Log[10]~italic(P)),
                    pCutoff = 0.001,
                    FCcutoff = 0.5,
                    pointSize = 2.0, 
                    labSize = 3.0,
                    col = c('black', 'red', 'orange', 'blue'),
                    legendPosition = '',
                    legendLabels = c('NS','Log2 FC', 'p.val','p.val & Log2 FC')) +
        labs(title = paste(name, "p < 0.001"),
             titleLabSize = 10,
             subtitle = "pre vs post",
             caption = "")
    
    
    p5 <- 
    EnhancedVolcano(results,
                    lab = results$SYMBOL,
                    # selectLab = 'Peptide9A',
                    x = 'logFC',
                    y = 'P.Value',
                    xlim = c(floor(min(results$logFC)), ceiling(max(abs(results$logFC)))),
                    ylim = c(0, ceiling(max(-log10(results$P.Value)))),
                    xlab = bquote(~Log[2]~ 'fold change'),
                    ylab = bquote(~-Log[10]~italic(P)),
                    pCutoff = 0.01,
                    FCcutoff = 0.5,
                    pointSize = 2.0, 
                    labSize = 3.0,
                    col = c('black', 'red', 'orange', 'blue'),
                    legendPosition = '',
                    legendLabels = c('NS','Log2 FC', 'p.val','p.val & Log2 FC')) +
        labs(title = paste(name, "p < 0.01"),
             titleLabSize = 10,
             subtitle = "pre vs post",
             caption = "")
    
    
    p6 <- 
    EnhancedVolcano(results,
                    lab = results$SYMBOL,
                    # selectLab = 'Peptide9A',
                    x = 'logFC',
                    y = 'P.Value',
                    xlim = c(floor(min(results$logFC)), ceiling(max(abs(results$logFC)))),
                    ylim = c(0, ceiling(max(-log10(results$P.Value)))),
                    xlab = bquote(~Log[2]~ 'fold change'),
                    ylab = bquote(~-Log[10]~italic(P)),
                    pCutoff = 0.05,
                    FCcutoff = 0.5,
                    pointSize = 2.0, 
                    labSize = 3.0,
                    col = c('black', 'red', 'orange', 'blue'),
                    legendPosition = '',
                    legendLabels = c('NS','Log2 FC', 'p.val','p.val & Log2 FC')) +
        labs(title = paste(name, "p < 0.05"),
             titleLabSize = 10,
             subtitle = "pre vs post",
             caption = "")
    
    png(filename = paste0("plot/volcano_adjp_", group, ".png"),
        width = 1280, height = 600)
    print(p1 + p2 + p3  + plot_layout(guides = "collect") & 
        theme(legend.position = "bottom", 
              legend.text = element_text(size = 10)))
    dev.off()
    
    png(filename = paste0("plot/volcano_p_", group, ".png"),
        width = 1280, height = 600)
    print(p4 + p5 + p6  + plot_layout(guides = "collect") & 
              theme(legend.position = "bottom", 
                    legend.text = element_text(size = 10)))
    dev.off()
}
### 特定のgeneを取り出したplot図
# topGene <- rownames(res)[663] #86: 128B, 525: SHMOOSE, 148: Humanin, 663: MOTSc, SHLP6, 644: SHLP2
# plotCounts(dds, topGene, intgroup = c("sarcopenia_status"), xlab = "Group", ylabel = "Normalized Counts") 



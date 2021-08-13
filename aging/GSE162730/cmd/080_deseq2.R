# ref: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
setwd("/Volumes/HDD24TB/MetaAnalysisProject_Apr2021/aging/GSE162730")
# young vs old

#### helpers ####
library("stringr")
library("dplyr")
library("tidyr")
library("magrittr")
library("ggplot2")

#### Importing Data ####
#### 8/13/2021 Rna
samples <- read.table("sampletable_rna.txt", sep = "\t", header = TRUE)
head(samples, 3)
#           run subject library time gender
# 1 SRR13202581     S21 RNA-Seq  pre female
# 2 SRR13202584     S21 RNA-Seq post female
# 3 SRR13202593     S23 RNA-Seq  pre   male

counts <- read.table("countmat_rna.txt", sep = "\t", header = TRUE)
counts[1:5,1:5]
#                 SRR13202581 SRR13202584 SRR13202593 SRR13202596 SRR13202605
# ENSG00000284662        0.86        0.31        0.25        0.34        0.45
# ENSG00000186827       20.63       10.25        5.17        5.27       14.17
# ENSG00000186891        5.33        4.00        0.00        2.03        1.00
# ENSG00000160072      425.39      197.80      129.01      165.94      225.18
# ENSG00000041988      158.07       76.25       56.45       61.71       80.95

dim(samples)  # 16  5
dim(counts)  # 60664    16


#### Constructing DESeqDataSet (DESeqDataSetFromMatrix) ####
library("DESeq2")
gender_id = unique(samples$gender)
dds <- list()
for (i in gender_id) {
    print(i)
    countdata <- counts[,samples$run[samples$gender == i]]
    coldata <- samples[samples$gender == i,]
    dds[[i]] <- DESeqDataSetFromMatrix(countData = round(countdata),
                                       colData = coldata, 
                                       design = ~ subject + time)
    print(dds[[i]])
    cat("\n")
}
# [1] "female"
# converting counts to integer mode
# class: DESeqDataSet 
# dim: 60664 12 
# metadata(1): version
# assays(1): counts
# rownames(60664): ENSG00000284662 ENSG00000186827 ... ENSG00000277475 ENSG00000275405
# rowData names(0):
#   colnames(12): 1 2 ... 13 14
# colData names(5): run subject library time gender
# 
# [1] "male"
# converting counts to integer mode
# class: DESeqDataSet 
# dim: 60664 4 
# metadata(1): version
# assays(1): counts
# rownames(60664): ENSG00000284662 ENSG00000186827 ... ENSG00000277475 ENSG00000275405
# rowData names(0):
#   colnames(4): 3 4 15 16
# colData names(5): run subject library time gender


#### Filtering low counts data ####
# minimal filtering
for (i in gender_id) {
    cat("dds", i, "\n")
    keep <- rowSums(counts(dds[[i]])) > 1
    dds[[i]] <- dds[[i]][keep,]
    cat(dim(dds[[i]]), "\n")
    cat(dim(colData(dds[[i]])), "\n\n")
}
# dds female 
# 47705 12 
# 12 5 
# 
# dds male 
# 43558 4 
# 4 5


# cf. at least 3 samples with a count of 10 or higher
# keep2 <- rowSums(counts(dds) >= 10) >= 3
# table(keep2)



#### Exploratory analysis and visualization ####
if (!file.exists("plot")) {
    dir.create("plot")
}

# 1. Normalization ---
for (i in gender_id) {
    gender_letter <- str_sub(i, 1, 1)
    dds[[i]] <- estimateSizeFactors(dds[[i]])
    
    # 2-1. box-plot ----
    png(paste0("plot/boxplot_check_normalization_", gender_letter, ".png"), width = 1280, height = 800)
    par(mfrow = c(1, 2), mar = c(5,5,4,2))
    boxplot(log2(counts(dds[[i]]) + 1), cex.axis = 0.7, ## color 
            las = 1, xlab = "log2(counts)", horizontal = TRUE, main = paste0(i, " Raw counts"))
    boxplot(log2(counts(dds[[i]], normalized = TRUE) + 1), cex.axis = 0.7,  ## color 
            las = 1, xlab = "log2(normalized counts)", horizontal = TRUE, main = "Normalized counts") 
    dev.off()
    
    # 2.2 density-plot ----
    library(affy)
    png(paste0("plot/densityplot_check_normalization_", gender_letter, ".png"), width = 1280, height = 800)
    par(mfrow = c(1, 2), mar = c(5,5,4,2))
    plotDensity(log2(counts(dds[[i]]) + 1),
                xlab = paste0(i, " log2(counts)"), cex.lab = 0.7, panel.first = grid()) 
    plotDensity(log2(counts(dds[[i]], normalized = TRUE) + 1),
                xlab = "log2(normalized counts)", cex.lab = 0.7, panel.first = grid())
    dev.off()
}


# 3. vst(): Variance Stabilizing Transformation ----
# recommended for n > 30 dataset
vsd = list()
for (i in gender_id) {
    vsd[[i]] <- vst(dds[[i]], blind = FALSE)
    cat("vsd ", i, "\n")
    print(head(assay(vsd[[i]]), 3))
    cat("\n")
}
# blind = FALSE: means that differences between cell lines and treatment 
# (the variables in the design) will not contribute to the expected variance-mean trend of the experiment

# these will be used soon
# vsd  female 
#                        1        2        5        6        7        8        9       10       11
# ENSG00000284662 5.299515 5.166957 5.166957 5.347628 5.406577 5.354960 5.166957 5.166957 5.166957
# ENSG00000186827 5.770222 5.757756 5.822295 5.675666 5.613970 5.786171 5.836940 5.924704 5.693526
# ENSG00000186891 5.462950 5.542169 5.343503 5.347628 5.406577 5.492130 5.166957 5.166957 5.431335
#                       12       13       14
# ENSG00000284662 5.166957 5.166957 5.166957
# ENSG00000186827 5.598819 5.832102 5.537712
# ENSG00000186891 5.520008 5.166957 5.166957
# 
# vsd  male 
#                        3        4       15       16
# ENSG00000186827 7.932091 7.910489 7.966084 7.952073
# ENSG00000186891 7.679663 7.825743 7.770371 7.679663
# ENSG00000160072 8.924499 8.967867 8.980761 9.048592



# cf. rlog()
# n < 30 dataset might go well with rlog()
# rld <- rlog(dds, blind = FALSE)
# head(assay(rld), 3)


# 4. PCA plot ----
pcaData = list()
for (i in gender_id) {
    gender_letter <- str_sub(i, 1, 1)
    cat("calculating PCA for dds", i, "\n")
    pcaData[[i]] <- plotPCA(vsd[[i]], intgroup = c("time"), returnData = TRUE)
    print(head(pcaData[[i]]))
    percentVar <- round(100 * attr(pcaData[[i]], "percentVar"))
    
    pca_df <- data.frame(pcaData[[i]],
                         time = samples$time[samples$gender == i]
    )
    
    ### pca plot rewite
    png(paste0("plot/pca_", gender_letter, ".png"))
    print(
        ggplot(pca_df, aes(x = PC1, y = PC2, color = time)) +
            geom_point(size = 3) +
            xlab(paste0("PC1: ", percentVar[1], "% variance")) +
            ylab(paste0("PC2: ", percentVar[2], "% variance")) +
            coord_fixed() +
            ggtitle(paste0("PCA with VST data ", i))
    )
    dev.off()
    
    cat("\n")
}
# calculating PCA for dds female 
#          PC1       PC2 group time name
# 1 -5.8925172  8.718160   pre  pre    1
# 2 -0.7584437  9.897606  post post    2
# 5 -5.1521836 -7.933238   pre  pre    5
# 6 -9.5585536 -8.249404  post post    6
# 7  2.5663501 -3.683785   pre  pre    7
# 8  7.5281312 -6.484747  post post    8
# 
# calculating PCA for dds male 
#           PC1       PC2 group time name
# 3  -11.324007  3.173347   pre  pre    3
# 4   -5.142259 -4.176294  post post    4
# 15   8.858219  5.462419   pre  pre   15
# 16   7.608048 -4.459473  post post   16



#### DESeq2 ####
res = list()
for (i in gender_id) {
    print(i)
    dds[[i]]$time <- relevel(as.factor(dds[[i]]$time), ref = "pre")
    dds[[i]] <- DESeq(dds[[i]])
    res[[i]] <- results(dds[[i]])  # results(dds, contrast = c("time", "post", "pre"))
    print(head(res[[i]]))
    cat("\n")
    
    summary(res[[i]])  # FDR < 0.1
    cat("\n")
}
# [1] "female"
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE       stat    pvalue      padj
# <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
# ENSG00000284662   0.205998     0.03509519  3.085594  0.0113739  0.990925  0.999995
# ENSG00000186827   5.917091    -0.15466780  0.519766 -0.2975718  0.766030  0.999995
# ENSG00000186891   0.938536     0.49205360  1.199766  0.4101247  0.681714  0.999995
# ENSG00000160072 125.967796     0.00880607  0.139090  0.0633120  0.949518  0.999995
# ENSG00000041988  47.036213    -0.08747999  0.188047 -0.4652036  0.641786  0.999995
# ENSG00000260179   3.643064    -0.05506193  0.686311 -0.0802289  0.936055  0.999995
# 
# out of 47705 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.0021%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# 
# 
# [1] "male"
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE        stat    pvalue      padj
# <numeric>      <numeric> <numeric>   <numeric> <numeric> <numeric>
# ENSG00000186827   6.747917     -0.1993291  1.403643 -0.14200842  0.887073        NA
# ENSG00000186891   0.728657      0.5081277  4.604596  0.11035231  0.912130        NA
# ENSG00000160072 178.552910      0.1316133  0.322689  0.40786418  0.683373        NA
# ENSG00000041988  71.720699      0.0446927  0.468669  0.09536094  0.924028        NA
# ENSG00000260179   4.392608      1.0323010  1.736686  0.59440840  0.552239        NA
# ENSG00000234396   0.431803     -0.0316639  4.887701 -0.00647828  0.994831        NA
# 
# out of 43558 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 11, 0.025%
# LFC < 0 (down)     : 26, 0.06%
# outliers [1]       : 0, 0%
# low counts [2]     : 38846, 89%
# (mean count < 388)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results



#### results memo ####
### head(res[[m]]) ----
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 6 rows and 6 columns
#                   baseMean log2FoldChange     lfcSE      stat    pvalue      padj
#                   <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# ENSG00000284662   5.305012      0.5280580 0.3683311  1.433650 0.1516721        NA
# ENSG00000186827   9.897404      0.6976768 0.3402539  2.050459 0.0403196        NA
# ENSG00000186891   0.686316     -0.5609034 1.1213943 -0.500184 0.6169456        NA
# ENSG00000160072 164.351008     -0.0483656 0.0869517 -0.556235 0.5780501  0.841077
# ENSG00000041988  97.857617      0.0437376 0.0953073  0.458911 0.6462980  0.872899
# ENSG00000260179   2.859081      0.2892322 0.5065239  0.571014 0.5679902        NA

### summary(res[[m]]) ----
# out of 43164 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 842, 2%
# LFC < 0 (down)     : 1036, 2.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 24269, 56%
# (mean count < 10)


### head(res[[f]]) ----
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 6 rows and 6 columns
#                   baseMean log2FoldChange     lfcSE       stat    pvalue      padj
#                   <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
# ENSG00000284662   3.900494     0.04528241  0.432045  0.1048094  0.916527        NA
# ENSG00000186827   8.581578    -0.09168495  0.384861 -0.2382285  0.811704        NA
# ENSG00000186891   0.810169     0.54240961  1.075564  0.5043024  0.614049        NA
# ENSG00000160072 107.053323    -0.00789172  0.110231 -0.0715924  0.942926  0.978687
# ENSG00000041988  80.218381     0.13116801  0.131364  0.9985100  0.318032  0.601431
# ENSG00000260179   1.486128     0.30203507  0.660123  0.4575435  0.647280        NA

### summary(res[[f]]) ----
# out of 42250 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1033, 2.4%
# LFC < 0 (down)     : 1382, 3.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 24574, 58%
# (mean count < 10)


# explore results ----
library(RColorBrewer)
resSig = list()
for (i in gender_id) {
    gender_letter <- str_sub(i, 1, 1)  
    # lower the FDR to 0.05
    cat("results for", i, "\n")
    results <- res[[i]]
    cat("p < 0.05:", nrow(subset(results, pvalue < 0.05)), "\n")
    cat("p < 0.01:", nrow(subset(results, pvalue < 0.01)), "\n")
    cat("p < 0.001:", nrow(subset(results, pvalue < 0.001)), "\n")
    cat("adj.P < 0.1:", nrow(subset(results, padj < 0.1)),"\n")
    cat("adj.P < 0.05:", nrow(subset(results, padj < 0.05)),"\n")
    cat("adj.P < 0.01:", nrow(subset(results, padj < 0.01)),"\n\n")
    
    # subset the res, sort it by the log2 fold change:
    resSig[[i]] <- subset(results, padj < 0.1)
    ressig <- resSig[[i]]
    # strongest down-regulation
    cat("strongest upregulation: adj.p < 0.1", "\n")
    print(head(ressig[order(ressig$log2FoldChange),]))
    cat("\n")
    
    # strongest upregulation
    cat("strongest upregulation: adj.p < 0.1", "\n")
    print(head(ressig[order(-ressig$log2FoldChange),]))
    cat("\n")
    
    ## histgram ##
    png(paste0("plot/hist_p_adjp_", gender_letter, ".png"), width = 1280, height = 800)
    par(mfrow = c(1,2))
    hist(results$pvalue, col = brewer.pal(3, name = "Set2")[1], 
         main = paste(i, "Pval"), xlab  = NULL)
    hist(results$padj, col = brewer.pal(3, name = "Set2")[2],
         main = paste(i, "adj.Pval"), xlab = NULL)
    dev.off()
}

## memo ##
# results for female 
# p < 0.05: 370 
# p < 0.01: 83 
# p < 0.001: 23 
# adj.P < 0.1: 1 
# adj.P < 0.05: 1 
# adj.P < 0.01: 1 
# 
# strongest upregulation: adj.p < 0.1 
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 1 row and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# ENSG00000141526   379.068      -0.605803  0.103204  -5.86996 4.35901e-09 0.000207946
# 
# strongest upregulation: adj.p < 0.1 
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 1 row and 6 columns
#                  baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#                 <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# ENSG00000141526   379.068      -0.605803  0.103204  -5.86996 4.35901e-09 0.000207946
# 
# results for male 
# p < 0.05: 345 
# p < 0.01: 129 
# p < 0.001: 45 
# adj.P < 0.1: 37 
# adj.P < 0.05: 13 
# adj.P < 0.01: 4 
# 
# strongest upregulation: adj.p < 0.1 
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 6 rows and 6 columns
#                  baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#                 <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# ENSG00000166147  2364.727       -1.22481  0.214592  -5.70763 1.14560e-08 5.39805e-05
# ENSG00000091986   520.469       -1.17224  0.262283  -4.46937 7.84507e-06 9.24149e-03
# ENSG00000038427   527.286       -1.13040  0.241148  -4.68756 2.76482e-06 4.34260e-03
# ENSG00000077009  1232.935       -1.06921  0.249428  -4.28666 1.81376e-05 1.70929e-02
# ENSG00000181418  1261.963       -1.05848  0.211312  -5.00910 5.46852e-07 1.28838e-03
# ENSG00000199990  1076.376       -1.01623  0.249392  -4.07482 4.60492e-05 2.39909e-02
# 
# strongest upregulation: adj.p < 0.1 
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 6 rows and 6 columns
#                  baseMean log2FoldChange     lfcSE      stat      pvalue      padj
#                 <numeric>      <numeric> <numeric> <numeric>   <numeric> <numeric>
# ENSG00000269736   866.097       1.284974  0.310532   4.13797 3.50395e-05 0.0235866
# ENSG00000283657   865.949       0.903511  0.250194   3.61124 3.04734e-04 0.0668574
# ENSG00000240474   554.033       0.859089  0.231913   3.70436 2.11922e-04 0.0554763
# ENSG00000215187   497.047       0.851491  0.236197   3.60500 3.12153e-04 0.0668574
# ENSG00000234594   863.119       0.821293  0.242148   3.39170 6.94602e-04 0.0935914
# ENSG00000123689  1826.808       0.818497  0.196593   4.16340 3.13537e-05 0.0235866



#### Annotation & Write out the results ####
library("AnnotationDbi")
library("Homo.sapiens")
# columns(Homo.sapiens)
for (i in gender_id) {
    gender_letter <- str_sub(i, 1, 1)  
    cat("annotating", i, "\n")  
    results <- resSig[[i]]  
    results$symbol <- mapIds(Homo.sapiens,
                             keys = row.names(results),
                             column = "SYMBOL",
                             keytype = "ENSEMBL" # change
    )
    cat("head(results, 3):", "\n")
    print(head(results, 3))
    cat("\n")
    
    cat("dim(results):", "\n")
    print(dim(results))
    cat("\n\n")
    # keep <- res$symbol[!duplicated(res$symbol)]
    # dup_symbol <- res$symbol[duplicated(res$symbol)]
    # length(keep)  # [1] 17752
    # length(dup_symbol)  # [1] 248
    # res <- subset(res,!duplicated(res$symbol))
    # dim(res)  # [1] 17752     7
    
    write.csv(results, paste0("resGSE162730", gender_letter, ".csv"))
}

## memo ##
# annotating female 
# 'select()' returned 1:1 mapping between keys and columns
# head(results, 3): 
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 1 row and 7 columns
#                  baseMean log2FoldChange     lfcSE      stat      pvalue        padj      symbol
#                 <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric> <character>
# ENSG00000141526   379.068      -0.605803  0.103204  -5.86996 4.35901e-09 0.000207946     SLC16A3
# 
# dim(results): 
# [1] 1 7
# 
# 
# annotating male 
# 'select()' returned 1:1 mapping between keys and columns
# head(results, 3): 
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 3 rows and 7 columns
#                  baseMean log2FoldChange     lfcSE      stat      pvalue      padj      symbol
#                 <numeric>      <numeric> <numeric> <numeric>   <numeric> <numeric> <character>
# ENSG00000123689  1826.808       0.818497  0.196593   4.16340 3.13537e-05 0.0235866        G0S2
# ENSG00000202031  1442.986      -0.745542  0.194097  -3.84108 1.22497e-04 0.0449423    SNORD38A
# ENSG00000207241   562.184      -0.803471  0.231473  -3.47112 5.18288e-04 0.0897698    SNORD45A
# 
# dim(results): 
# [1] 37  7



#### Volcano plot ####
library(EnhancedVolcano)
library(patchwork)

# str_view(list.files(), "res.*csv")
result_files <- list.files(pattern = "res.*csv")
result_files
# "resGSE162730f.csv" "resGSE162730m.csv"
for (i in 1:length(result_files)) {
    file <- result_files[i]
    results <- read.csv(file)
    name <- str_sub(file, 4, 12)
    group <- str_sub(file, -5, -5)
    
    p1 <-
        EnhancedVolcano(results,
                        lab = results$symbol,
                        # selectLab = 'Peptide9A',
                        x = 'log2FoldChange',
                        y = 'padj',
                        xlim = c(floor(min(results$log2FoldChange)), ceiling(max(abs(results$log2FoldChange)))),
                        ylim = c(0, ceiling(max(-log10(results$padj)))),
                        xlab = bquote(~Log[2]~ 'fold change'),
                        ylab = bquote(~-Log[10]~adjusted~italic(P)),
                        pCutoff = 0.01,
                        FCcutoff = 0.5,
                        pointSize = 2.0, 
                        labSize = 3.0,
                        col = c('black', 'red', 'orange', 'blue'),
                        legendPosition = '',
                        legendLabels = c('NS','Log2 FC',
                                         'padj','padj & Log2 FC')) +
        
        labs(title = paste(name, "adj.p < 0.01"),
             titleLabSize = 10,
             subtitle = "pre vs post",
             caption = "")
    
    
    p2 <- 
        EnhancedVolcano(results,
                        lab = results$symbol,
                        # selectLab = 'Peptide9A',
                        x = 'log2FoldChange',
                        y = 'padj',
                        xlim = c(floor(min(results$log2FoldChange)), ceiling(max(abs(results$log2FoldChange)))),
                        ylim = c(0, ceiling(max(-log10(results$padj)))),
                        xlab = bquote(~Log[2]~ 'fold change'),
                        ylab = bquote(~-Log[10]~adjusted~italic(P)),
                        pCutoff = 0.05,
                        FCcutoff = 0.5,
                        pointSize = 2.0, 
                        labSize = 3.0,
                        col = c('black', 'red', 'orange', 'blue'),
                        legendPosition = '',
                        legendLabels = c('NS','Log2 FC',
                                         'padj','padj & Log2 FC')) +
        
        labs(title = paste(name, "adj.p < 0.05"),
             titleLabSize = 10,
             subtitle = "pre vs post",
             caption = "")
    
    
    p3 <- 
        EnhancedVolcano(results,
                        lab = results$symbol,
                        # selectLab = 'Peptide9A',
                        x = 'log2FoldChange',
                        y = 'padj',
                        xlim = c(floor(min(results$log2FoldChange)), ceiling(max(abs(results$log2FoldChange)))),
                        ylim = c(0, ceiling(max(-log10(results$padj)))),
                        xlab = bquote(~Log[2]~ 'fold change'),
                        ylab = bquote(~-Log[10]~adjusted~italic(P)),
                        pCutoff = 0.1,
                        FCcutoff = 0.5,
                        pointSize = 2.0, 
                        labSize = 3.0,
                        col = c('black', 'red', 'orange', 'blue'),
                        legendPosition = '',
                        legendLabels = c('NS','Log2 FC',
                                         'padj','padj & Log2 FC')) +
        
        labs(title = paste(name, "adj.p < 0.1"),
             titleLabSize = 10,
             subtitle = "pre vs post",
             caption = "")
    
    
    p4 <- 
        EnhancedVolcano(results,
                        lab = results$symbol,
                        # selectLab = 'Peptide9A',
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        xlim = c(floor(min(results$log2FoldChange)), ceiling(max(abs(results$log2FoldChange)))),
                        ylim = c(0, ceiling(max(-log10(results$pvalue)))),
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
                        lab = results$symbol,
                        # selectLab = 'Peptide9A',
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        xlim = c(floor(min(results$log2FoldChange)), ceiling(max(abs(results$log2FoldChange)))),
                        ylim = c(0, ceiling(max(-log10(results$pvalue)))),
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
                        lab = results$symbol,
                        # selectLab = 'Peptide9A',
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        xlim = c(floor(min(results$log2FoldChange)), ceiling(max(abs(results$log2FoldChange)))),
                        ylim = c(0, ceiling(max(-log10(results$pvalue)))),
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
        width = 1280, height = 800)
    print(p1 + p2 + p3  + plot_layout(guides = "collect") & 
              theme(legend.position = "bottom", 
                    legend.text = element_text(size = 10)))
    dev.off()
    
    png(filename = paste0("plot/volcano_p_", group, ".png"),
        width = 1280, height = 800)
    print(p4 + p5 + p6  + plot_layout(guides = "collect") & 
              theme(legend.position = "bottom", 
                    legend.text = element_text(size = 10)))
    dev.off()
}


################ Explore Results ################
#table(df_cor$p < 0.1)
# dfSig <- subset(res, pvalue < 0.01)
# dfSigOrdered <- dfSig[order(dfSig$pvalue),]
# head(dfSigOrdered)
# dim(dfSigOrdered)  # [1] 3262    7
# 
# dfSigadj5 <- subset(res, padj < 0.05)
# dfSig5Ordered <- dfSigadj5[order(dfSigadj5$padj),]
# head(dfSig5Ordered)
# dim(dfSig5Ordered)  # [1] 3018    7
# 
# dfSigadj01 <- subset(res, padj < 0.01)
# dfSig01Ordered <- dfSigadj01[order(dfSigadj01$padj),]
# head(dfSig01Ordered)
# dim(dfSig01Ordered)  # [1] 1856    7
###################################################
###################################################



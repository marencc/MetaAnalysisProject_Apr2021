# ref: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
setwd("/Volumes/HDD24TB/MetaAnalysisProject_Apr2021/inactive/GSE113165")

#### helpers ####
library("stringr")
library("dplyr")
library("tidyr")
library("magrittr")
library("ggplot2")

#### Importing Data ####
# samples <- read.table("pdata.txt", sep = "\t", header = TRUE)
samples <- read.table("sampletable.txt", sep = "\t", header = TRUE)
head(samples, 3)
#                   run subject gender age time
# SRR7007949 SRR7007949   S3005   male old  pre
# SRR7007950 SRR7007950   S3005   male old post
# SRR7007951 SRR7007951   S3008   male old  pre

counts <- read.table("countmat.fastp.txt", sep = "\t", header = TRUE)
counts[1:5, 1:5]
#                 SRR7007949 SRR7007950 SRR7007951 SRR7007952 SRR7007953
# ENSG00000284662          9          4          2          4         12
# ENSG00000186827          3         23         14         11          7
# ENSG00000186891          1          2          0          0          1
# ENSG00000160072        154        168        153        198        195
# ENSG00000041988         88         78         99        132        115
dim(counts)  # 60664    56



#### Constructing DESeqDataSet (DESeqDataSetFromMatrix) ####
library("DESeq2")
samples$age <- relevel(factor(samples$age), ref = "young")
gender_id = unique(samples$gender)
dds <- list()
for (i in gender_id) {
  print(i)
  coldata <- samples[samples$gender == i,]
  countdata <- counts[,rownames(coldata)]
  dds[[i]] <- DESeqDataSetFromMatrix(countData = countdata,
                                     colData = coldata, 
                                     design = ~ subject + age)
  print(dds[[i]])
  cat("\n")
}
# [1] "male"
# class: DESeqDataSet 
# dim: 60664 26 
# metadata(1): version
# assays(1): counts
# rownames(60664): ENSG00000284662 ENSG00000186827 ... ENSG00000277475 ENSG00000275405
# rowData names(0):
#   colnames(26): SRR7007949 SRR7007950 ... SRR7007993 SRR7007994
# colData names(8): run age ... time color
# 
# [1] "female"
# class: DESeqDataSet 
# dim: 60664 30 
# metadata(1): version
# assays(1): counts
# rownames(60664): ENSG00000284662 ENSG00000186827 ... ENSG00000277475 ENSG00000275405
# rowData names(0):
#   colnames(30): SRR7007959 SRR7007960 ... SRR7008003 SRR7008004
# colData names(8): run age ... time color


#### Filtering low counts data ####
# minimal filtering
for (i in gender_id) {
  cat("dds", i, "\n")
  keep <- rowSums(counts(dds[[i]])) > 1
  dds[[i]] <- dds[[i]][keep,]
  cat(dim(dds[[i]]), "\n")  # 45771    56
  cat(dim(colData(dds[[i]])), "\n\n")  # 56  8
}
# dds male 
# 43164 26 
# 26 5 
# 
# dds female 
# 42250 30 
# 30 5 


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
  boxplot(log2(counts(dds[[i]]) + 1), cex.axis = 0.7, 
          las = 1, xlab = "log2(counts)", horizontal = TRUE, main = paste0(i, " Raw counts"))
  boxplot(log2(counts(dds[[i]], normalized = TRUE) + 1), cex.axis = 0.7, 
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
  print(assay(vsd[[i]])[1:5, 1:5])
  cat("\n")
}
# blind = FALSE: means that differences between cell lines and treatment 
# (the variables in the design) will not contribute to the expected variance-mean trend of the experiment

# these will be used for PCA plot
# vsd  male 
# SRR7007949 SRR7007950 SRR7007951 SRR7007952 SRR7007953
# ENSG00000284662   5.888296   5.625707   5.483992   5.606308   5.918535
# ENSG00000186827   5.578088   6.271543   6.026989   5.903236   5.739142
# ENSG00000186891   5.396487   5.486257   5.146958   5.146958   5.372128
# ENSG00000160072   7.839329   7.850530   7.738680   7.936914   7.871279
# ENSG00000041988   7.288550   7.114665   7.314042   7.520123   7.345189
# 
# vsd  female 
# SRR7007959 SRR7007960 SRR7007971 SRR7007972 SRR7007973
# ENSG00000284662   5.388084   5.357162   5.453702   5.372651   5.514264
# ENSG00000186827   5.388084   5.678747   5.340906   5.745226   5.803522
# ENSG00000186891   5.067202   5.272407   5.067202   5.067202   5.383952
# ENSG00000160072   7.874223   7.727235   7.403795   7.360128   7.297004
# ENSG00000041988   7.162992   6.882574   7.090246   7.198588   7.199933


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
                       time = samples$time[samples$gender == i],
                       age = samples$age[samples$gender == i]
  )
  
  ### pca plot rewite
  png(paste0("plot/pca_", gender_letter, ".png"))
  print(
    ggplot(pca_df, aes(x = PC1, y = PC2, color = time, shape = age)) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed() +
      ggtitle(paste0("PCA with VST data ", gender_letter))
  )
  dev.off()
  
  cat("\n")
}
# calculating PCA for dds male 
# PC1       PC2 group time       name
# SRR7007949 -13.356014  9.680155   pre  pre SRR7007949
# SRR7007950 -10.164726  6.537911  post post SRR7007950
# SRR7007951  -9.444039  3.697176   pre  pre SRR7007951
# SRR7007952  -8.871973  4.686951  post post SRR7007952
# SRR7007953  -7.102402 -3.803158   pre  pre SRR7007953
# SRR7007954  -9.305616  8.834624  post post SRR7007954
# 
# calculating PCA for dds female 
# PC1        PC2 group time       name
# SRR7007959 -1.986666 23.2517898   pre  pre SRR7007959
# SRR7007960 -2.794442 27.4865627  post post SRR7007960
# SRR7007971  6.818909 -6.4074718   pre  pre SRR7007971
# SRR7007972  4.612592  4.8911358  post post SRR7007972
# SRR7007973 -1.360130 -4.0128047   pre  pre SRR7007973
# SRR7007974  6.573547 -0.7215968  post post SRR7007974



#### DESeq2 ####
res = list()
for (i in gender_id) {
  gender_letter <- str_sub(i, 1, 1)
  print(i)
  dds[[i]]$time <- relevel(as.factor(dds[[i]]$time), ref = "pre")
  dds[[i]] <- DESeq(dds[[i]])
  res[[i]] <- results(dds[[i]])  # results(dds, contrast = c("time", "post", "pre"))
  print(head(res[[i]]))
  cat("\n")
  
  summary(res[[i]])  # FDR < 0.1  ex.
  # out of 47412 with nonzero total read count
  # adjusted p-value < 0.1
  # LFC > 0 (up)       : 2209, 4.7%
  # LFC < 0 (down)     : 2380, 5%
  # outliers [1]       : 0, 0%
  # low counts [2]     : 23900, 50%
  # (mean count < 3)
  # [1] see 'cooksCutoff' argument of ?results
  # [2] see 'independentFiltering' argument of ?results
  
  cat("\n")
}

#### results memo ####
### head(res[[male]]) ----
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

### summary(res[[male]]) ----
# out of 43164 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 842, 2%
# LFC < 0 (down)     : 1036, 2.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 24269, 56%
# (mean count < 10)


### head(res[[female]]) ----
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

### summary(res[[female]]) ----
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
  cat("strongest upregulation", "\n")
  print(head(ressig[order(ressig$log2FoldChange),]))
  cat("\n")
  
  # strongest upregulation
  cat("strongest upregulation", "\n")
  print(head(ressig[order(-ressig$log2FoldChange),]))
  cat("\n")
  
  ## histgram ##
  png(paste0("plot/hist_p_adjp_", i, ".png"), width = 1280, height = 800)
  par(mfrow = c(1,2))
  hist(results$pvalue, col = brewer.pal(3, name = "Set2")[1], 
       main = paste(i, "Pval"), xlab  = NULL)
  hist(results$padj, col = brewer.pal(3, name = "Set2")[2],
       main = paste(i, "adj.Pval"), xlab = NULL)
  dev.off()
}

## memo ##
# results for male 
# p < 0.05: 3860 
# p < 0.01: 1955 
# p < 0.001: 963 
# adj.P < 0.1: 1878 
# adj.P < 0.05: 1371 
# adj.P < 0.01: 762 
# 
# strongest upregulation 
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 6 rows and 6 columns
#                  baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#                  <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# ENSG00000227496  116.5837       -3.29135  0.374483  -8.78906 1.50819e-18 2.19209e-15
# ENSG00000174502   33.6753       -2.76794  0.420414  -6.58385 4.58415e-11 1.60403e-08
# ENSG00000286986   18.1806       -2.26824  0.384207  -5.90369 3.55456e-09 7.71992e-07
# ENSG00000129988   28.2626       -2.11128  0.319719  -6.60356 4.01410e-11 1.54788e-08
# ENSG00000106366   16.8828       -1.96804  0.268462  -7.33082 2.28756e-13 1.60656e-10
# ENSG00000041982   46.8292       -1.88652  0.356130  -5.29729 1.17531e-07 1.48050e-05
# 
# strongest upregulation 
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 6 rows and 6 columns
#                  baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#                  <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# ENSG00000125740   13.2386        2.63516  0.451009   5.84282 5.13255e-09 1.05412e-06
# ENSG00000170345  386.0794        1.87227  0.579366   3.23158 1.23106e-03 2.32376e-02
# ENSG00000136244   10.9666        1.79353  0.593338   3.02278 2.50465e-03 3.84134e-02
# ENSG00000120738  101.5929        1.56216  0.442272   3.53214 4.12213e-04 1.01284e-02
# ENSG00000248359   18.7395        1.55260  0.367012   4.23038 2.33300e-05 1.10205e-03
# ENSG00000170525 2484.2871        1.50181  0.273723   5.48660 4.09747e-08 6.04856e-06
# 
# results for female 
# p < 0.05: 4374 
# p < 0.01: 2281 
# p < 0.001: 1095 
# adj.P < 0.1: 2415 
# adj.P < 0.05: 1701 
# adj.P < 0.01: 875 
# 
# strongest upregulation 
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 6 rows and 6 columns
#                  baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#                  <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# ENSG00000227496   87.5117       -2.91255  0.350301  -8.31440 9.22145e-17 1.08666e-13
# ENSG00000129988   23.6534       -2.45880  0.524685  -4.68623 2.78284e-06 1.75677e-04
# ENSG00000006327   67.6641       -2.35358  0.410713  -5.73048 1.00147e-08 1.56655e-06
# ENSG00000286986   11.2858       -2.06700  0.353898  -5.84066 5.19942e-09 9.37805e-07
# ENSG00000119508  102.5864       -1.84395  0.566986  -3.25219 1.14518e-03 1.82526e-02
# ENSG00000167772   58.4354       -1.79976  0.583276  -3.08561 2.03137e-03 2.72027e-02
# 
# strongest upregulation 
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 6 rows and 6 columns
#                  baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#                  <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# ENSG00000205215   40.2066        1.81933  0.398505   4.56538 4.98593e-06 2.81570e-04
# ENSG00000205266   37.1814        1.78999  0.349877   5.11606 3.11980e-07 2.94896e-05
# ENSG00000226549   18.1802        1.78891  0.468998   3.81433 1.36555e-04 3.84354e-03
# ENSG00000170525 2680.8799        1.68915  0.314090   5.37793 7.53465e-08 8.32168e-06
# ENSG00000248461   19.5651        1.64191  0.325261   5.04798 4.46508e-07 3.90717e-05
# ENSG00000205312   55.0392        1.63692  0.375065   4.36437 1.27489e-05 5.94589e-04



#### Annotation & Write out the results ####
library("AnnotationDbi")
library("Homo.sapiens")
# columns(Homo.sapiens)
for (i in gender_id) {
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
  
  gender_letter <- str_sub(i, 1, 1)
  write.csv(results, paste0("resGSE113165", gender_letter, ".csv"))
}

## memo ##
# annotating male 
# 'select()' returned 1:many mapping between keys and columns
# head(results, 3): 
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 3 rows and 7 columns
#                 baseMean log2FoldChange     lfcSE      stat      pvalue
#                 <numeric>      <numeric> <numeric> <numeric>   <numeric>
# ENSG00000224051   258.402      -0.626479 0.1228123  -5.10111 3.37661e-07
# ENSG00000162494   149.988       0.337795 0.1303520   2.59141 9.55849e-03
# ENSG00000159423   538.822      -0.310330 0.0949877  -3.26706 1.08672e-03
#                     padj      symbol
#                   <numeric> <character>
# ENSG00000224051 0.000035643        CPTP
# ENSG00000162494 0.097678532      LRRC38
# ENSG00000159423 0.021256329     ALDH4A1
# 
# dim(results): 
# [1] 1878    7


# annotating female 
# 'select()' returned 1:many mapping between keys and columns
# head(results, 3): 
# log2 fold change (MLE): time post vs pre 
# Wald test p-value: time post vs pre 
# DataFrame with 3 rows and 7 columns
#                 baseMean log2FoldChange     lfcSE      stat     pvalue      padj
#                 <numeric>      <numeric> <numeric> <numeric>  <numeric> <numeric>
# ENSG00000169972   25.5158       0.554527  0.186821   2.96823 0.00299518 0.0360155
# ENSG00000162591   70.7050       0.368828  0.140824   2.61906 0.00881714 0.0750639
# ENSG00000142583  575.7197      -0.259946  0.104325  -2.49170 0.01271332 0.0955852
#                       symbol
#                   <character>
# ENSG00000169972       PUSL1
# ENSG00000162591       MEGF6
# ENSG00000142583      SLC2A5
# 
# dim(results): 
# [1] 2415    7



#### Volcano plot ####
library(EnhancedVolcano)
library(patchwork)

result_files <- list.files(pattern = ".csv")
result_files
for (i in 1:length(result_files)) {
  file <- result_files[i]
  results <- read.csv(file)
  name <- str_sub(file, 4, 13)
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

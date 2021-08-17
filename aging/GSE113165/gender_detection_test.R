#### gender detection ####
# identify subjects' gender by "XIST (ENSG00000229807)"
# importing count table (Rna-seq) ----
count_data <- "countsRNA_prerRnaRemoval.txt"
counts <- read.table(count_data, sep = "\t", header = TRUE)
View(counts)
rownames(counts) <- counts$Geneid
counts <- counts[,7:ncol(counts)]
dim(counts)  # 60664    90
# colnames(counts) <- sub(".fastp.fastq.sort.bam", "", colnames(counts))
colnames(counts) <- sub("X.Volumes.HDD24TB.MetaAnalysisProject_Apr2021.inactive.GSE162730.fastq_Rnaseq.star.", "", colnames(counts))
colnames(counts) <- sub(".Aligned.sortedByCoord.out.bam", "", colnames(counts))
head(counts)

# summary(count_rna["ENSG00000229807",])
# ncol(count_rna["ENSG00000229807",])

samples_rna <- filter(samples, library == "RNA-Seq")
keep <- match(samples_rna$run, colnames(counts))
counts <- counts[,keep]
counts[1:5, 1:5]
#                 SRR13202581 SRR13202584 SRR13202605 SRR13202608 SRR13202617
# ENSG00000284662        0.86        0.31        0.45        1.20        1.60
# ENSG00000186827       20.63       10.25       14.17        8.14        7.38
# ENSG00000186891        5.33        4.00        1.00        1.00        2.00
# ENSG00000160072      425.39      197.80      225.18      258.94      232.75
# ENSG00000041988      158.07       76.25       80.95       70.18      101.74

samples_rna$run[1:5]
# "SRR13202581" "SRR13202584" "SRR13202605" "SRR13202608" "SRR13202617"

dim(counts)  # 60664  16


# checking the expression of XIST (ENSG00000229807) ----
counts_xist <- counts %>%
    filter(rownames(.) == "ENSG00000229807") %>% 
    t()
counts_xist
#             ENSG00000229807
# SRR13202581        35559.82
# SRR13202584        17581.34
# SRR13202605        16094.29
# SRR13202608        14841.51
# SRR13202617        17185.43
# SRR13202620        13848.68
# SRR13202653         4276.99
# SRR13202683         7712.66
# SRR13202593            5.75
# SRR13202596            5.76
# SRR13202665         4371.83
# SRR13202677         7983.59
# SRR13202701         3128.71
# SRR13202713         3935.49
# SRR13202749           10.99
# SRR13202752            5.31

summary(counts_xist)
#  ENSG00000229807   
# Min.   :    5.31  
# 1st Qu.: 2349.28  
# Median : 6042.24  
# Mean   : 9159.26  
# 3rd Qu.:15154.70  
# Max.   :35559.82 

# extracting male samples ----
male <- counts_xist < 100
counts_male <- counts_xist[male,]
counts_male
# SRR13202593 SRR13202596 SRR13202749 SRR13202752 
# 5.75        5.76       10.99        5.31

head(samples_rna)
samples_rna_male <- samples_rna[samples_rna$run %in% names(counts_male),]
sample_table_male <- data.frame(ENSG00000229807 = counts_male, samples_rna_male)
sample_table_male
#             ENSG00000229807         run subject library time
# SRR13202593            5.75 SRR13202593     S23 RNA-Seq  pre
# SRR13202596            5.76 SRR13202596     S23 RNA-Seq post
# SRR13202749           10.99 SRR13202749     S20 RNA-Seq  pre
# SRR13202752            5.31 SRR13202752     S20 RNA-Seq post

# memo: from avobe, patient 20 & 23 are asuumed to be male subjects ----
# adding gender column to sample table ----
samples$gender <- ifelse(samples$subject %in% c("S20", "S23"), "male", "female")
samples
#             run subject library time gender
# 1   SRR13202575     S20   OTHER  pre   male
# 4   SRR13202581     S21 RNA-Seq  pre female
# 5   SRR13202584     S21 RNA-Seq post female
# 7   SRR13202587     S21   OTHER  pre female
# 15  SRR13202602     S23   OTHER post   male
# 16  SRR13202605     S25 RNA-Seq  pre female
# 18  SRR13202608     S25 RNA-Seq post female
# 20  SRR13202611     S25   OTHER  pre female
# 21  SRR13202614     S25   OTHER post female
# 22  SRR13202617     S26 RNA-Seq  pre female
# 23  SRR13202620     S26 RNA-Seq post female
# 25  SRR13202623     S26   OTHER  pre female
# 27  SRR13202626     S26   OTHER post female
# 41  SRR13202653     S17 RNA-Seq  pre female
# 55  SRR13202683     S18 RNA-Seq post female
# 94  SRR13202578     S20   OTHER post   male
# 99  SRR13202590     S21   OTHER post female
# 100 SRR13202593     S23 RNA-Seq  pre   male
# 101 SRR13202596     S23 RNA-Seq post   male
# 103 SRR13202599     S23   OTHER  pre   male
# 118 SRR13202629     S17   OTHER  pre female
# 124 SRR13202641     S17   OTHER post female
# 135 SRR13202665     S17 RNA-Seq post female
# 142 SRR13202677     S18 RNA-Seq  pre female
# 149 SRR13202689     S18   OTHER  pre female
# 152 SRR13202695     S18   OTHER post female
# 156 SRR13202701     S19 RNA-Seq  pre female
# 163 SRR13202713     S19 RNA-Seq post female
# 168 SRR13202725     S19   OTHER  pre female
# 173 SRR13202737     S19   OTHER post female
# 179 SRR13202749     S20 RNA-Seq  pre   male
# 180 SRR13202752     S20 RNA-Seq post   male

dim(samples)  # 32  5

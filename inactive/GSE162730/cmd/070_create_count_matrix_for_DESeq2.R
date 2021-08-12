# creating count matrix that can be used for DESeq2
setwd("/Volumes/HDD24TB/MetaAnalysisProject_Apr2021/inactive/GSE162730")

#### helpers ####
library("stringr")
library("dplyr")
library("tidyr")


#### Importing Data ####
#### Sample table ####
# tyding up sample table ----
samples <- read.table("SraRunTable.txt", sep = ",", header = TRUE)
head(samples)
samples <- select(samples, c("Run", "source_name", "Assay.Type", "Time_point"))
samples <- samples[,c("Run", "source_name", "Assay.Type", "Time_point")]
names(samples) <- c("run", "subject", "library", "time")
head(samples, 3)
#           run    subject library                                        time
# 1 SRR13202575 Patient 20   OTHER             Pre reduced activity 0 baseline
# 2 SRR13202579 Patient 20   OTHER  Post reduced activity 60 min after Leucine
# 3 SRR13202580 Patient 20   OTHER Post reduced activity 180 min after Leucine

# renaming subjects ----
unique(samples$subject)
# "Patient 20" "Patient 21" "Patient 23" "Patient 25"
# "Patient 26" "Patient 17" "Patient 18" "Patient 19"
samples$subject <- paste0("S", str_sub(samples$subject, -2, -1))

# removing duplicated samples ----
# choose the first run number from each subject & time
samples %>% 
    filter(subject == "S20")
#            run subject library                                        time
# 1  SRR13202575     S20   OTHER             Pre reduced activity 0 baseline
# 2  SRR13202579     S20   OTHER  Post reduced activity 60 min after Leucine
# 3  SRR13202580     S20   OTHER Post reduced activity 180 min after Leucine
# 4  SRR13202750     S20 RNA-Seq   Pre reduced activity 60 min after Leucine
# 5  SRR13202751     S20 RNA-Seq  Pre reduced activity 180 min after Leucine
# 6  SRR13202753     S20 RNA-Seq  Post reduced activity 60 min after Leucine
# 7  SRR13202754     S20 RNA-Seq Post reduced activity 180 min after Leucine
# 8  SRR13202576     S20   OTHER   Pre reduced activity 60 min after Leucine
# 9  SRR13202577     S20   OTHER  Pre reduced activity 180 min after Leucine
# 10 SRR13202578     S20   OTHER            Post reduced activity 0 baseline
# 11 SRR13202749     S20 RNA-Seq             Pre reduced activity 0 baseline
# 12 SRR13202752     S20 RNA-Seq            Post reduced activity 0 baseline

# memo: we need "time" containing "baseline" ----
samples %>% 
    filter(str_detect(time, "baseline")) %>% 
    select(subject) %>% 
    table()
# S17 S18 S19 S20 S21 S23 S25 S26 
#  16   8  16   4   4   4   4   4

# memo: S17, 18, 19 should be reduced to 4 ----
samples %>% 
    filter(str_detect(time, "baseline") & subject %in% c("S17", "S18", "S19")) %>% 
    arrange(subject, run)
### memo: use the first entry from each subject's pre/post
#            run subject library                             time
# 1  SRR13202629     S17   OTHER  Pre reduced activity 0 baseline  ## SRR13202629
# 2  SRR13202630     S17   OTHER  Pre reduced activity 0 baseline
# 3  SRR13202631     S17   OTHER  Pre reduced activity 0 baseline
# 4  SRR13202632     S17   OTHER  Pre reduced activity 0 baseline
# 5  SRR13202641     S17   OTHER Post reduced activity 0 baseline  ## SRR13202641
# 6  SRR13202642     S17   OTHER Post reduced activity 0 baseline
# 7  SRR13202643     S17   OTHER Post reduced activity 0 baseline
# 8  SRR13202644     S17   OTHER Post reduced activity 0 baseline
# 9  SRR13202653     S17 RNA-Seq  Pre reduced activity 0 baseline  ## SRR13202653
# 10 SRR13202654     S17 RNA-Seq  Pre reduced activity 0 baseline
# 11 SRR13202655     S17 RNA-Seq  Pre reduced activity 0 baseline
# 12 SRR13202656     S17 RNA-Seq  Pre reduced activity 0 baseline
# 13 SRR13202665     S17 RNA-Seq Post reduced activity 0 baseline  ## SRR13202665
# 14 SRR13202666     S17 RNA-Seq Post reduced activity 0 baseline
# 15 SRR13202667     S17 RNA-Seq Post reduced activity 0 baseline
# 16 SRR13202668     S17 RNA-Seq Post reduced activity 0 baseline
# 17 SRR13202677     S18 RNA-Seq  Pre reduced activity 0 baseline  ## SRR13202677
# 18 SRR13202678     S18 RNA-Seq  Pre reduced activity 0 baseline
# 19 SRR13202683     S18 RNA-Seq Post reduced activity 0 baseline  ## SRR13202683
# 20 SRR13202684     S18 RNA-Seq Post reduced activity 0 baseline
# 21 SRR13202689     S18   OTHER  Pre reduced activity 0 baseline  ## SRR13202689
# 22 SRR13202690     S18   OTHER  Pre reduced activity 0 baseline
# 23 SRR13202695     S18   OTHER Post reduced activity 0 baseline  ## SRR13202695
# 24 SRR13202696     S18   OTHER Post reduced activity 0 baseline
# 25 SRR13202701     S19 RNA-Seq  Pre reduced activity 0 baseline  ## SRR13202701
# 26 SRR13202702     S19 RNA-Seq  Pre reduced activity 0 baseline
# 27 SRR13202703     S19 RNA-Seq  Pre reduced activity 0 baseline
# 28 SRR13202704     S19 RNA-Seq  Pre reduced activity 0 baseline
# 29 SRR13202713     S19 RNA-Seq Post reduced activity 0 baseline  ## SRR13202713
# 30 SRR13202714     S19 RNA-Seq Post reduced activity 0 baseline
# 31 SRR13202715     S19 RNA-Seq Post reduced activity 0 baseline
# 32 SRR13202716     S19 RNA-Seq Post reduced activity 0 baseline
# 33 SRR13202725     S19   OTHER  Pre reduced activity 0 baseline  ## SRR13202725
# 34 SRR13202726     S19   OTHER  Pre reduced activity 0 baseline
# 35 SRR13202727     S19   OTHER  Pre reduced activity 0 baseline
# 36 SRR13202728     S19   OTHER  Pre reduced activity 0 baseline
# 37 SRR13202737     S19   OTHER Post reduced activity 0 baseline  ## SRR13202737
# 38 SRR13202738     S19   OTHER Post reduced activity 0 baseline
# 39 SRR13202739     S19   OTHER Post reduced activity 0 baseline
# 40 SRR13202740     S19   OTHER Post reduced activity 0 baseline

# extracting unique rows ----
samples_unique <- samples %>% 
    filter(str_detect(time, "baseline")) %>% 
    arrange(subject, run) %>% 
    group_by(subject, library, time) %>% 
    summarise(runs_to_use = head(run, 1), count = n_distinct(run)) %>% 
    arrange(subject, runs_to_use)
samples_unique
# A tibble: 32 × 5
# Groups:   subject, library [16]
#   subject library time                             runs_to_use count
#   <chr>   <chr>   <chr>                            <chr>       <int>
# 1 S17     OTHER   Pre reduced activity 0 baseline  SRR13202629     4
# 2 S17     OTHER   Post reduced activity 0 baseline SRR13202641     4
# 3 S17     RNA-Seq Pre reduced activity 0 baseline  SRR13202653     4
# 4 S17     RNA-Seq Post reduced activity 0 baseline SRR13202665     4
# 5 S18     RNA-Seq Pre reduced activity 0 baseline  SRR13202677     2
# 6 S18     RNA-Seq Post reduced activity 0 baseline SRR13202683     2
# 7 S18     OTHER   Pre reduced activity 0 baseline  SRR13202689     2
# 8 S18     OTHER   Post reduced activity 0 baseline SRR13202695     2
# 9 S19     RNA-Seq Pre reduced activity 0 baseline  SRR13202701     4
# 10 S19     RNA-Seq Post reduced activity 0 baseline SRR13202713     4
# … with 22 more rows

# removing duplicated rows from all-sample table ----
samples <- samples[samples$run %in% samples_unique$runs_to_use,]
samples
#             run subject library                             time
# 1   SRR13202575     S20   OTHER  Pre reduced activity 0 baseline
# 4   SRR13202581     S21 RNA-Seq  Pre reduced activity 0 baseline
# 5   SRR13202584     S21 RNA-Seq Post reduced activity 0 baseline
# 7   SRR13202587     S21   OTHER  Pre reduced activity 0 baseline
# 15  SRR13202602     S23   OTHER Post reduced activity 0 baseline
# 16  SRR13202605     S25 RNA-Seq  Pre reduced activity 0 baseline
# 18  SRR13202608     S25 RNA-Seq Post reduced activity 0 baseline
# 20  SRR13202611     S25   OTHER  Pre reduced activity 0 baseline
# 21  SRR13202614     S25   OTHER Post reduced activity 0 baseline
# 22  SRR13202617     S26 RNA-Seq  Pre reduced activity 0 baseline
# 23  SRR13202620     S26 RNA-Seq Post reduced activity 0 baseline
# 25  SRR13202623     S26   OTHER  Pre reduced activity 0 baseline
# 27  SRR13202626     S26   OTHER Post reduced activity 0 baseline
# 41  SRR13202653     S17 RNA-Seq  Pre reduced activity 0 baseline
# 55  SRR13202683     S18 RNA-Seq Post reduced activity 0 baseline
# 94  SRR13202578     S20   OTHER Post reduced activity 0 baseline
# 99  SRR13202590     S21   OTHER Post reduced activity 0 baseline
# 100 SRR13202593     S23 RNA-Seq  Pre reduced activity 0 baseline
# 101 SRR13202596     S23 RNA-Seq Post reduced activity 0 baseline
# 103 SRR13202599     S23   OTHER  Pre reduced activity 0 baseline
# 118 SRR13202629     S17   OTHER  Pre reduced activity 0 baseline
# 124 SRR13202641     S17   OTHER Post reduced activity 0 baseline
# 135 SRR13202665     S17 RNA-Seq Post reduced activity 0 baseline
# 142 SRR13202677     S18 RNA-Seq  Pre reduced activity 0 baseline
# 149 SRR13202689     S18   OTHER  Pre reduced activity 0 baseline
# 152 SRR13202695     S18   OTHER Post reduced activity 0 baseline
# 156 SRR13202701     S19 RNA-Seq  Pre reduced activity 0 baseline
# 163 SRR13202713     S19 RNA-Seq Post reduced activity 0 baseline
# 168 SRR13202725     S19   OTHER  Pre reduced activity 0 baseline
# 173 SRR13202737     S19   OTHER Post reduced activity 0 baseline
# 179 SRR13202749     S20 RNA-Seq  Pre reduced activity 0 baseline
# 180 SRR13202752     S20 RNA-Seq Post reduced activity 0 baseline

nrow(samples)  # 32
table(samples$subject)
# S17 S18 S19 S20 S21 S23 S25 S26 
#   4   4   4   4   4   4   4   4

with(samples, table(subject, library))
# library
# subject OTHER RNA-Seq
# S17     2       2
# S18     2       2
# S19     2       2
# S20     2       2
# S21     2       2
# S23     2       2
# S25     2       2
# S26     2       2

with(samples, table(subject, time))
# time
# subject Post reduced activity 0 baseline Pre reduced activity 0 baseline
# S17                                2                               2
# S18                                2                               2
# S19                                2                               2
# S20                                2                               2
# S21                                2                               2
# S23                                2                               2
# S25                                2                               2
# S26                                2                               2


# renaming time ----
samples$time <- ifelse(str_detect(samples$time, "Pre"), "pre", "post")
head(samples)
#            run subject library time
# 1  SRR13202575     S20   OTHER  pre
# 4  SRR13202581     S21 RNA-Seq  pre
# 5  SRR13202584     S21 RNA-Seq post
# 7  SRR13202587     S21   OTHER  pre
# 15 SRR13202602     S23   OTHER post
# 16 SRR13202605     S25 RNA-Seq  pre



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

# writing out sample tables ----
samples <- samples %>% arrange(run)
write.table(samples, "sampletable_all.txt", sep = "\t")

samples_rna <- filter(samples, library == "RNA-Seq")
samples_ribo <- filter(samples, library == "OTHER")
write.table(samples_rna, "sampletable_rna.txt", sep = "\t")
write.table(samples_ribo, "sampletable_ribo.txt", sep = "\t")


# creating SRR_Acc_List containing avobe samples run number ----
write.table(samples_rna$run, "SRR_Acc_List_rna_unique.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(samples_ribo$run, "SRR_Acc_List_ribo_unique.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


# writing out count matrix ----
counts <- counts[,match(samples_rna$run, colnames(counts))]
write.table(counts, "countmat_rna.txt", sep = "\t", col.names = TRUE)


# checking the data ----
s <- read.table("sampletable_rna.txt", sep = "\t", header = TRUE)
head(s)
#           run subject library time gender
# 1 SRR13202581     S21 RNA-Seq  pre female
# 2 SRR13202584     S21 RNA-Seq post female
# 3 SRR13202593     S23 RNA-Seq  pre   male
# 4 SRR13202596     S23 RNA-Seq post   male
# 5 SRR13202605     S25 RNA-Seq  pre female
# 6 SRR13202608     S25 RNA-Seq post female
dim(s)  # 16  5

c <- read.table("countmat_rna.txt", sep = "\t", header = TRUE)
c[1:5,1:6]
#                 SRR13202581 SRR13202584 SRR13202593 SRR13202596 SRR13202605 SRR13202608
# ENSG00000284662        0.86        0.31        0.25        0.34        0.45        1.20
# ENSG00000186827       20.63       10.25        5.17        5.27       14.17        8.14
# ENSG00000186891        5.33        4.00        0.00        2.03        1.00        1.00
# ENSG00000160072      425.39      197.80      129.01      165.94      225.18      258.94
# ENSG00000041988      158.07       76.25       56.45       61.71       80.95       70.18
dim(c)  # 60664    16

# creating count matrix that can be used for DESeq2
setwd("/Volumes/HDD24TB/MetaAnalysisProject_Apr2021/inactive/GSE162730")

#### helpers ####
library("stringr")
library("dplyr")
library("tidyr")

#### Importing Data ####
counts <- read.table("counts.txt", sep = "\t", header = TRUE)
head(counts)
View(counts)
rownames(counts) <- counts$Geneid
colnames(counts) <- sub(".fastp.fastq.sort.bam", "", colnames(counts))

# samples <- read.table("pdata.txt", sep = "\t", header = TRUE)
samples <- read.table("SraRunTable.txt", sep = ",", header = TRUE)
head(samples)
samples <- select(samples, c("Run", "source_name", "Assay.Type", "Time_point"))
samples <- samples[,c("Run", "source_name", "Assay.Type", "Time_point")]
names(samples) <- c("run", "subject", "library", "time")
head(samples, 3)
#           run    subject library                                        time gender
# 1 SRR13202575 Patient 20   OTHER             Pre reduced activity 0 baseline   male
# 2 SRR13202579 Patient 20   OTHER  Post reduced activity 60 min after Leucine   male
# 3 SRR13202580 Patient 20   OTHER Post reduced activity 180 min after Leucine   male


# subsetting: Rna-seq and Ribo-seq
unique(samples$assay.type)  # "OTHER"   "RNA-Seq"
samples_rna <- subset(samples, library == "RNA-Seq")
samples_ribo <- subset(samples, library == "OTHER")
head(samples_rna, 3)
#           run source_name assay.type                                  time_point
# 4 SRR13202581  Patient 21    RNA-Seq             Pre reduced activity 0 baseline
# 5 SRR13202584  Patient 21    RNA-Seq            Post reduced activity 0 baseline
# 6 SRR13202586  Patient 21    RNA-Seq Post reduced activity 180 min after Leucine

head(samples_ribo, 3)
#           run source_name assay.type                                  time_point
# 1 SRR13202575  Patient 20      OTHER             Pre reduced activity 0 baseline
# 2 SRR13202579  Patient 20      OTHER  Post reduced activity 60 min after Leucine
# 3 SRR13202580  Patient 20      OTHER Post reduced activity 180 min after Leucine



### tring if XIST could identify gender
# Rna-seq
count_rna <- read.table("countsRNA_prerRnaRemoval.txt", sep = "\t", header = TRUE)
head(count_rna)
rownames(count_rna) <- count_rna$Geneid
colnames(count_rna) <- sub("X.Volumes.HDD24TB.MetaAnalysisProject_Apr2021.inactive.GSE162730.fastq_Rnaseq.star.", "", colnames(count_rna))
colnames(count_rna) <- sub(".Aligned.sortedByCoord.out.bam", "", colnames(count_rna))
names(count_rna)
# summary(count_rna["ENSG00000229807",])
# ncol(count_rna["ENSG00000229807",])

count_rna <- count_rna[,7:ncol(count_rna)]
dim(count_rna)  # 60664    90

keep <- match(samples_rna$run, colnames(count_rna))
count_rna <- count_rna[,keep]
head(count_rna)
head(samples_rna)


# tiding up sample information
head(samples_rna)
#           run     subject    library                                        time
# 1 SRR13202575  Patient 20      OTHER             Pre reduced activity 0 baseline
# 2 SRR13202579  Patient 20      OTHER  Post reduced activity 60 min after Leucine
# 3 SRR13202580  Patient 20      OTHER Post reduced activity 180 min after Leucine
# 4 SRR13202581  Patient 21    RNA-Seq             Pre reduced activity 0 baseline
# 5 SRR13202584  Patient 21    RNA-Seq            Post reduced activity 0 baseline
# 6 SRR13202586  Patient 21    RNA-Seq Post reduced activity 180 min after Leucine


unique(samples_rna$subject)
# "Patient 20" "Patient 21" "Patient 23" "Patient 25"
# "Patient 26" "Patient 17" "Patient 18" "Patient 19"
samples_rna$subject <- paste0("S", str_sub(samples_rna$subject, -2, -1))


#### removing duplicates ####
# chose the first run number from each subject & time
samples_rna %>% 
    filter(subject == "S20" & library == "RNA-Seq")
#           run subject library                                        time
# 1 SRR13202750     S20 RNA-Seq   Pre reduced activity 60 min after Leucine
# 2 SRR13202751     S20 RNA-Seq  Pre reduced activity 180 min after Leucine
# 3 SRR13202753     S20 RNA-Seq  Post reduced activity 60 min after Leucine
# 4 SRR13202754     S20 RNA-Seq Post reduced activity 180 min after Leucine
# 5 SRR13202749     S20 RNA-Seq             Pre reduced activity 0 baseline
# 6 SRR13202752     S20 RNA-Seq            Post reduced activity 0 baseline

### memo: we need "time" containing "baseline"
samples_rna %>% 
    filter(library == "RNA-Seq" & str_detect(time, "baseline")) %>% 
    select(subject) %>% 
    table()
# S17 S18 S19 S20 S21 S23 S25 S26 
# 8   4   8   2   2   2   2   2

### memo: S17, 18, 19 should be reduced to 2
samples_rna %>% 
    filter(library == "RNA-Seq" & str_detect(time, "baseline") & subject %in% c("S17", "S18", "S19")) %>% 
    arrange(subject, run)
### memo: use the first entry from each subject's pre/post
#            run subject library                             time
# 1  SRR13202653     S17 RNA-Seq  Pre reduced activity 0 baseline  ##
# 2  SRR13202654     S17 RNA-Seq  Pre reduced activity 0 baseline
# 3  SRR13202655     S17 RNA-Seq  Pre reduced activity 0 baseline
# 4  SRR13202656     S17 RNA-Seq  Pre reduced activity 0 baseline
# 5  SRR13202665     S17 RNA-Seq Post reduced activity 0 baseline  ##
# 6  SRR13202666     S17 RNA-Seq Post reduced activity 0 baseline
# 7  SRR13202667     S17 RNA-Seq Post reduced activity 0 baseline
# 8  SRR13202668     S17 RNA-Seq Post reduced activity 0 baseline
# 9  SRR13202677     S18 RNA-Seq  Pre reduced activity 0 baseline  ##
# 10 SRR13202678     S18 RNA-Seq  Pre reduced activity 0 baseline
# 11 SRR13202683     S18 RNA-Seq Post reduced activity 0 baseline  ##
# 12 SRR13202684     S18 RNA-Seq Post reduced activity 0 baseline
# 13 SRR13202701     S19 RNA-Seq  Pre reduced activity 0 baseline  ##
# 14 SRR13202702     S19 RNA-Seq  Pre reduced activity 0 baseline
# 15 SRR13202703     S19 RNA-Seq  Pre reduced activity 0 baseline
# 16 SRR13202704     S19 RNA-Seq  Pre reduced activity 0 baseline
# 17 SRR13202713     S19 RNA-Seq Post reduced activity 0 baseline  ##
# 18 SRR13202714     S19 RNA-Seq Post reduced activity 0 baseline
# 19 SRR13202715     S19 RNA-Seq Post reduced activity 0 baseline
# 20 SRR13202716     S19 RNA-Seq Post reduced activity 0 baseline


samples_rna_rm_dup <- samples_rna %>% 
    filter(library == "RNA-Seq" & str_detect(time, "baseline")) %>% 
    arrange(subject, run) %>% 
    group_by(subject, time) %>% 
    summarise(runs_to_use = head(run, 1), count = n_distinct(run)) %>% 
    arrange(subject, runs_to_use)
samples_rna_rm_dup
# A tibble: 16 × 4
# Groups:   subject [8]
# subject   time                             runs_to_use count
# <chr>     <chr>                            <chr>       <int>
# 1 S17     Pre reduced activity 0 baseline  SRR13202653     4
# 2 S17     Post reduced activity 0 baseline SRR13202665     4
# 3 S18     Pre reduced activity 0 baseline  SRR13202677     2
# 4 S18     Post reduced activity 0 baseline SRR13202683     2
# 5 S19     Pre reduced activity 0 baseline  SRR13202701     4
# 6 S19     Post reduced activity 0 baseline SRR13202713     4
# 7 S20     Pre reduced activity 0 baseline  SRR13202749     1
# 8 S20     Post reduced activity 0 baseline SRR13202752     1
# 9 S21     Pre reduced activity 0 baseline  SRR13202581     1
# 10 S21     Post reduced activity 0 baseline SRR13202584     1
# 11 S23     Pre reduced activity 0 baseline  SRR13202593     1
# 12 S23     Post reduced activity 0 baseline SRR13202596     1
# 13 S25     Pre reduced activity 0 baseline  SRR13202605     1
# 14 S25     Post reduced activity 0 baseline SRR13202608     1
# 15 S26     Pre reduced activity 0 baseline  SRR13202617     1
# 16 S26     Post reduced activity 0 baseline SRR13202620     1


samples_rna <- samples_rna[samples_rna$run %in% samples_rna_rm_dup$runs_to_use,]
samples_rna
#             run subject library                             time
# 4   SRR13202581     S21 RNA-Seq  Pre reduced activity 0 baseline
# 5   SRR13202584     S21 RNA-Seq Post reduced activity 0 baseline
# 16  SRR13202605     S25 RNA-Seq  Pre reduced activity 0 baseline
# 18  SRR13202608     S25 RNA-Seq Post reduced activity 0 baseline
# 22  SRR13202617     S26 RNA-Seq  Pre reduced activity 0 baseline
# 23  SRR13202620     S26 RNA-Seq Post reduced activity 0 baseline
# 41  SRR13202653     S17 RNA-Seq  Pre reduced activity 0 baseline
# 55  SRR13202683     S18 RNA-Seq Post reduced activity 0 baseline
# 100 SRR13202593     S23 RNA-Seq  Pre reduced activity 0 baseline
# 101 SRR13202596     S23 RNA-Seq Post reduced activity 0 baseline
# 135 SRR13202665     S17 RNA-Seq Post reduced activity 0 baseline
# 142 SRR13202677     S18 RNA-Seq  Pre reduced activity 0 baseline
# 156 SRR13202701     S19 RNA-Seq  Pre reduced activity 0 baseline
# 163 SRR13202713     S19 RNA-Seq Post reduced activity 0 baseline
# 179 SRR13202749     S20 RNA-Seq  Pre reduced activity 0 baseline
# 180 SRR13202752     S20 RNA-Seq Post reduced activity 0 baseline

dim(samples_rna)  # 16  4

samples_rna$time <- ifelse(str_detect(samples_rna$time, "Pre"), "pre", "post")
samples_rna$time <- relevel(as.factor(samples_rna$time), ref = "pre")


#### gender detection ####
count_rna_new <- count_rna %>%
    filter(rownames(.) == "ENSG00000229807") %>% 
    t()
head(count_rna_new)
#             ENSG00000229807
# SRR13202581        35559.82
# SRR13202584        17581.34
# SRR13202586        13726.25
# SRR13202594            4.46
# SRR13202595            5.79
# SRR13202598            1.26

summary(count_rna_new)
# ENSG00000229807   
# Min.   :    1.26  
# 1st Qu.: 3153.09  
# Median : 3941.03  
# Mean   : 6533.63  
# 3rd Qu.: 9726.54  
# Max.   :35559.82 

male <- count_rna_new < 100
count_rna_male <- count_rna_new[male,]
count_rna_male
length(count_rna_male)  # 12

names(count_rna_male)
head(samples_rna_rna)
samples_rna_male <- samples_rna[samples_rna$run %in% names(count_rna_male),]
sample_table_male <- data.frame(ENSG00000229807 = count_rna_male, samples_rna_male)
sample_table_male
#             ENSG00000229807         run     subject    library                                        time
# SRR13202594            4.46 SRR13202594  Patient 23    RNA-Seq   Pre reduced activity 60 min after Leucine
# SRR13202595            5.79 SRR13202595  Patient 23    RNA-Seq  Pre reduced activity 180 min after Leucine
# SRR13202598            1.26 SRR13202598  Patient 23    RNA-Seq Post reduced activity 180 min after Leucine
# SRR13202750            8.23 SRR13202750  Patient 20    RNA-Seq   Pre reduced activity 60 min after Leucine
# SRR13202751            9.89 SRR13202751  Patient 20    RNA-Seq  Pre reduced activity 180 min after Leucine
# SRR13202753            4.63 SRR13202753  Patient 20    RNA-Seq  Post reduced activity 60 min after Leucine
# SRR13202754            3.69 SRR13202754  Patient 20    RNA-Seq Post reduced activity 180 min after Leucine
# SRR13202593            5.75 SRR13202593  Patient 23    RNA-Seq             Pre reduced activity 0 baseline
# SRR13202596            5.76 SRR13202596  Patient 23    RNA-Seq            Post reduced activity 0 baseline
# SRR13202597            4.12 SRR13202597  Patient 23    RNA-Seq  Post reduced activity 60 min after Leucine
# SRR13202749           10.99 SRR13202749  Patient 20    RNA-Seq             Pre reduced activity 0 baseline
# SRR13202752            5.31 SRR13202752  Patient 20    RNA-Seq            Post reduced activity 0 baseline

#### assuming that patient 20 & 23 would be male subjects
# storing male data as a csv file
write.csv(sample_table_male, "sample_table_male.csv", row.names = TRUE)

samples_rna$gender <- ifelse(samples_rna$subject %in% c("S20", "S23"), "male", "female")
head(samples_rna)
#             run subject library time gender
# 4   SRR13202581     S21 RNA-Seq  pre female
# 5   SRR13202584     S21 RNA-Seq post female
# 16  SRR13202605     S25 RNA-Seq  pre female
# 18  SRR13202608     S25 RNA-Seq post female
# 22  SRR13202617     S26 RNA-Seq  pre female
# 23  SRR13202620     S26 RNA-Seq post female
# 41  SRR13202653     S17 RNA-Seq  pre female
# 55  SRR13202683     S18 RNA-Seq post female
# 100 SRR13202593     S23 RNA-Seq  pre   male
# 101 SRR13202596     S23 RNA-Seq post   male
# 135 SRR13202665     S17 RNA-Seq post female
# 142 SRR13202677     S18 RNA-Seq  pre female
# 156 SRR13202701     S19 RNA-Seq  pre female
# 163 SRR13202713     S19 RNA-Seq post female
# 179 SRR13202749     S20 RNA-Seq  pre   male
# 180 SRR13202752     S20 RNA-Seq post   male

dim(samples_rna)  # 16  5

write.csv(samples_rna, "sampletable_rna.csv", row.names = samples_rna$run)


# writing out tables
counts_rna <- counts_rna[,match(samples_rna$run, colnames(counts))]
# write.table(samples_rna, "sampletable.txt", sep = "\t", col.names = TRUE, row.names = samples_rna$run)  # already there
write.table(counts_rna, "countmat_rna.txt", sep = "\t", col.names = TRUE)

# checking the data
c <- read.table("countmat_rna.txt", sep = "\t", header = TRUE)
head(c)
dim(c)  # 60664    56

s <- read.table("sampletable_rna.txt", sep = "\t", header = TRUE)
head(s)
dim(s)  # 56  8

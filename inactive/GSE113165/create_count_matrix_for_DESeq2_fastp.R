# creating count matrix that can be used for DESeq2
setwd("/Volumes/HDD24TB/MetaAnalysisProject_Apr2021/inactive/GSE113165")

#### helpers ####
library("stringr")
library("dplyr")
library("tidyr")

#### Importing Data ####
# count matrix ----
counts <- read.table("counts.fastp.txt", sep = "\t", header = TRUE)
head(counts)
View(counts)
rownames(counts) <- counts$Geneid
counts <- counts[,7:ncol(counts)]
colnames(counts) <- sub(".fastp.fastq.sort.bam", "", colnames(counts))
counts[1:5,1:5]
#                 SRR7007949 SRR7007950 SRR7007951 SRR7007952 SRR7007953
# ENSG00000284662          9          4          2          4         12
# ENSG00000186827          3         23         14         11          7
# ENSG00000186891          1          2          0          0          1
# ENSG00000160072        154        168        153        198        195
# ENSG00000041988         88         78         99        132        115


# sample table ----
samples <- read.table("SraRunTable.txt", sep = ",", header = TRUE)
samples <- samples[,c("Run", "Subject", "sex", "age", "Time")]
names(samples) <- tolower(names(samples))
names(samples)[3] <- "gender"
samples$subject <- paste0("S", samples$subject)
samples$time <- ifelse(str_detect(samples$time, "pre"), "pre", "post")
head(samples, 3)
#          run subject gender age time
# 1 SRR7007949   S3005   male old  pre
# 2 SRR7007950   S3005   male old post
# 3 SRR7007951   S3008   male old  pre

table(samples$gender, samples$time)
#         post pre
# female   15  15
# male     13  13

table(samples$gender, samples$age)
#         old young
# female  16    14
# male    22     4



#### writing out tables ####
counts <- counts[,match(samples$run, colnames(counts))]
write.table(samples, "sampletable.txt", sep = "\t", col.names = TRUE, row.names = samples$run)
write.table(counts, "countmat.fastp.txt", sep = "\t", col.names = TRUE)

#### checking the data ####
s <- read.table("sampletable.txt", sep = "\t", header = TRUE)
head(s)
#                   run subject gender age time
# SRR7007949 SRR7007949   S3005   male old  pre
# SRR7007950 SRR7007950   S3005   male old post
# SRR7007951 SRR7007951   S3008   male old  pre
# SRR7007952 SRR7007952   S3008   male old post
# SRR7007953 SRR7007953   S3010   male old  pre
# SRR7007954 SRR7007954   S3010   male old post
dim(s)  # 56  5

c <- read.table("countmat.fastp.txt", sep = "\t", header = TRUE)
c[1:5, 1:5]
#                 SRR7007949 SRR7007950 SRR7007951 SRR7007952 SRR7007953
# ENSG00000284662          9          4          2          4         12
# ENSG00000186827          3         23         14         11          7
# ENSG00000186891          1          2          0          0          1
# ENSG00000160072        154        168        153        198        195
# ENSG00000041988         88         78         99        132        115
dim(c)  # 60664    56

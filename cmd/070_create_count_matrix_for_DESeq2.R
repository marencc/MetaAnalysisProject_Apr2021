# creating count matrix that can be used for DESeq2
setwd("/Volumes/HDD24TB/MetaAnalysisProject_Apr2021/inactive/GSE113165")

#### helpers ####
library("stringr")
library("dplyr")
library("tidyr")

#### Importing Data ####
counts <- read.table("counts.fastp.txt", sep = "\t", header = TRUE)
head(counts)
View(counts)
rownames(counts) <- counts$Geneid
colnames(counts) <- sub(".fastp.fastq.sort.bam", "", colnames(counts))

# samples <- read.table("pdata.txt", sep = "\t", header = TRUE)
samples <- read.table("SraRunTable.txt", sep = ",", header = TRUE)
samples <- samples[,c("Run", "age", "Sample.Name", "sex", "Subject", "susceptibility", "Time")]
names(samples) <- tolower(names(samples))
head(samples, 3)
#          run age sample.name  sex subject susceptibility          time
# 1 SRR7007949 old  GSM3098323 male    3005            low  pre bed rest
# 2 SRR7007950 old  GSM3098324 male    3005            low post bed rest
# 3 SRR7007951 old  GSM3098325 male    3008            low  pre bed rest

# samples$Time <- ifelse(str_detect(samples$Time, "pre"), "pre", "post")

# tiding up sample information
samples$age <- relevel(as.factor(samples$age), ref = "young")
samples$sex <- str_sub(samples$sex, 1, 1)
samples$sex <- relevel(as.factor(samples$sex), ref = "m")
samples$subject <- paste0("S", samples$subject)
samples$susceptibility <- str_replace(samples$susceptibility, "Not used in Susceptibility Study", "NA")
samples$susceptibility <- as.factor(samples$susceptibility)
samples$time <- ifelse(str_detect(samples$time, "pre"), "pre", "post")
samples$time <- relevel(as.factor(samples$time), ref = "pre")
summary(samples)
# run                age        sample.name        sex    subject            susceptibility time   
# Length:56          young:18   Length:56          m:26   Length:56          high:24        pre :28  
# Class :character   old  :38   Class :character   f:30   Class :character   low :28        post:28  
# Mode  :character              Mode  :character          Mode  :character   NA  : 4 

head(samples, 3)
#          run age sample.name sex subject susceptibility time
# 1 SRR7007949 old  GSM3098323   m   S3005            low  pre
# 2 SRR7007950 old  GSM3098324   m   S3005            low post
# 3 SRR7007951 old  GSM3098325   m   S3008            low  pre

table(samples$sex, samples$time)
#   pre post
# m  13   13
# f  15   15

col.time <- c("pre"="blue","post"="orange")
samples$color <- col.time[as.vector(samples$time)]
head(samples)


# writing out tables
counts <- counts[,match(samples$run, colnames(counts))]
# write.table(samples, "sampletable.txt", sep = "\t", col.names = TRUE, row.names = samples$run)  # already there
write.table(counts, "countmat.fastp.txt", sep = "\t", col.names = TRUE)

# checking the data
c <- read.table("countmat.fastp.txt", sep = "\t", header = TRUE)
head(c)
dim(c)  # 60664    56

s <- read.table("sampletable.txt", sep = "\t", header = TRUE)
head(s)
dim(s)  # 56  8

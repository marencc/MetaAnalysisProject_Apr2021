# Get the annotation GPL id (see Annotation: GPL10558)
gpl <- getGEO('GPL20880', destdir=".")
Meta(gpl)$title
# "[HuGene-1_1-st] Affymetrix Human Gene 1.1 ST Array [HuGene11st_Hs_ENTREZG_18.0.0]"

# Inspect the table of the gpl annotation object
colnames(Table(gpl))
# [1] "ID"          "Description" "SPOT_ID"  

# Get the gene symbol and entrez ids to be used for annotations
Table(gpl)[1:10,]
#              ID                                           Description   SPOT_ID
# 1          1_at                                alpha-1-B glycoprotein         1
# 2         10_at N-acetyltransferase 2 (arylamine N-acetyltransferase)        10
# 3        100_at                                   adenosine deaminase       100
# 4       1000_at             cadherin 2, type 1, N-cadherin (neuronal)      1000
# 5      10000_at         v-akt murine thymoma viral oncogene homolog 3     10000
# 6  100009676_at                                ZBTB11 antisense RNA 1 100009676
# 7      10001_at                            mediator complex subunit 6     10001
# 8      10002_at       nuclear receptor subfamily 2, group E, member 3     10002
# 9      10003_at        N-acetylated alpha-linked acidic dipeptidase 2     10003
# 10 100033413_at                    small nucleolar RNA, C/D box 116-1 100033413

dim(Table(gpl))
# [1] 19654     3

# Get the gene expression data for all the probes with a gene symbol
# geneProbes <- which(!is.na(Table(gpl)$Symbol))
# probeids <- as.character(Table(gpl)$ID[geneProbes])

probes <- intersect(probeids, rownames(exprs(eset)))
length(probes)

geneMatrix <- exprs(eset)[probes, ]

inds <- which(Table(gpl)$ID %in% probes)
# Check you get the same probes
head(probes)
head(as.character(Table(gpl)$ID[inds]))

# Create the expression matrix with gene ids
geneMatTable <- cbind(geneMatrix, Table(gpl)[inds, c(1, 6, 9, 12)])
head(geneMatTable)

# Save a copy of the expression matrix as a csv file
write.csv(geneMatTable, paste(GEO_DATASETS[1], "_DataMatrix.csv", sep=""), row.names=T)

#Laura Poulton
#MSCS REU
#6/21/17

#assignment - find patients in TCGA from "TCGA-KIRP" cancer type that have:
# gene expression, miRNA expression, DNA methylation, AND copy number variation
# data available

library(TCGAbiolinks)
library(dplyr)
library(DT)

#lists data.category, data.type, platform, and workflow types
datatable(readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=2046985454"), filter = 'top', options = list(scrollX = TRUE, keys = TRUE, pageLength = 20), rownames = FALSE)

#getting DNA methylation info
#Ilumina human methylation 450 platform
query.met.1 <- GDCquery(project = "TCGA-KIRP", data.category = "DNA methylation", legacy = TRUE, data.type="Methylation beta value", platform = "Illumina Human Methylation 450")

#getting DNA methylation info
#Ilumina human methylation 27 platform
query.met.2 <- GDCquery(project = "TCGA-KIRP", data.category = "DNA methylation", legacy = TRUE, data.type="Methylation beta value", platform = "Illumina Human Methylation 27")

#getting gene expression data
#HTSeq-Counts workflow
query.exp.1 <- GDCquery(project = "TCGA-KIRP", data.category="Transcriptome Profiling", data.type = "Gene Expression Quantification", legacy = FALSE, workflow.type = "HTSeq - Counts")

#getting gene expression data
#HTSeq-FPKM workflow
query.exp.2 <- GDCquery(project = "TCGA-KIRP", data.category="Transcriptome Profiling", data.type ="Gene Expression Quantification", legacy = FALSE, workflow.type ="HTSeq - FPKM")

#getting gene expression data
#HTSeq-FPKM-UQ workflow
query.exp.3 <- GDCquery(project = "TCGA-KIRP", data.category="Transcriptome Profiling", data.type ="Gene Expression Quantification", legacy = FALSE, workflow.type ="HTSeq - FPKM-UQ")

#getting miRNA expression data
query.miR <- GDCquery(project = "TCGA-KIRP", data.category="Transcriptome Profiling", data.type = "miRNA Expression Quantification", legacy = FALSE)

#getting copy number variation data
#Copy number segment type
query.num.1 <- GDCquery(project = "TCGA-KIRP", data.category="Copy Number Variation", data.type = "Copy Number Segment", legacy = FALSE)

#getting copy number variation data
#masked Copy number segment type
query.num.2 <- GDCquery(project = "TCGA-KIRP", data.category="Copy Number Variation", data.type = "Masked Copy Number Segment", legacy = FALSE)

#now to find patients with all types of data
#note that the cases column contains the barcode
#and the patient ID is the first 12 characters of the barcode
#third section separated by dashes
#also note that intersect only takes two arguments at a time

#patients with dna methylation and gene expression data
met.exp <- intersect(substr(getResults(c(query.met.1, query.met.2), cols = "cases"), 1, 12), substr(getResults(c(query.exp.1, query.exp.2, query.exp.3), cols = "cases"), 1, 12))

#patients with miRNA expression and copy number variation data
miR.num <- intersect(substr(getResults(query.miR, cols = "cases"), 1, 12), substr(getResults(c(query.num.1, query.num.2), cols = "cases"), 1, 12))

#patients with all
common.patients <- intersect(met.exp, miR.num)

#second part of assignment
#find number of normal samples
#and solid tumor samples
#that satisfy the above conditions

#finding dna methylation solid tumor samples
tumors.met <- GDCquery(project= "TCGA-KIRP", data.category= "DNA Methylation", barcode = common.patients, sample.type = c("Primary solid Tumor", "Recurrent Solid Tumor"))

#finding gene expression solid tumor samples
tumors.exp <- GDCquery(project= "TCGA-KIRP", data.category= "Transcriptome Profiling", barcode = common.patients, sample.type = c("Primary solid Tumor", "Recurrent Solid Tumor"))

#finding copy number variation solid tumor samples
tumors.num <- GDCquery(project= "TCGA-KIRP", data.category= "Copy Number Variation", barcode = common.patients, sample.type = c("Primary solid Tumor", "Recurrent Solid Tumor"))

#finding dna methylation normal samples
normal.met <- GDCquery(project= "TCGA-KIRP", data.category= "DNA Methylation", barcode = common.patients, sample.type = c("Blood Derived Normal", "Solid Tissue Normal", "Buccal Cell Normal", "EBV Immortalized Normal", "Bone Marrow Normal"))

#finding gene expression normal samples
normal.exp <- GDCquery(project= "TCGA-KIRP", data.category= "Transcriptome Profiling", barcode = common.patients, sample.type = c("Blood Derived Normal", "Solid Tissue Normal", "Buccal Cell Normal", "EBV Immortalized Normal", "Bone Marrow Normal"))

#finding copy number variation normal samples
normal.num <- GDCquery(project= "TCGA-KIRP", data.category= "Copy Number Variation", barcode = common.patients, sample.type = c("Blood Derived Normal", "Solid Tissue Normal", "Buccal Cell Normal", "EBV Immortalized Normal", "Bone Marrow Normal"))


#now finding intersections
tumors.exp.met <- intersect(substr(getResults(tumors.met, cols = "cases"), 1, 15), substr(getResults(tumors.exp, cols = "cases"), 1, 15))

tumors <- intersect(substr(getResults(tumors.num, cols = "cases"), 1, 15), tumors.exp.met)

#now for normal samples
normal.exp.met <- intersect(substr(getResults(normal.met, cols = "cases"), 1, 15), substr(getResults(normal.exp, cols = "cases"), 1, 15)) 

normal <- intersect(substr(getResults(normal.num, cols = "cases"), 1, 15), normal.exp.met)

#results: 274 tumor samples and 22 normal samples

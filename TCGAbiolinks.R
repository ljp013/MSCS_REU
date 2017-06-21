#Laura Poulton
#MSCS REU
#6/21/17

#TCGA biolinks tutorial
#from: https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/index.html


#to install
#source("https://bioconductor.org/biocLite.R")
#biocLite("TCGAbiolinks")

#libraries needed
#library(TCGAbiolinks)
#library(dplyr)
#library(DT) #may need to install using biocLite("DT") first


main <- function() {
    # #creates data table of TCGA data options 
    # #(copy number variation, gene expression, miRNA expression, etc)
    # #all columns in this data table can be arguments in a GDC query
    # datatable(readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=2046985454"), filter = 'top', options = list(scrollX = TRUE, keys = TRUE, pageLength = 20), rownames = FALSE)
    # 
    # #creates data table of DNA methylation data available for glioblastoma and low grade glioma samples
    # #TCGA-GBM - glioblastoma multiform
    # #TCGA-LGG - low grade gliomas
    # #legacy = FALSE means we are searching the harmonized (current) database
    # #GDC is Genomic Data Commons
    # #Illumnia Human Methylation 450 is one type of DNA methylation data available
    # 
    # #will search GDC database for query 
    # 
    # query <- GDCquery(project = c("TCGA-GBM", "TCGA-LGG"),
    #                   data.category = "DNA Methylation",
    #                   legacy = FALSE,
    #                   platform = c("Illumina Human Methylation 450"),
    #                   sample.type = "Recurrent Solid Tumor")
    # 
    # #and then create data table from results
    # #filter = 'top' means put column labels at the top of the table
    # 
    # datatable(getResults(query), 
    #           filter = 'top',
    #           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    #           rownames = FALSE)
    # 
    # 
    # #find gene expression and methylation data on Colon Adenocarcinoma tumor
    # #TCGA-COAD
    # #look up workflow TSeq - FPKM-UQ
    # query.met <- GDCquery(project = "TCGA-COAD",
    #                       data.category = "DNA Methylation",
    #                       legacy = FALSE,
    #                       platform = c("Illumina Human Methylation 450"))
    # query.exp <- GDCquery(project = "TCGA-COAD",
    #                       data.category = "Transcriptome Profiling",
    #                       data.type = "Gene Expression Quantification", 
    #                       workflow.type = "HTSeq - FPKM-UQ")
    # 
    # # Get all patients that have DNA methylation and gene expression.
    # #the cases column has the barcodes
    # common.patients <- intersect(substr(getResults(query.met, cols = "cases"), 1, 12),
    #                              substr(getResults(query.exp, cols = "cases"), 1, 12))
    # 
    # # Only select the first 5 patients
    # # for dna methylation
    # query.met <- GDCquery(project = "TCGA-COAD",
    #                       data.category = "DNA Methylation",
    #                       legacy = FALSE,
    #                       platform = c("Illumina Human Methylation 450"),
    #                       barcode = common.patients[1:5])
    # # and gene expression
    # query.exp <- GDCquery(project = "TCGA-COAD",
    #                       data.category = "Transcriptome Profiling",
    #                       data.type = "Gene Expression Quantification", 
    #                       workflow.type = "HTSeq - FPKM-UQ",
    #                       barcode = common.patients[1:5])
    # 
    # #now create data table for dna methylation data
    # #only showing data type and barcode
    # datatable(getResults(query.met, cols = c("data_type","cases")),
    #           filter = 'top',
    #           options = list(scrollX = TRUE, keys = TRUE, pageLength = 6), 
    #           rownames = FALSE)
    # #and gene expression data
    # datatable(getResults(query.exp, cols = c("data_type","cases")), 
    #           filter = 'top',
    #           options = list(scrollX = TRUE, keys = TRUE, pageLength = 6), 
    #           rownames = FALSE)
    # 
    # 
    # #now how to make sure you have the right barcodes associated with the right file names
    # #this is for breast cancer
    # #TCGA-BRCA
    # query <- GDCquery(project = c("TCGA-BRCA"),
    #                   data.category = "Raw Sequencing Data",  
    #                   sample.type = "Primary solid Tumor")
    # # Only first 10 to make render faster
    # datatable(getResults(query, rows = 1:10,cols = c("file_name","cases")), 
    #           filter = 'top',
    #           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    #           rownames = FALSE)
    #
    # #now to find dna methylation in legacy database with two different types of methylation
    # #illumina 450 and ilumina 27
    # query <- GDCquery(project = c("TCGA-GBM","TCGA-LGG"),
    #                   legacy = TRUE,
    #                   data.category = "DNA methylation",
    #                   platform = c("Illumina Human Methylation 450", "Illumina Human Methylation 27"))
    # #once again only getting first 10 rows
    # datatable(getResults(query, rows = 1:10), 
    #           filter = 'top',
    #           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    #           rownames = FALSE)
    #
    # 
    # # actually downloading data
    # # Gene expression aligned against hg38
    # # dowloading data from two samples
    # query <- GDCquery(project = "TCGA-GBM",
    #                   data.category = "Transcriptome Profiling",
    #                   data.type = "Gene Expression Quantification", 
    #                   workflow.type = "HTSeq - FPKM-UQ",
    #                   barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"))
    # GDCdownload(query)
    # genedata <- GDCprepare(query)
    # #data on two patient samples
    # datatable(as.data.frame(colData(genedata)),
    #           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
    #           rownames = FALSE)
    # #gene id and expression data for both samples
    # datatable(assay(genedata)[1:100,],
    #          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
    #          rownames = TRUE)
    # 
    # rowRanges(genedata)
}



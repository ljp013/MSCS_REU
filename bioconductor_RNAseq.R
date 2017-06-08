#Laura Poulton
#MSCS REU - Discovering Significant Pathways of Gene Regulation
#6/08/17

#going through the workflow at https://www.bioconductor.org/help/workflows/RNAseq123/
#to install Bioconductor 3.5:
## source("https://bioconductor.org/biocLite.R")
## biocLite()
#May need to type "defaults write org.R-project.R force.LANG en_US.UTF-8" into terminal
#to get packages to download properly

#to install this workflow:
## source("http://bioconductor.org/workflows.R")
## workflowInstall("RNAseq123")
## biocLite(c("package names")) to install any packages not initially installed

#this runs everything at once
main <- function() {
    #setup libraries
    library.setup()
    
    #data packaging
    #download tar file
    ## load.data()
    
    genecounts <- initialize.data()
    
    #first remove first 11 characters from each file name (GSM ids) to simplify things
    colnames(genecounts) <- simplify.names(genecounts)
    
    #assign cell types and lane numbers
    genecounts$samples <- group.data(genecounts)
    
    #now assign gene ids, symbols, and chromosomes
    genecounts$genes <- id.genes(genecounts)
    genecounts
    
    #data preprocessing
}


#Setup libraries
library.setup <- function() {
    library(limma) #modeling functions
    library(Glimma) #interactive modeling functions
    library(edgeR) 
    library(Mus.musculus) #mouse data
    #should load all required packages
}


#Downloading initial data set
load.data <- function() {
    url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
    # :: means in the utils package use this function (in this case download.file)
    utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
    utils::untar("GSE63310_RAW.tar", exdir = ".")
    files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
               "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
               "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
    for(i in paste(files, ".gz", sep=""))
        R.utils::gunzip(i, overwrite=TRUE)
}

#all the text files have "raw gene level counts"
#now to get out the data that we want

initialize.data <- function() {
    files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", 
               "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
               "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
               "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
               "GSM1545545_JMS9-P8c.txt")
    read.delim(files[1], nrow=5)
    #this will show the first five rows of the first file
    
    #we want all the files in one matrix, so we use the readDGE function from edgeR
    #we also only want first and third columns, gene id and count
    genecounts <- readDGE(files, columns=c(1,3))
    #try out class(genecounts) and dim(genecounts) (has 27179 rows and 9 columns)
    #9 files combined
    dim(genecounts)
    genecounts
}

simplify.names <- function(x) {
    samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
    samplenames
}

#want further analysis to include amy "experimental variables, both biological and technical, that could have an effect on expression levels"
#in this sample case, they are comparing gene count and cell type
#possible cell types: basal, LP, and ML
#the DGE-list object obtained using the readDGE function above as a "samples data frame
# that stores both cell type and batch (sequencing lane) information"
# 3 levels for each (lane 4, 6, or 8) (group Basal, LP, or ML)
group.data <- function (x) {
    
    #now assign cell types
    group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP"))
    #dollar signs go to that data frame
    x$samples$group <- group
    
    #now assign lanes
    #rep means repeat first argument number of times specifed in second argument
    lane <- as.factor(rep(c("L004", "L006", "L008"), c(3,4,2)))
    x$samples$lane <- lane
    x$samples
}

#it is also possible to extract gene ids, symbols, and chromosome names
# from the Mus.musculus package in the genes data frame of the DGEList-object
id.genes <- function(genecounts) {
    geneid <- rownames(genecounts)
    genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype= "ENTREZID")

    head(genes)
    #check for and remove duplicated gene ids
    genes <- genes[!duplicated(genes$ENTREZID),]
}


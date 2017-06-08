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
    #cpm is counts per million
    #lcpm is log2 counts per million
    cpm <- cpm(genecounts)
    lcpm <- cpm(genecounts, log=TRUE)
    
    #removing genes are not expressed in all samples for all cell types
    #"typically a CPM value of 1 is used [as a cutoff] in our analyses as it separates
    # expressed genes from unexpressed genes well for most datasets"
    genecounts <- del.unexpressed(genecounts, cpm)
    dim(genecounts)
    
    #now graph raw (current lcpm) vs filtered (lcpm <- cpm(new genecounts, log = TRUE)) data
    raw.vs.filtered(genecounts, lcpm, samplenames=colnames(genecounts))
    
    #normalising distributions
    #calcNormFactors function is from edgeR
    genecounts <- calcNormFactors(genecounts, method="TMM")
    genecounts$samples$norm.factors
    
    #understanding the power of normalising
    norm.plots(genecounts)
    
    #unsupervised clustering of samples
    #using multi-dimensional scaling (MDS) plots
    MDS.plots(genecounts)
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

#x is genecount data, cpm is counts per million of data
del.unexpressed <- function(x, cpm) {
    keep.exprs <- rowSums(cpm>1)>=3
    x <- x[keep.exprs,, keep.lib.sizes=FALSE]
    x
}

#create plot of gene expression density with raw and filtered data
#x is filtered gene count data, lcpm is the log2 counts per million of unfiltered data
raw.vs.filtered <- function(x, lcpm, samplenames) {
    #uses specific color library
    library(RColorBrewer)
    #number samples is the number of columns in gene count data
    nsamples <- ncol(x)
    #set colors from "paired" set
    col <- brewer.pal(nsamples, "Paired")
    #want two plots side by side (so 1 row by 2 col matrix)
    par(mfrow=c(1,2))
    #plot density of unfiltered lcpm data for first column (first samples)
    #use first color, line width of 2, not setting title or x label yet
    #las=2 sets tick marks perpendicular to the axis
    plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
         main="", xlab="")
    #Now set title and x label
    title(main="A. Raw data", xlab="Log-cpm")
    #add a line, lty=3 means dotted, v=0 means at 0 (cpm = 1 or lcpm = 0 was the cutoff)
    abline(v=0, lty=3)
    #plot density of unfiltered lcpm data for rest of samples
    #already did 1 so start 2:nsamples
    for (i in 2:nsamples){
        den <- density(lcpm[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
    }
    #legend has samplenames, same colors as line colors
    #bty is type of box, 'n' none or 'o' on
    legend("topright", samplenames, text.col=col, bty="n")
    #now find lcpm of filtered data
    lcpm <- cpm(x, log=TRUE)
    #and plot that first sample
    plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
         main="", xlab="")
    title(main="B. Filtered data", xlab="Log-cpm")
    #add cutoff line
    abline(v=0, lty=3)
    #and the rest of the samples
    for (i in 2:nsamples){
        den <- density(lcpm[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
    }
    legend("topright", samplenames, text.col=col, bty="n")
}

#x is gene count data
#shows what calcNormFactors function can do
norm.plots <- function(x) {
    x2 <- x
    x2$samples$norm.factors <- 1
    #first set two sample counts to be anomalous
    #first sample to 5% of original
    #second sample to 5x original
    x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
    x2$counts[,2] <- x2$counts[,2]*5
    
    #set colors
    col <- brewer.pal(ncol(x), "Paired")
    
    #now plot that vs. normalised version
    #want 2 graphs side by side
    par(mfrow=c(1,2))
    #find lcpm
    lcpm <- cpm(x2, log=TRUE)
    #they are using a boxplot this time
    boxplot(lcpm, las=2, col=col, main="")
    title(main="A. Example: Unnormalised data",ylab="Log-cpm")
    
    #now calculate normalised and plot that
    x2 <- calcNormFactors(x2)  
    lcpm <- cpm(x2, log=TRUE)
    boxplot(lcpm, las=2, col=col, main="")
    title(main="B. Example: Normalised data",ylab="Log-cpm")
}

#creates MDS plots grouped by cell types and sequencing lanes
#uses plotMDS function in limma
#x is gene count data
MDS.plots <- function(x) {
    #find lcpms of data
    lcpm <- cpm(x, log=TRUE)
    #want two graphs side by side
    par(mfrow=c(1,2))
    #find cell types
    col.group <- x$samples$group
    #use colors from set 1
    levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
    col.group <- as.character(col.group)
    #find lanes
    col.lane <- x$samples$lane
    #use colors from set 2
    levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
    col.lane <- as.character(col.lane)
    #if you don't specify dimensions assumes 1 and 2
    plotMDS(lcpm, labels=x$samples$group, col=col.group)
    title(main="A. Sample groups")
    #for lanes they specified dimensions 3 and 4
    plotMDS(lcpm, labels=x$samples$lane, col=col.lane, dim=c(3,4))
    title(main="B. Sequencing lanes")
    
    #Glimma package has interactive option
    #launch=TRUE opens webpage
    glMDSPlot(lcpm, labels=paste(x$samples$group, x$samples$lane, sep="_"), 
              groups=x$samples[,c(2,5)], launch=TRUE)
}

library(data.table)
library(edgeR)
workingdirectory <- Sys.getenv("WORKDIR")
sampledatafile <- Sys.getenv("SAMPLEDATAFILE")
normalisationmethod <- Sys.getenv("NORMALISATIONMETHOD")
normalisationscale <- as.numeric(Sys.getenv("NORMALISATIONSCALE"))
normalisationdecimals <- as.numeric(Sys.getenv("NORMALISATIONDECIMALS"))

#Load files
hit.table <- data.frame(fread(paste(workingdirectory,"/GeneTables/GeneHitTable.csv",sep=""),sep=",",header=TRUE),row.names=1)

#Define groups
sample.data <- read.table(sampledatafile,header=FALSE)
if (ncol(sample.data) > 3){
samples <- colnames(hit.table)
groups <- sample.data[sample.data[,1] %in% samples,4]
}else{
print("Group column does not exist.")
}

#Create DGEL list object
dgList <- DGEList(counts=hit.table, lib.size = colSums(hit.table), genes=rownames(hit.table), group=groups, remove.zeros=TRUE)

#Calculate and apply normalisation factors

#Method TMM
dgList.tmm <- calcNormFactors(dgList, method="TMM")
tmm.nf <- dgList.tmm$samples$norm.factors
hit.table.tmm <- round(sweep(hit.table, 2, tmm.nf, FUN="*"))

#Method RLE
dgList.RLE <- calcNormFactors(dgList, method="RLE")
RLE.nf <- dgList.RLE$samples$norm.factors
hit.table.RLE <- round(sweep(hit.table, 2, RLE.nf, FUN="*"))

#Method UQ - PROBABLY NOT WORKING!!
dgList.UQ <- calcNormFactors(dgList, method="upperquartile")
UQ.nf <- dgList.UQ$samples$norm.factors
hit.table.UQ <- round(sweep(hit.table, 2, UQ.nf, FUN="*"))

#Method TSS
hit.table.tss <- round(sweep(hit.table, 2, totals, FUN="/") * normalisationscale,normalisationdecimals)

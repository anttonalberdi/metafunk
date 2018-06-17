library(data.table)
library(edgeR)
workingdirectory <- Sys.getenv("WORKDIR")
sampledatafile <- Sys.getenv("SAMPLEDATAFILE")
normalisationscale <- as.numeric(Sys.getenv("NORMALISATIONSCALE"))
normalisationdecimals <- as.numeric(Sys.getenv("NORMALISATIONDECIMALS"))

#Load files
hit.table <- data.frame(fread(paste(workingdirectory,"/GeneTables/GeneHitTable.csv",sep=""),sep=",",header=TRUE),row.names=1)

#Define groups
sample.data <- read.table(sampledatafile,header=FALSE)
groups <- sample.data[,4]

#Select reference samples


#Create DGEL list object
dgList <- DGEList(counts=hit.table, lib.size = colSums(hit.table), genes=rownames(hit.table), group=groups, remove.zeros=TRUE)

#Calculate and apply normalisation factors
dgList.tmm <- calcNormFactors(dgList, method="TMM")
tmm.nf <- dgList.tmm$samples$norm.factors
hit.table.tmm <- round(sweep(hit.table, 2, tmm.nf, FUN="*"))

dgList.RLE <- calcNormFactors(dgList, method="RLE")
RLE.nf <- dgList.RLE$samples$norm.factors
hit.table.RLE <- round(sweep(hit.table, 2, RLE.nf, FUN="*"))

dgList.UQ <- calcNormFactors(dgList, method="upperquartile")
UQ.nf <- dgList.UQ$samples$norm.factors
hit.table.UQ <- round(sweep(hit.table, 2, UQ.nf, FUN="*"))


hit.table.tss <- round(sweep(hit.table, 2, totals, FUN="/") * normalisationscale,normalisationdecimals)

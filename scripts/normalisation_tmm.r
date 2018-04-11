library(data.table)
library(edgeR)
workingdirectory <- Sys.getenv("WORKDIR")
normalisationscale <- as.numeric(Sys.getenv("NORMALISATIONSCALE"))
normalisationdecimals <- as.numeric(Sys.getenv("NORMALISATIONDECIMALS"))

#Load files
hit.table <- data.frame(fread(paste(workingdirectory,"/GeneTables/GeneHitTable.csv",sep=""),sep=",",header=TRUE),row.names=1)

#Define groups
groups <- c(1,1,1,1,1,1,2,2,2,3,3,3,3,3,3)

#Create DGEL list object
dgList <- DGEList(counts=hit.table, genes=rownames(hit.table), group=groups)

#Calculate normalisation factors
dgList <- calcNormFactors(dgList, method="TMM")

pdf("nmds.pdf")
plotMDS(dgList)
dev.off()

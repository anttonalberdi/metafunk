library(data.table)
workingdirectory <- Sys.getenv("WORKDIR")
normalisationscale <- as.numeric(Sys.getenv("NORMALISATIONSCALE"))
normalisationdecimals <- as.numeric(Sys.getenv("NORMALISATIONDECIMALS"))

#Load files
hit.table <- data.frame(fread(paste(workingdirectory,"/GeneTables/GeneHitTable.csv",sep=""),sep=",",header=TRUE),row.names=1)
gene.lengths <- data.frame(fread(paste(workingdirectory,"/GenePrediction/assembly.genes.lengths",sep=""),sep="\t",header=FALSE),row.names=1)
gene.lengths <- gene.lengths[rownames(hit.table),]
totals <- colSums(hit.table)
hit.table.tss <- round(sweep(hit.table, 2, totals, FUN="/") * normalisationscale,normalisationdecimals)
cov.table.tss <- sweep(hit.table.tss, 1, gene.lengths, FUN="/")

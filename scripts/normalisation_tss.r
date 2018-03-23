library(data.table)
workingdirectory <- Sys.getenv("WORKDIR")
normalisationscale <- Sys.getenv("NORMALISATIONSCALE")
normalisationdecimals <- Sys.getenv("NORMALISATIONDECIMALS")

#Load files
hit.table <- data.frame(fread(paste(workingdirectory,"/GeneHitTable.csv",sep=""),sep=",",header=TRUE),row.names=1)
gene.lengths <- data.frame(fread(paste(workingdirectory,"GenePrediction/assembly.genes.lengths",sep=""),sep="\t",header=FALSE),row.names=1)
totals <- colSums(hit.table)
hit.table.tss <- round(sweep(hit.table, 2, weigthing.values, FUN="/"),normalisationdecimals) * normalisationscale
cov.table.tss <- round(sweep(hit.table.tss, 1, gene.lengths, FUN="/"),normalisationdecimals)

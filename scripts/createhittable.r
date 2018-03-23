library(data.frame)
workingdirectory <- Sys.getenv("WORKDIR")
cov.table <- data.frame(fread(paste(workingdirectory,"CoverageTable.csv",sep=""),sep=",",header=TRUE),row.names=1)
gene.lengths <- data.frame(fread(paste(workingdirectory,"GenePrediction/assembly.genes.lengths",sep=""),sep="\t",header=FALSE),row.names=1)
gene.lengths <- gene.lengths[rownames(cov.table),]
hit.table <- round(sweep(cov.table, 1, gene.lengths, FUN="*"),0)
write.table(hit.table,paste(workingdirectory,"HitTable.csv",sep=""),sep=",",col.names=TRUE,row.names=FALSE,quote=FALSE)

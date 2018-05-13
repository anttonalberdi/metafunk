library(data.table)
workingdirectory <- Sys.getenv("WORKDIR")
cov.table <- data.frame(fread(paste(workingdirectory,"/GeneTables/GeneCoverageTable.csv",sep=""),sep=",",header=TRUE),row.names=1)
gene.lengths <- data.frame(fread(paste(workingdirectory,"/GenePrediction/assembly.genes.lengths",sep=""),sep="\t",header=FALSE),row.names=1)
gene.lengths <- gene.lengths[rownames(cov.table),]
hit.table <- round(sweep(cov.table, 1, gene.lengths, FUN="*"),0)
contig <- rownames(hit.table)
hit.table2 <- cbind(contig,hit.table)
write.table(hit.table2,paste(workingdirectory,"/GeneTables/GeneHitTable.csv",sep=""),sep=",",col.names=TRUE,row.names=FALSE,quote=FALSE)

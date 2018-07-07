library(data.table)
library(edgeR)
workingdirectory <- Sys.getenv("WORKDIR")
sampledatafile <- Sys.getenv("SAMPLEDATAFILE")
normalisationmethod <- Sys.getenv("NORMALISATIONMETHOD")

#Load files
hit.table <- data.frame(fread(paste(workingdirectory,"/GeneTables/GeneHitTable.csv",sep=""),sep=",",header=TRUE),row.names=1)
gene.lengths <- data.frame(fread(paste(workingdirectory,"/GenePrediction/assembly.genes.lengths",sep=""),sep="\t",header=FALSE))
colnames(gene.lengths) <- c("gene","length")
gene.lengths <- gene.lengths[order(gene.lengths[,1]),] 

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
if (grepl("tmm", normalisationmethod) == TRUE){
dgList.tmm <- calcNormFactors(dgList, method="TMM")
tmm.nf <- dgList.tmm$samples$norm.factors
hit.table.tmm <- round(sweep(hit.table, 2, tmm.nf, FUN="*"))
hit.table.tmm <- hit.table.tmm[order(row.names(hit.table.tmm)),] 
gene.lengths.subset <- gene.lengths[gene.lengths[,1] %in% rownames(hit.table.tmm),]
coverage.table.tmm <- sweep(hit.table.tmm, 1, gene.lengths.subset[,2], FUN="/")
write.table(coverage.table.tmm,paste(workingdirectory,"/GeneTables/GeneCoverageTable.tmm.csv",sep=""),row.names=TRUE, col.names=TRUE,sep=",",quote=FALSE)
write.table(hit.table.tmm,paste(workingdirectory,"/GeneTables/GeneHitTable.tmm.csv",sep=""),row.names=TRUE, col.names=TRUE,sep=",",quote=FALSE)
}

#Method RLE
if (grepl("rle", normalisationmethod) == TRUE){
dgList.RLE <- calcNormFactors(dgList, method="RLE")
rle.nf <- dgList.RLE$samples$norm.factors
hit.table.rle <- round(sweep(hit.table, 2, rle.nf, FUN="*"))
hit.table.rle <- hit.table.rle[order(row.names(hit.table.rle)),] 
gene.lengths.subset <- gene.lengths[gene.lengths[,1] %in% rownames(hit.table.rle),]
coverage.table.rle <- sweep(hit.table.rle, 1, gene.lengths.subset[,2], FUN="/")
write.table(coverage.table.rle,paste(paste(workingdirectory,"/GeneTables/GeneCoverageTable.rle.csv",sep=""),row.names=TRUE, col.names=TRUE,sep=",",quote=FALSE)
write.table(hit.table.rle,paste(workingdirectory,"/GeneTables/GeneHitTable.rle.csv",sep=""),row.names=TRUE, col.names=TRUE,sep=",",quote=FALSE)
}

#Method UQ - PROBABLY NOT WORKING!!
if (grepl("uq", normalisationmethod) == TRUE){
dgList.UQ <- calcNormFactors(dgList, method="upperquartile")
UQ.nf <- dgList.UQ$samples$norm.factors
hit.table.UQ <- round(sweep(hit.table, 2, UQ.nf, FUN="*"))
}
  
#Method TSS
if (grepl("tss", normalisationmethod) == TRUE){
#Identify the sample with the seq depth closest to the mean
seqdepths <- colSums(hit.table)
y <- which.min(abs(seqdepths - mean(seqdepths))) 
reference <- seqdepths[y]
# 
tss.nf <- round(reference/seqdepths, 5)
hit.table.tss <- round(sweep(hit.table, 2, tss.nf, FUN="*"))
hit.table.tss <- hit.table.tss[order(row.names(hit.table.tss)),] 
gene.lengths.subset <- gene.lengths[gene.lengths[,1] %in% rownames(hit.table.tss),]
coverage.table.tss <- sweep(hit.table.tss, 1, gene.lengths.subset[,2], FUN="/")
write.table(coverage.table.tss,paste(paste(workingdirectory,"/GeneTables/GeneCoverageTable.tss.csv",sep=""),row.names=TRUE, col.names=TRUE,sep=",",quote=FALSE)
write.table(hit.table.tss,paste(workingdirectory,"/GeneTables/GeneHitTable.tss.csv",sep=""),row.names=TRUE, col.names=TRUE,sep=",",quote=FALSE)
}

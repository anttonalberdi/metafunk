
library(data.table)
workingdirectory <- Sys.getenv("WORKDIR")
metafunkdirectory <- Sys.getenv("METAFUNKDIR")
sampledatafile <- Sys.getenv("SAMPLEDATAFILE")
normalisationmethod <- Sys.getenv("NORMALISATIONMETHOD")
keggthreshold <- Sys.getenv("KEGGTHRESHOLD")
timestamp <- Sys.getenv("TIMESTAMP")

#Choose normalisation methods
methods <- strsplit(normalisationmethod, ",")[[1]]

#Load and prepare annotation file
annot.table <- fread(paste(workingdirectory,"/GeneAnnotationKEGG/assembly.genes.KEGG.annotated.",keggthreshold,".txt",sep=""),sep="\t",header=FALSE,,colClasses=list(character=c("V4")))
colnames(annot.table) <- c("contig","KEGG","KO","Path","e-value","Identity")
KO.table <- annot.table[,c("contig","KO")]
KO.table <- KO.table[!duplicated(KO.table), ]
KOtoPath <- annot.table[,c("KO","Path")]
KOtoPath <- KOtoPath[!duplicated(KOtoPath), ]
KOtoPath <- KOtoPath[complete.cases(KOtoPath), ]

#Loop across normalisation methods
for (method in methods){
text=paste("      Aggregating KO and Path values for normalisation method: ",method,sep="")  
write(text,file=paste(workingdirectory,"/run_",timestamp,".log",sep=""),append=TRUE)
cov.table <- fread(paste(workingdirectory,"/GeneTables/GeneCoverageTable.",method,".csv",sep=""),sep=",")
colnames(cov.table)[1] <- "contig"
  
#KO aggregation
cov.KO.table <- merge(cov.table,KO.table,by="contig")
cov.KO.table.aggregated <- aggregate(subset(cov.KO.table, select = -c(contig,KO)),by=list(cov.KO.table$KO),FUN=sum)
colnames(cov.KO.table.aggregated)[1]<-"KO"
write.table(cov.KO.table.aggregated,paste(workingdirectory,"/GeneTables/GeneCoverageTable.",method,".KO.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

#Path aggregation
cov.KO.Path.table <- merge(cov.KO.table,KOtoPath,by="KO",allow.cartesian=TRUE)
cov.Path.table.aggregated <- aggregate(subset(cov.KO.Path.table, select = -c(KO,contig,Path)),by=list(cov.KO.Path.table$Path),FUN=sum)
colnames(cov.Path.table.aggregated)[1]<-"Path"
write.table(cov.Path.table.aggregated,paste(workingdirectory,"/GeneTables/GeneCoverageTable.",method,".Path.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

#Limit Path aggregation to metabolism
KEGG.metabolism <- data.frame(fread(paste(metafunkdirectory,"/files/KEGG.metabolism.csv",sep=""),header=TRUE,sep=",",colClasses="character"))
cov.Path.table.aggregated.metabolism <- cov.Path.table.aggregated[cov.Path.table.aggregated$Path %in% KEGG.metabolism$Pathway,]
write.table(cov.Path.table.aggregated.metabolism,paste(workingdirectory,"/GeneTables/GeneCoverageTable.",method,".Path.metabolism.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)


#Metabolic domain aggregation
cov.Path.table.aggregated.metabolism.domain <- merge(cov.Path.table.aggregated.metabolism,KEGG.metabolism[,c(1,3)],by.x="Path",by.y="Pathway")
cov.domain.table.aggregated <- aggregate(subset(cov.Path.table.aggregated.metabolism.domain, select = -c(Path,Class)),by=list(cov.Path.table.aggregated.metabolism.domain$Class),FUN=sum)
colnames(cov.domain.table.aggregated)[1] <- "Domain"
write.table(cov.domain.table.aggregated,paste(workingdirectory,"/GeneTables/GeneCoverageTable.",method,".Domain.metabolism.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)
}

##### Run statistical analyses ######

cov.Path.metabolism <- data.frame(fread(paste(workingdirectory,"/GeneTables/GeneCoverageTable.",method,".Path.metabolism.csv",sep=""),sep=",",header=TRUE,colClasses=list(character=c("Path"))),row.names=1)
cov.Domain.metabolism <- data.frame(fread(paste(workingdirectory,"/GeneTables/GeneCoverageTable.",method,".Domain.metabolism.csv",sep=""),sep=",",header=TRUE,colClasses=list(character=c("Domain"))),row.names=1)
cov.KO <- data.frame(fread(paste(workingdirectory,"/GeneTables/GeneCoverageTable.",method,".KO.csv",sep=""),sep=",",header=TRUE,colClasses=list(character=c("KO"))),row.names=1)

#Define number of groups
sampledata <- read.table(sampledatafile,row.names=1)
sample.group <- cbind(rownames(sampledata),sampledata[,3])
colnames(sample.group) <- c("Sample","Group")
groups <- sampledata[,3]
groupnames <- unique(groups)
groupnumber <- length(groupnames)

if (groupnumber == 1){
"It is not possible to contrast groups, as only one group has been defined"
}
if (groupnumber == 2){

#Define groups
group1 <- rownames(sampledata[sampledata[,3] == groupnames[1],])
group2 <- rownames(sampledata[sampledata[,3] == groupnames[2],])

#Domain level Wilcoxon
    domains <- rownames(cov.Domain.metabolism)
    domain.table <- c()
    for (domain in domains){
    x <- as.numeric(cov.Domain.metabolism[domain,group1])
    y <- as.numeric(cov.Domain.metabolism[domain,group2])
    pvalue <- wilcox.test(x,y)$p.value
    row <- cbind(domain,pvalue)
    domain.table <- rbind(domain.table,row)
    }
    write.table(domain.table,paste(workingdirectory,"/FunctionalStats/KEGG.domain.metabolism.",method,".wilcoxontest.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

  
#Pathway level Wilcoxon
    pathways <- rownames(cov.Path.metabolism)
    pathway.table <- c()
    for (path in pathways){
    x <- as.numeric(cov.Path.metabolism[path,group1])
    y <- as.numeric(cov.Path.metabolism[path,group2])
    pvalue <- wilcox.test(x,y)$p.value
    row <- cbind(path,pvalue)
    pathway.table <- rbind(pathway.table,row)
    }
    write.table(pathway.table,paste(workingdirectory,"/FunctionalStats/KEGG.pathway.metabolism.",method,".wilcoxontest.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

#KO level Wilcoxon
    KOs <- rownames(cov.KO)
    KO.table <- c()
    for (KO in KOs){
    x <- as.numeric(cov.KO[KO,group1])
    y <- as.numeric(cov.KO[KO,group2])
    pvalue <- wilcox.test(x,y)$p.value
    row <- cbind(KO,pvalue)
    KO.table <- rbind(KO.table,row)
    }
    write.table(KO.table,paste(workingdirectory,"/FunctionalStats/KEGG.KO.",method,".wilcoxontest.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

#Domain level Jitterplot
domain.melt <- cbind(rep(domains,length(sampledata[,3])),melt(cov.Domain.metabolism))
domain.melt2 <- merge(domain.melt,sample.group,by.x="variable",by.y="Sample")
colnames(domain.melt) <- c("Domain","Sample","Value","Group")
  
}
if (groupnumber > 2){

}

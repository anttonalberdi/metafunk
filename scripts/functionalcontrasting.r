
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
}

##### Run statistical analyses ######

#Define number of groups
sampledata <- read.table(sampledatafile,row.names=1)
groupnumber <- length(unique(sampledata[,3]))

if (groupnumber == 1){
"It is not possible to contrast groups, as only one group has been defined"
}
if (groupnumber == 2){

}

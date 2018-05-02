library(data.table)
workingdirectory <- Sys.getenv("WORKDIR")
metafunkdirectory <- Sys.getenv("METAFUNKDIR")

#Load coverage and annotation tables
cov.table <- data.frame(fread(paste(workingdirectory,"/GeneTables/GeneCoverageTable.csv",sep=""),sep=",",header=TRUE),row.names=1)
gene.annotations <- data.frame(fread(paste(workingdirectory,"/GeneAnnotationKEGG/assembly.genes.KEGG.annotated.txt",sep=""),sep="\t",header=TRUE))

#Merge coverage and annotation tables
cov.table.annotated <- merge(cov.table,gene.annotations,by.x="row.names",by.y="Gene")

#Load USiCG files
USiCGs <- read.table(paste(metafunkdirectory,"/files/USiCGs.txt",sep=""),sep="\t",header=TRUE)
USiCGs <- USiCGs[,1]

#Define non-coverage columns
drops <- c("Row.names","KEGG","KO","Pathway","evalue","pident")

#Loop through USiCGs
estimationtable <- c()

for (i in USiCGs){
subset <- cov.table.annotated[which(cov.table.annotated$KO == i),]
subset <- subset[!duplicated(subset$Row.names),]
aggregated <- aggregate(subset[,!(colnames(subset) %in% drops)],by=list(subset$KO),FUN=sum)
estimationtable <- rbind(estimationtable,aggregated)
}

#Create summary table
sample <- colnames(estimationtable[,-1])
mean <- apply(estimationtable[,-1], 2, mean)
sd <- apply(estimationtable[,-1], 2, sd)
finaltable <- cbind(sample,mean,sd)

#Write table
write.table(finaltable,paste(workingdirectory,"/GenomeEstimation/estimation_kegg.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

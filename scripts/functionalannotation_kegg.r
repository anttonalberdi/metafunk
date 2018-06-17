library(data.table)
workingdirectory <- Sys.getenv("WORKDIR")
genes_ko_file <- Sys.getenv("GENES_KO")
genes_path_file <- Sys.getenv("GENES_PATH")

#####
# 1) Merge Entry, KO and Path codes
#####

#Load and modify KEGG entry list
KEGG.entry <- fread(paste(workingdirectory,"/GeneAnnotationKEGG/assembly.genes.KEGG.entrylist.txt",sep=""),sep="\t",head=FALSE)
names(KEGG.entry) <- "Gene"

#Load and modify Entry to KO mapping file
KEGG.entry_KO <- fread(genes_ko_file,sep="\t",head=FALSE)
colnames(KEGG.entry_KO) <- c("Gene","KO")
KEGG.entry_KO$KO <- gsub("^ko:","",KEGG.entry_KO$KO)

#Load and modify Entry to Pathway mapping file
KEGG.entry_path <- fread(genes_path_file,sep="\t",head=FALSE)
colnames(KEGG.entry_path) <- c("Gene","Pathway")
KEGG.entry_path$Pathway <- gsub("[^0-9]","",KEGG.entry_path$Pathway)

#Merge Entry, KO and Path
KEGG.entry.KO <- merge(KEGG.entry,KEGG.entry_KO,by="Gene",all.x=TRUE)
KEGG.entry.KO.path <- merge(KEGG.entry.KO,KEGG.entry_path,by="Gene",all.x=TRUE)
colnames(KEGG.entry.KO.path) <- c("KEGG","KO","Pathway")

#####
# 2) Append functional information to Diamond results
#####

#Load and modify Diamond blastp result file
KEGG.table <- fread(paste(workingdirectory,"/GeneAnnotationKEGG/assembly.genes.KEGG.txt",sep=""),sep="\t",header=FALSE)
colnames(KEGG.table) <- c("Gene","KEGG","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

#Merge Diamond result file and annotation
KEGG.merged <- merge(KEGG.table[,c("Gene","KEGG","evalue","pident")],KEGG.entry.KO.path,by="KEGG",all.x=TRUE,allow.cartesian=TRUE)
KEGG.merged.sorted <- KEGG.merged[order(KEGG.merged$Gene),c("Gene","KEGG","KO","Pathway","evalue","pident")]

#Output annotated table
write.table(KEGG.merged.sorted,paste(workingdirectory,"/GeneAnnotationKEGG/assembly.genes.KEGG.annotated.txt",sep=""),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

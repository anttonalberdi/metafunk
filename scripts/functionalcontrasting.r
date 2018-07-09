
library(data.table)
workingdirectory <- Sys.getenv("WORKDIR")
sampledatafile <- Sys.getenv("SAMPLEDATAFILE")
normalisationmethod <- Sys.getenv("NORMALISATIONMETHOD")
keggthreshold <- Sys.getenv("KEGGTHRESHOLD")

#Choose normalisation methods
methods <- strsplit(normalisationmethod, ",")[[1]]

#Load annotation fikle
annot.table <- fread(paste(workingdirectory,"/GeneAnnotationKEGG/assembly.genes.KEGG.annotated.",keggthreshold,".txt",sep=""),sep="\t",header=FALSE,,colClasses=list(character=c("V4")))

#Loop across normalisation methods
for (method in methods){

cov.table <- fread(paste(workingdirectory,"/GeneTables/GeneCoverageTable.",method,".csv",sep=""),sep=",",header=TRUE)

}

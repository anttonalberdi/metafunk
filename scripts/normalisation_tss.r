library(data.frame)
workingdirectory <- Sys.getenv("WORKDIR")

#Load files
hit.table <- data.frame(fread(paste(workingdirectory,"/Hit.table.csv",sep=""),sep=",",header=TRUE),row.names=1)

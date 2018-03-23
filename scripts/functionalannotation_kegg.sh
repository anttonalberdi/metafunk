#Source dependencies
source "$metafunkdirectory/settings.sh"

mkdir -p ${workingdirectory}/${project}/GeneAnnotationKEGG

#Perform Diamond blastp
diamond blastp -d ${keggdatabase} -p ${threads} -q ${workingdirectory}/${project}/GenePrediction/assembly.genes.faa --out ${workingdirectory}/${project}/GeneAnnotationKEGG/assembly.KEGG.txt --outfmt 6 --max-target-seqs 1 --evalue ${keggevalue}

#Source settings file
source $settingsfile

#Create Coassembly directory
mkdir -p ${workdir}/GeneAnnotationKEGG

#Create diamond database
if [ ! -f ${keggdatabase}.dmnd ]; then
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 	Creating diamond database"
  diamond makedb -p ${threads} --in ${keggdatabase} -d ${keggdatabase}
  else
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 	Diamond database already exists"
fi

#Perform Diamond blastp
diamond blastp -d ${keggdatabase} -p ${threads} -q ${workdir}/GenePrediction/assembly.genes.faa --out ${workdir}/GeneAnnotationKEGG/assembly.genes.KEGG.txt --outfmt 6 --max-target-seqs 1 --evalue 0.01

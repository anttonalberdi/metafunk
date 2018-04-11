#Source settings file
source $settingsfile

#Create Coassembly directory
mkdir -p ${workdir}/GeneAnnotationKEGG

#Create diamond database
keggdatabaseext=$(echo ${keggdatabase} | awk -F . '{print $NF}')

if [[ $keggdatabaseext != "dmnd" ]]; then
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 	Creating diamond database from $keggdatabaseext file"
  diamond makedb -p ${threads} --in ${keggdatabase} -d ${keggdatabase}
  else
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 	Diamond database already exists"
fi

#Perform Diamond blastp
#diamond blastp -d ${keggdatabase} -p ${threads} -q ${workdir}/GenePrediction/assembly.genes.faa --out ${workdir}/GeneAnnotationKEGG/assembly.genes.KEGG.txt --outfmt 6 --max-target-seqs 1 --evalue 0.01

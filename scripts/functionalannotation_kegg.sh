#Source settings file
source $settingsfile

#Create Coassembly directory
mkdir -p ${workdir}/GeneAnnotationKEGG

#Check if
#Create diamond database
keggdatabaseext=$(echo ${keggdatabase} | awk -F . '{print $NF}')

if [[ $keggdatabaseext != "dmnd" ]]; then
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 	Creating diamond database from .$keggdatabaseext file" >>  ${workdir}/run_${timestamp}.log
  diamond makedb -p ${threads} --in ${keggdatabase} -d ${keggdatabase}
  keggdatabase=${keggdatabase}.dmnd
  else
  if [[ -s ${keggdatabase} ]]; then
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now | 	Diamond database already exists" >>  ${workdir}/run_${timestamp}.log
    else
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now | 	The indicated Diamond database does not exist" >>  ${workdir}/run_${timestamp}.log
  fi
fi

#Perform Diamond blastp
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | 	Performing Diamond blastp" >>  ${workdir}/run_${timestamp}.log
diamond blastp -d ${keggdatabase} -p ${threads} -q ${workdir}/GenePrediction/assembly.genes.faa --out ${workdir}/GeneAnnotationKEGG/assembly.genes.KEGG.txt --outfmt 6 --max-target-seqs 1 --evalue 0.01
if [[ -s ${workdir}/GeneAnnotationKEGG/assembly.genes.KEGG.txt ]]; then
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 	There was an error during the Diamond blastp" >>  ${workdir}/run_${timestamp}.log
else
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 	Diamond blastp was succesfully finished" >>  ${workdir}/run_${timestamp}.log
fi

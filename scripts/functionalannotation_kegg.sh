#Source settings file
source $settingsfile

#Create Coassembly directory
mkdir -p ${workdir}/GeneAnnotationKEGG

#Check if
#Create diamond database
keggdatabaseext=$(echo ${keggdatabase} | awk -F . '{print $NF}')

if [[ $keggdatabaseext != "dmnd" ]]; then
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 	Creating diamond database from $keggdatabaseext file" >>  ${workdir}/run_${timestamp}.log
  diamond makedb -p ${threads} --in ${keggdatabase} -d ${keggdatabase}
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
#diamond blastp -d ${keggdatabase} -p ${threads} -q ${workdir}/GenePrediction/assembly.genes.faa --out ${workdir}/GeneAnnotationKEGG/assembly.genes.KEGG.txt --outfmt 6 --max-target-seqs 1 --evalue 0.01

#Source settings file
source $settingsfile

#Create Coassembly directory
mkdir -p ${workdir}/GeneAnnotationEggNog

#Check if
#Create diamond database
eggnogdatabaseext=$(echo ${eggnogdatabase} | awk -F . '{print $NF}')

if [[ $eggnogdatabase != "dmnd" ]]; then
  if [[ -s ${eggnogdatabase}.dmnd ]]; then
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now | 	Diamond database already exists" >>  ${workdir}/run_${timestamp}.log
    eggnogdatabase=${eggnogdatabase}.dmnd
  else
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now | 	Creating diamond database from .$eggnogdatabaseext file" >>  ${workdir}/run_${timestamp}.log
    eggnogdatabasedir=$(echo $eggnogdatabase | sed 's%/[^/]*$%/%' )
    mkdir ${eggnogdatabasedir}/extract
    tar -xvf ${eggnogdatabase} -C ${eggnogdatabasedir}/extract
    cat ${eggnogdatabasedir}/extract/*/* > ${eggnogdatabasedir}.fasta
    diamond makedb -p ${threads} --in ${eggnogdatabasedir}.fasta -d ${eggnogdatabasedir}
    eggnogdatabase=${eggnogdatabase}.dmnd
    fi
  else
  if [[ -s ${eggnogdatabase} ]]; then
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
diamond blastp -d ${eggnogdatabase} -p ${threads} -q ${workdir}/GenePrediction/assembly.genes.faa --out ${workdir}/GeneAnnotationEggNog/assembly.genes.EggNog.txt --outfmt 6 --max-target-seqs 1 --evalue 0.01
if [[ -s ${workdir}/GeneAnnotationEggNog/assembly.genes.EggNog.txt ]]; then
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 	Diamond blastp was succesfully finished" >>  ${workdir}/run_${timestamp}.log
else
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 	There was an error during the Diamond blastp" >>  ${workdir}/run_${timestamp}.log
  exit
fi

#Get unique EggNog entry list
#while read line; do
#query=$(echo $line | cut -d$' ' -f2)
#   result=$(echo $line | cut -d$' ' -f2 | grep -f - ${eggnogmembers} | cut -f2)
#   echo -e "$query\t$result"
#done < ${workdir}/GeneAnnotationEggNog/assembly.genes.EggNog.txt > ${workdir}/GeneAnnotationEggNog/assembly.genes.EggNog.entrylist.txt



#cut -f2 ${workdir}/GeneAnnotationEggNog/assembly.genes.EggNog.txt | sort | uniq > ${workdir}/GeneAnnotationEggNog/assembly.genes.KEGG.entrylist.txt

#Assign KO and Pathway codes to KEGG entries
#echo "$now |    Annotating KEGG entries" >> ${workdir}/run_${timestamp}.log
#export WORKDIR="${workdir}"
#export ANNOTATIONS="${eggnogannotations}"
#export MEMBERS="${eggnogmembers}"
#
#if [[ -s ${workdir}/GeneAnnotationEggNog/assembly.genes.KEGG.annotated.txt ]]; then
#  now=$(date +"%Y-%m-%d %H:%M:%S")
#  echo "$now | 	Genes have been succesfully annotated" >>  ${workdir}/run_${timestamp}.log
#  else
#  now=$(date +"%Y-%m-%d %H:%M:%S")
#  echo "$now | 	There was a problem when annotating the genes" >>  ${workdir}/run_${timestamp}.log
#fi

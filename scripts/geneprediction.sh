#Source dependencies
source "$metafunkdirectory/settings.sh"

#Create Gene Prediction directory
mkdir -p ${workingdirectory}/${project}/GenePrediction

#Check if assembly file exists
if [ ! -f "${workingdirectory}/${project}/CoAssembly/Megahit/final.contigs.fa" ]
then
  now=$(date +"%Y-%d-%m %H:%M:%S")
  echo "$now | ERROR: Assembly file does not exist." >> ${workingdirectory}/${project}/run.log
fi

#Predict genes
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Predicting genes in the assembly" >> ${workingdirectory}/${project}/run.log
prodigal -p meta -q -i ${workingdirectory}/${project}/CoAssembly/Megahit/final.contigs.fa -f gff -o ${workingdirectory}/${project}/GenePrediction/assembly.genes.gff -a ${workingdirectory}/${project}/GenePrediction/assembly.genes.faa -d ${workingdirectory}/${project}/GenePrediction/assembly.genes.fna
genenumber=$(grep -c ">" ${workingdirectory}/${project}/GenePrediction/assembly.genes.faa)
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | $genenumber genes were predicted" >> ${workingdirectory}/${project}/run.log

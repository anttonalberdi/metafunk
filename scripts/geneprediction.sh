#Source settings file
source $settingsfile

#Create Gene Prediction directory
mkdir -p ${workdir}/GenePrediction

#Check if assembly file exists
if [ ! -f "${workdir}/CoAssembly/Megahit/final.contigs.fa" ]; then
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | ERROR: Assembly file does not exist." >> ${workdir}/run_${timestamp}.log
  exit
fi

#Predict genes
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Predicting genes in the assembly" >> ${workdir}/run.log
prodigal -p meta -q -i ${workdir}/CoAssembly/Megahit/final.contigs.fa -f gff -o ${workdir}/GenePrediction/assembly.genes.gff -a ${workdir}/GenePrediction/assembly.genes.faa -d ${workdir}/GenePrediction/assembly.genes.fna
genenumber=$(grep -c ">" ${workdir}/GenePrediction/assembly.genes.faa)
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | $genenumber genes were predicted" >> ${workdir}/run_${timestamp}.log

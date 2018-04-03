#Source settings file
source $settingsfile

#Create Gene Prediction directory
mkdir -p ${workdir}/GenePrediction

#Check if assembly file exists
if [ ! -f "${workdir}/CoAssembly/Megahit/final.contigs.fa" ]; then
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | ERROR: Assembly file does not exist." >> ${workdir}/run_${timestamp}.log
  exit
else
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | Splitting assembly file for parallelising" >> ${workdir}/run_${timestamp}.log
fi

#Split file to parallelize gene prediction

#Check if number of threads is even (if not, convert to even)
rows=$(cat ${workdir}/CoAssembly/Megahit/final.contigs.fa | wc -l)
subsetrows=$((rows/threads))
if [[ $((subsetrows%2)) -eq 0 ]];then
  subsetrows2=$subsetrows
  else
  subsetrows2=$((subsetrows+1))
fi
split -l $subsetrows2 ${workdir}/CoAssembly/Megahit/final.contigs.fa ${workdir}/CoAssembly/Megahit/final.contigs.fa.

#Predict genes
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Predicting genes in the assembly" >> ${workdir}/run_${timestamp}.log
#Declare function
function geneprediction() {

  sample=${1}
  settingsfile=${2}
  suffix=$(echo $sample | awk '{print $NF}' FS=.)

source $settingsfile

prodigal -p meta -q -i ${sample} -f gff -o ${workdir}/GenePrediction/assembly.genes.gff.${suffix} -a ${workdir}/GenePrediction/assembly.genes.faa.${suffix} -d ${workdir}/GenePrediction/assembly.genes.fna.${suffix}
}

#Create file list
filelist=$(ls ${workdir}/CoAssembly/Megahit/final.contigs.fa.*)
#Export function lowcompjob
export -f geneprediction
#Loop in parallel across samples specified in sample.data.txt
parallel -j ${threads} -k geneprediction {} ${settingsfile} < ${filelist}

#Merge all files
cat ${workdir}/GenePrediction/assembly.genes.gff.* > ${workdir}/GenePrediction/assembly.genes.gff
cat ${workdir}/GenePrediction/assembly.genes.faa.* > ${workdir}/GenePrediction/assembly.genes.faa
cat ${workdir}/GenePrediction/assembly.genes.fna.* > ${workdir}/GenePrediction/assembly.genes.fna

#Count number of predicted genes
genenumber=$(grep -c ">" ${workdir}/GenePrediction/assembly.genes.faa)

now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | $genenumber genes were predicted" >> ${workdir}/run_${timestamp}.log

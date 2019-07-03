#Source settings file
source $settingsfile

#Create DuplicatesRemoved directory
mkdir -p ${workdir}/DuplicatesRemoved

#Select source folder from which data will be retrieved (check if directories contain files)
if [[ "$(ls -A ${workdir}/QualityFiltered/)" ]]; then
sourcefolder="QualityFiltered"
else
sourcefolder="RawData"
fi
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Removing duplicates from files in directory $sourcefolder" >> ${workdir}/run_${timestamp}.log

#Declare function
function remdupjob() {

sample=${1}
settingsfile=${2}
sourcefolder=${3}

source $settingsfile

#Obtain data from sample.data.txt columns
samplename=$(echo $sample | cut -d ' ' -f1 )
sampleinfo=$(echo $sample | cut -d ' ' -f2 )

now=$(date +"%Y-%m-%d %H:%M:%S")

if [[ $sampleinfo =~ "/" ]]; then
  #It is PE
  #Decompress if needed
  if [[ -f ${workdir}/${sourcefolder}/${samplename}_1.fastq.gz ]]; then
  gunzip ${workdir}/${sourcefolder}/${samplename}_1.fastq.gz
  fi
  if [[ -f ${workdir}/${sourcefolder}/${samplename}_2.fastq.gz ]]; then
  gunzip ${workdir}/${sourcefolder}/${samplename}_2.fastq.gz
  fi
  echo "$now | 		Removing duplicates from sample $samplename" >> ${workdir}/run_${timestamp}.log
  cat ${workdir}/${sourcefolder}/${samplename}_1.fastq | seqkit rmdup -s -d ${workdir}/DuplicatesRemoved/${samplename}_1.duplicates.fastq -o ${workdir}/DuplicatesRemoved/${samplename}_1.fastq #2>> ${workdir}/run_${timestamp}.log
  cat ${workdir}/${sourcefolder}/${samplename}_2.fastq | seqkit rmdup -s -d ${workdir}/DuplicatesRemoved/${samplename}_2.duplicates.fastq -o ${workdir}/DuplicatesRemoved/${samplename}_2.fastq #2>> ${workdir}/run_${timestamp}.log
  #Repair paired-end reads using BBMap script repair.sh
  repair.sh in=${workdir}/DuplicatesRemoved/${samplename}_1.fastq in2=${workdir}/${sourcefolder}/${samplename}_2.fastq out=${workdir}/DuplicatesRemoved/${samplename}_1.fastq out2=${workdir}/DuplicatesRemoved/${samplename}_2.fastq overwrite=t
  #Get statistics
  before1_1=$(cat ${workdir}/${sourcefolder}/${samplename}_1.fastq | wc -l)
  before1_2=$((before1_1 / 4))
  after1_1=$(cat ${workdir}/DuplicatesRemoved/${samplename}_1.fastq | wc -l)
  after1_2=$((after1_1 / 4))
  difference1=$((before1_2 - after1_2))
  percentage1=$((100-(after1_2 * 100 / before1_2 )))
  #Print statistics
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 		From sample $samplename, $difference1 duplicated reads out of $before1_2 (${percentage1}%) were removed " >> ${workdir}/run_${timestamp}.log
  #Compress source files
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 		Compressing files ${sourcefolder}/${samplename}_1.fastq and ${sourcefolder}/${samplename}_2.fastq" >> ${workdir}/run_${timestamp}.log
  pigz -p ${threads} ${workdir}/${sourcefolder}/${samplename}_1.fastq
  pigz -p ${threads} ${workdir}/${sourcefolder}/${samplename}_2.fastq
else
  #It is SR
  #Decompress if needed
  if [[ -f ${workdir}/${sourcefolder}/${samplename}.fastq.gz ]]; then
  gunzip ${workdir}/${sourcefolder}/${samplename}.fastq.gz
  fi
  echo "$now | 		Removing duplicates from sample $samplename" >> ${workdir}/run_${timestamp}.log
  cat ${workdir}/${sourcefolder}/${samplename}.fastq | seqkit rmdup -s -o ${workdir}/DuplicatesRemoved/${samplename}.fastq
  #Get statistics
  before1=$(cat ${workdir}/${sourcefolder}/${samplename}.fastq | wc -l)
  before2=$((before1 / 4))
  after1=$(cat ${workdir}/DuplicatesRemoved/${samplename}.fastq | wc -l)
  after2=$((after1 / 4))
  difference=$((before2 - after2))
  percentage=$((100-(after2 * 100 / before2 )))
  #Print statistics
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 		From sample $samplename, $difference duplicated reads (${percentage}%) were removed " >> ${workdir}/run_${timestamp}.log
  #Compress source file
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 		Compressing file ${sourcefolder}/${samplename}.fastq" >> ${workdir}/run_${timestamp}.log
  pigz -p ${threads} ${workdir}/${sourcefolder}/${samplename}.fastq
fi
}

#Loop in parallel across samples specified in sample.data.txt
export -f remdupjob
parallel -j ${threads} -k remdupjob {} ${settingsfile} ${sourcefolder} <${sampledatafile}

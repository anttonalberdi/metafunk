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
  echo "$now | 		Removing duplicates from sample $samplename" >> ${workdir}/run_${timestamp}.log
  cat ${workdir}/${sourcefolder}/${samplename}_1.fastq | seqkit rmdup -s -o ${workdir}/DuplicatesRemoved/${samplename}_1.fastq
  cat ${workdir}/${sourcefolder}/${samplename}_2.fastq | seqkit rmdup -s -o ${workdir}/DuplicatesRemoved/${samplename}_2.fastq
  #Get statistics
  before1_1=$(cat ${workdir}/${sourcefolder}/${samplename}_1.fastq | wc -l)
  before1_2=$((before1_1 / 4))
  before2_1=$(cat ${workdir}/${sourcefolder}/${samplename}_2.fastq | wc -l)
  before2_2=$((before2_1 / 4))
  after1_1=$(cat ${workdir}/DuplicatesRemoved/${samplename}_1.fastq | wc -l)
  after1_2=$((after1_1 / 4))
  after2_1=$(cat ${workdir}/DuplicatesRemoved/${samplename}_2.fastq | wc -l)
  after2_2=$((after2_1 / 4))
  difference1=$((before1_2 - after1_2))
  difference2=$((before2_2 - after2_2))
  percentage1=$((100-(after1_2 * 100 / before1_2 )))
  percentage2=$((100-(after2_2 * 100 / before2_2 )))
  #Print statistics
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | 		From sample $samplename, $difference1 (PE1) and $difference2 (PE2) duplicated reads (${percentage1}% and ${percentage2}%) were removed " >> ${workdir}/run_${timestamp}.log
else
  #It is SR
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
fi
}

#Loop in parallel across samples specified in sample.data.txt
export -f remdupjob
parallel -j ${threads} -k remdupjob {} ${settingsfile} ${sourcefolder} <${sampledatafile}

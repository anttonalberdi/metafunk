#Source dependencies
source "$metafunkdirectory/settings.sh"

#Create DuplicatesRemoved directory
mkdir -p ${workingdirectory}/${project}/DuplicatesRemoved

#Select source folder from which data will be retrieved
if [[ "$(ls -A ${workingdirectory}/${project}/QualityFiltered/)" ]]; then
sourcefolder="QualityFiltered"
else
sourcefolder="RawData"
fi
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Duplicates will be removed from files in directory $sourcefolder" >> ${workingdirectory}/${project}/run.log

#Loop across samples specified in sample.data.txt
while read sample; do

  #Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  samplefile=$(echo $sample | cut -d ' ' -f2 )

  now=$(date +"%Y-%d-%m %H:%M:%S")

  if [[ $samplefile =~ "/" ]]; then
  #It is PE
  echo "$now | 		Removing duplicates from sample $samplename" >> ${workingdirectory}/${project}/run.log
  cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}_1.fastq | seqkit rmdup -s -o ${workingdirectory}/${project}/DuplicatesRemoved/${samplename}_1.fastq
  cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}_2.fastq | seqkit rmdup -s -o ${workingdirectory}/${project}/DuplicatesRemoved/${samplename}_2.fastq
  #Get statistics
  before1_1=$(cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}_1.fastq | wc -l)
  before1_2=$((before1_1 / 4))
  before2_1=$(cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}_2.fastq | wc -l)
  before2_2=$((before2_1 / 4))
  after1_1=$(cat ${workingdirectory}/${project}/DuplicatesRemoved/${samplename}_1.fastq | wc -l)
  after1_2=$((after1_1 / 4))
  after2_1=$(cat ${workingdirectory}/${project}/DuplicatesRemoved/${samplename}_2.fastq | wc -l)
  after2_2=$((after2_1 / 4))
  difference1=$((before1_2 - after1_2))
  difference2=$((before2_2 - after2_2))
  percentage1=$((100-(after1_2 * 100 / before1_2 )))
  percentage2=$((100-(after2_2 * 100 / before2_2 )))
  #Print statistics
  now=$(date +"%Y-%d-%m %H:%M:%S")
  echo "$now | 		From sample $samplename, $difference1 (PE1) and $difference2 (PE2) duplicated reads (${percentage1}% and ${percentage2}%) were removed " >> ${workingdirectory}/${project}/run.log
  else
  #It is SR
  echo "$now | 		Removing duplicates from sample $samplename" >> ${workingdirectory}/${project}/run.log
  cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}.fastq | seqkit rmdup -s -o ${workingdirectory}/${project}/DuplicatesRemoved/${samplename}.fastq
  #Get statistics
  before1=$(cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}.fastq | wc -l)
  before2=$((before1 / 4))
  after1=$(cat ${workingdirectory}/${project}/DuplicatesRemoved/${samplename}.fastq | wc -l)
  after2=$((after1 / 4))
  difference=$((before2 - after2))
  percentage=$((100-(after2 * 100 / before2 )))
  #Print statistics
  now=$(date +"%Y-%d-%m %H:%M:%S")
  echo "$now | 		From sample $samplename, $difference duplicated reads (${percentage}%) were removed " >> ${workingdirectory}/${project}/run.log

  fi




done < ${metafunkdirectory}/sample.data.txt

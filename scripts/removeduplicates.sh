#Source dependencies
source "$metafunkdirectory/settings.sh"

mkdir -p ${workingdirectory}/${project}/DuplicatesRemoved

#Select source folder from which data will be retrieved
if [[ $qualityfiltering == "yes" ]]; then
sourcefolder="QualityFiltered"
else
sourcefolder="RawData"
fi

#Loop across samples specified in sample.data.txt
while read sample; do

  #Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  sampleread=$(echo $sample | cut -d ' ' -f2 )

  #Get file names
  if [[ $seqtype == "SR" ]]; then
  samplefile=${samplename}; fi
  if [[ $seqtype == "PE" ]]; then
  samplefile=$(echo ${samplename}_${sampleread}); fi

  #Do the actual duplicate removal
  now=$(date +"%Y-%d-%m %H:%M:%S")
  echo "$now | 		Removing duplicates from sample $samplefile" >> ${workingdirectory}/${project}/run.log
  cat ${workingdirectory}/${project}/RawData/${samplefile}.fastq | seqkit rmdup -s -o ${workingdirectory}/${project}/DuplicatesRemoved/${samplefile}.fastq

  #Get statistics
  before1=$(cat ${workingdirectory}/${project}/RawData/${samplefile}.fastq | wc -l)
  before2=$((before1 / 4))
  after1=$(cat ${workingdirectory}/${project}/DuplicatesRemoved/${samplefile}.fastq | wc -l)
  after2=$((after1 / 4))
  difference=$((before2 - after2))
  percentage=$((100-(after2 * 100 / before2 )))

  #Print statistics
  now=$(date +"%Y-%d-%m %H:%M:%S")
  echo "$now | 		$difference duplicated reads ($percentage %) were removed from sample $samplefile" >> ${workingdirectory}/${project}/run.log
done < ${metafunkdirectory}/sample.data.txt

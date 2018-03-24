#Source dependencies
source "$metafunkdirectory/settings.sh"

#Create LowComplexFiltered directory
mkdir -p ${workingdirectory}/${project}/LowComplexFiltered

#Select source folder from which data will be retrieved
if [[ $qualityfiltering == "yes" ]]; then
sourcefolder="QualityFiltered"
elif [[ $removeduplicates == "yes" ]]; then
sourcefolder="DuplicatesRemoved"
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

  #Perform the actual low complexity filtering
  now=$(date +"%Y-%d-%m %H:%M:%S")
  echo "$now | 		Filtering sample ${samplefile}" >> ${workingdirectory}/${project}/run.log
  prinseq-lite.pl -lc_method "dust" -lc_threshold ${dustvalue} -fastq  ${workingdirectory}/${project}/${sourcefolder}/${samplefile}.fastq -out_good ${workingdirectory}/${project}/LowComplexFiltered/${samplefile}
  now=$(date +"%Y-%d-%m %H:%M:%S")
  echo "$now | 		${samplefile} was successfully filtered" >> ${workingdirectory}/${project}/run.log

done < ${metafunkdirectory}/sample.data.txt

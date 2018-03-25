#Source dependencies
source "$metafunkdirectory/settings.sh"

mkdir -p ${workingdirectory}/${project}/QualityFiltered

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

now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Quality filtering sample $samplefile" >> ${workingdirectory}/${project}/run.log
AdapterRemoval --file1 ${workingdirectory}/${project}/RawData/${samplefile}.fastq --basename ${workingdirectory}/${project}/QualityFiltered/${samplefile} --minquality ${minavgquality} --minlength ${minseqlength} --trimqualities --trimns --maxns 5 --qualitymax ${qualitymax} --threads ${threads}
mv ${workingdirectory}/${project}/QualityFiltered/${samplefile}.truncated ${workingdirectory}/${project}/QualityFiltered/${samplefile}.fastq
now=$(date +"%Y-%d-%m %H:%M:%S")
#MESSAGE TO MYSELF: ASS STATISTICS!!!
echo "$now | 		Sample $samplefile was succesfully filtered" >> ${workingdirectory}/${project}/run.log
done < ${metafunkdirectory}/sample.data.txt

#Source settings
source "$metafunkdirectory/settings.sh"

#Create QualityFiltered directory
mkdir -p ${workingdirectory}/${project}/QualityFiltered

#Loop across samples specified in sample.data.txt
while read sample; do

  #Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  samplefile=$(echo $sample | cut -d ' ' -f2 )
  now=$(date +"%Y-%d-%m %H:%M:%S")

  if [[ $samplefile =~ "/" ]]; then
  #It is PE
  echo "$now | 		Quality filtering sample $samplename" >> ${workingdirectory}/${project}/run.log
  AdapterRemoval --file1 ${workingdirectory}/${project}/RawData/${samplename}_1.fastq --file2 ${workingdirectory}/${project}/RawData/${samplename}_2.fastq --basename ${workingdirectory}/${project}/QualityFiltered/${samplename} --minquality ${minavgquality} --minlength ${minseqlength} --trimqualities --trimns --maxns 5 --qualitymax ${qualitymax} --threads ${threads}
  mv ${workingdirectory}/${project}/QualityFiltered/${samplename}.pair1.truncated ${workingdirectory}/${project}/QualityFiltered/${samplename}_1.fastq
  mv ${workingdirectory}/${project}/QualityFiltered/${samplename}.pair2.truncated ${workingdirectory}/${project}/QualityFiltered/${samplename}_2.fastq
  else
  #It is SR
  echo "$now | 		Quality filtering sample $samplename" >> ${workingdirectory}/${project}/run.log
  AdapterRemoval --file1 ${workingdirectory}/${project}/RawData/${samplename}.fastq --basename ${workingdirectory}/${project}/QualityFiltered/${samplename} --minquality ${minavgquality} --minlength ${minseqlength} --trimqualities --trimns --maxns 5 --qualitymax ${qualitymax} --threads ${threads}
  mv ${workingdirectory}/${project}/QualityFiltered/${samplename}.truncated ${workingdirectory}/${project}/QualityFiltered/${samplename}.fastq
  fi

now=$(date +"%Y-%d-%m %H:%M:%S")
#MESSAGE TO MYSELF: ASS STATISTICS!!!
echo "$now | 		Sample $samplename was succesfully filtered" >> ${workingdirectory}/${project}/run.log

done < ${metafunkdirectory}/sample.data.txt

#Source settings
source "$metafunkdirectory/settings.sh"

#Create QualityFiltered directory
mkdir -p ${workingdirectory}/${project}/QualityFiltered

#Loop across samples specified in sample.data.txt
while read sample; do

  #Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  sampleread=$(echo $sample | cut -d ' ' -f2 )
  now=$(date +"%Y-%d-%m %H:%M:%S")

  #Identify PE or SR
  if [[ $samplefile =~ "/" ]]; then
  #It is PE
  echo "$now | 		Quality filtering sample $samplefile" >> ${workingdirectory}/${project}/run.log
  AdapterRemoval --file1 ${workingdirectory}/${project}/RawData/${samplefile}_1.fastq --file2 ${workingdirectory}/${project}/RawData/${samplefile}_2.fastq --basename ${workingdirectory}/${project}/QualityFiltered/${samplefile} --minquality ${minavgquality} --minlength ${minseqlength} --trimqualities --trimns --maxns 5 --qualitymax ${qualitymax} --threads ${threads}
  mv ${workingdirectory}/${project}/QualityFiltered/${samplefile}.pair1.truncated ${workingdirectory}/${project}/QualityFiltered/${samplefile}_1.fastq
  mv ${workingdirectory}/${project}/QualityFiltered/${samplefile}.pair2.truncated ${workingdirectory}/${project}/QualityFiltered/${samplefile}_2.fastq
  else
  #It is SR
  echo "$now | 		Quality filtering sample $samplefile" >> ${workingdirectory}/${project}/run.log
  AdapterRemoval --file1 ${workingdirectory}/${project}/RawData/${samplefile}.fastq --basename ${workingdirectory}/${project}/QualityFiltered/${samplefile} --minquality ${minavgquality} --minlength ${minseqlength} --trimqualities --trimns --maxns 5 --qualitymax ${qualitymax} --threads ${threads}
  mv ${workingdirectory}/${project}/QualityFiltered/${samplefile}.truncated ${workingdirectory}/${project}/QualityFiltered/${samplefile}.fastq
  fi

now=$(date +"%Y-%d-%m %H:%M:%S")
#MESSAGE TO MYSELF: ASS STATISTICS!!!
echo "$now | 		Sample $samplefile was succesfully filtered" >> ${workingdirectory}/${project}/run.log

done < ${metafunkdirectory}/sample.data.txt

source "$metafunkdirectory/settings.sh"

mkdir -p ${workingdirectory}/${project}/RawData

#Loop across samples specified in sample.data.txt
while read sample; do

  #Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  sampleread=$(echo $sample | cut -d ' ' -f2 )
  samplefile=$(echo $sample | cut -d ' ' -f3 )

  #If uncompressed SR files
  if [[ $compressed == "no" && $seqtype == "SR" ]]; then
    if [ ! -f "${datadirectory}/${samplefile}" ]; then
      echo "ERROR: File ${samplefile} was not found. Check the settings are correct."
    else
      echo "Transferring file ${samplefile}"
    cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/
    mv ${workingdirectory}/${project}/RawData/${samplename}.* ${workingdirectory}/${project}/RawData/${samplename}.fastq
    fi
  fi

  #If uncompressed PE files
  if [[ $compressed == "no" && $seqtype == "PE" ]]; then
    if [ ! -f "${datadirectory}/${samplefile}" ]; then
      echo "File ${samplefile} was not found. Check the settings are correct."
    else
      echo "Transferring file ${samplefile}"
      if [[ $sampleread == 1 ]]; then
      cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}_1.fastq
      fi
      if [[ $sampleread == 2 ]]; then
      cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}_2.fastq
      fi
    fi
  fi

  #If compressed SR files
  if [[ $compressed == "yes" && $seqtype == "SR" ]]; then
    if [ ! -f "${datadirectory}/${samplefile}" ]; then
      echo "File ${samplefile} was not found. Check the settings are correct."
    else
      echo "Transferring and uncompressing file ${samplefile}"
    cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}.fastq
    pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}.fastq
    fi
  fi

  #If compressed PE files
  if [[ $compressed == "yes" && $seqtype == "SR" ]]; then
    if [ ! -f "${datadirectory}/${samplefile}" ]; then
      echo "File ${samplefile} was not found. Check the settings are correct."
    else
      echo "Transferring and uncompressing file ${samplefile}"
      #Process if read 1
      if [[ $sampleread == 1 ]]; then
      cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}_1.fastq.gz
      pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}_1.fastq.gz
      fi
      #Process if read 2
      if [[ $sampleread == 2 ]]; then
      cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}_2.fastq.gz
      pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}_2.fastq.gz
      fi
    fi
  fi
done < ${metafunkdirectory}/sample.data.txt

#Check if files were succesfully transferred
if [ -z "$(ls -A ${workingdirectory}/${project})" ]; then
  echo "ERROR: The data were not transferred."
  exit
else
  #Print stats
  filenumber=$(ls ${workingdirectory}/${project}/RawData/| wc -l)
  now=$(date +"%Y-%d-%m %H:%M:%S")
  echo "$now | $filenumber files were copied and uncompressed" >> ${workingdirectory}/${project}/run.log

fi

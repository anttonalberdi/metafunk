
  now=$(date +"%Y-%d-%m %H:%M:%S")

  #If uncompressed SR files
  if [[ $compressed == "no" && $seqtype == "SR" ]]; then
    if [ ! -f "${datadirectory}/${samplefile}" ]; then
      echo "$now |    ERROR: File ${samplefile} was not found. Check the settings are correct"  >> ${workingdirectory}/${project}/run.log
    else
      echo "$now |    Transferring file ${samplefile}"  >> ${workingdirectory}/${project}/run.log
    cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/
    mv ${workingdirectory}/${project}/RawData/${samplename}.* ${workingdirectory}/${project}/RawData/${samplename}.fastq
    fi
  fi

  #If uncompressed PE files
  if [[ $compressed == "no" && $seqtype == "PE" ]]; then
    if [ ! -f "${datadirectory}/${samplefile}" ]; then
      echo "$now |    ERROR: File ${samplefile} was not found. Check the settings are correct"  >> ${workingdirectory}/${project}/run.log
    else
      echo "$now |    Transferring file ${samplefile}"  >> ${workingdirectory}/${project}/run.log
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
      echo "$now |    ERROR: File ${samplefile} was not found. Check the settings are correct"  >> ${workingdirectory}/${project}/run.log
    else
      echo "$now |    Transferring and uncompressing file ${samplefile}" >> ${workingdirectory}/${project}/run.log
    cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}.fastq
    pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}.fastq
    fi
  fi

  #If compressed PE files
  if [[ $compressed == "yes" && $seqtype == "PE" ]]; then
    if [ ! -f "${datadirectory}/${samplefile}" ]; then
      echo "$now |    ERROR: File ${samplefile} was not found. Check the settings are correct" >> ${workingdirectory}/${project}/run.log
    else
      echo "$now |    Transferring and uncompressing file ${samplefile}" >> ${workingdirectory}/${project}/run.log
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

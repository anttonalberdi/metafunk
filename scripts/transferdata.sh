#Source settings file
source $settingsfile

###### NOTE FOR MYSELF: SR multifile IS NOT WORKING PROPERLY!! 26/03/2018

mkdir -p ${workdir}/RawData

#Loop across samples specified in sample.data.txt
while read sample; do

  #Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  sampleinfo=$(echo $sample | cut -d ' ' -f2 )
  now=$(date +"%Y-%d-%m %H:%M:%S")

  echo "$now |    Processing sample $samplename" >> ${workdir}/run_${timestamp}.log

  if [[ $sampleinfo =~ "/" && ! $sampleinfo =~ ";" ]]; then
  #It is PE single file
    #Get file names
    samplefile1=$(echo $sampleinfo | cut -d'/' -f1)
    samplefile2=$(echo $sampleinfo | cut -d'/' -f2)

    #Transfer both files
      #PE1
      if [[ $samplefile1 == *.fastq.gz || $samplefile1 == *.fq.gz ]]; then
      cp ${datadir}/${samplefile1} ${workdir}/RawData/${samplename}_1.fastq.gz
      pigz -d -p ${threads} ${workdir}/RawData/${samplename}_1.fastq.gz
      elif [[ $samplefile1 == *.fastq || $samplefile1 == *.fq ]]; then
      cp ${datadir}/${samplefile1} ${workdir}/RawData/${samplename}_1.fastq.gz
      else
      echo "$now |    ERROR: The extension of file $samplefile1 is not recognised" >> ${workdir}/run_${timestamp}.log
      fi

      #PE2
      if [[ $samplefile2 == *.fastq.gz || $samplefile2 == *.fq.gz ]]; then
      cp ${datadir}/${samplefile2} ${workdir}/RawData/${samplename}_2.fastq.gz
      pigz -d -p ${threads} ${workdir}/RawData/${samplename}_2.fastq.gz
      elif [[ $samplefile1 == *.fastq || $samplefile1 == *.fq ]]; then
      cp ${datadir}/${samplefile1} ${workdir}/RawData/${samplename}_2.fastq.gz
      else
      echo "$now |    ERROR: The extension of file $samplefile2 is not recognised" >> ${workdir}/run_${timestamp}.log
      fi

  elif [[ $sampleinfo =~ "/" && $sampleinfo =~ ";" ]]; then
  #It is PE multi-file
    #Get file names
    samplefile1=$(echo $sampleinfo | cut -d'/' -f1)
    samplefile2=$(echo $sampleinfo | cut -d'/' -f2)

    #PE1
    IFS='; ' read -r -a array <<< $samplefile1
    n=0
    for samplefile in "${array[@]}"; do
      n=$((n+1))
      if [[ $samplefile == *.fastq.gz || $samplefile == *.fq.gz ]]; then
      cp ${datadir}/${samplefile} ${workdir}/RawData/${samplename}_1_${n}.fastq.gz
      pigz -d -p ${threads} ${workdir}/RawData/${samplename}_1_${n}.fastq.gz
      elif [[ $samplefile == *.fastq || $samplefile == *.fq ]]; then
      cp ${datadir}/${samplefile} ${workdir}/RawData/${samplename}_1_${n}.fastq
      else
      echo "$now |    ERROR: The extension of file $samplefile is not recognised" >> ${workdir}/run_${timestamp}.log
      fi
    done
    #Merge all files
    now=$(date +"%Y-%d-%m %H:%M:%S")
    echo "$now |      Merging PE1 files" >> ${workdir}/run_${timestamp}.log
    cat ${workdir}/RawData/${samplename}_1_* > ${workdir}/RawData/${samplename}_1.fastq
    rm ${workdir}/RawData/${samplename}_1_*

    #PE2
    IFS='; ' read -r -a array <<< $samplefile2
    n=0
    for samplefile in "${array[@]}"; do
      n=$((n+1))
      if [[ $samplefile == *.fastq.gz || $samplefile == *.fq.gz ]]; then
      cp ${datadir}/${samplefile} ${workdir}/RawData/${samplename}_2_${n}.fastq.gz
      pigz -d -p ${threads} ${workdir}/RawData/${samplename}_2_${n}.fastq.gz
      elif [[ $samplefile == *.fastq || $samplefile == *.fq ]]; then
      cp ${datadir}/${samplefile} ${workdir}/RawData/${samplename}_2_${n}.fastq
      else
      echo "$now |    ERROR: The extension of file $samplefile is not recognised" >> ${workdir}/run_${timestamp}.log
      fi
    done
    #Merge all files
    echo "$now |      Merging PE2 files" >> ${workdir}/run_${timestamp}.log
    cat ${workdir}/RawData/${samplename}_2_* > ${workdir}/RawData/${samplename}_2.fastq
    rm ${workdir}/RawData/${samplename}_2_*

  elif [[ ! $sampleinfo =~ "/" && $sampleinfo =~ ";" ]]; then
  #It is SR multifile
  #Get file names
  IFS='; ' read -r -a array <<< $sampleinfo
  n=0
  for samplefile in "${array[@]}"; do
    n=$((n+1))
    if [[ $samplefile == *.fastq.gz || $samplefile == *.fq.gz ]]; then
    cp ${datadir}/${samplefile} ${workdir}/RawData/${samplename}_${n}.fastq.gz
    pigz -d -p ${threads} ${workdir}/RawData/${samplename}_${n}.fastq.gz
    elif [[ $samplefile == *.fastq || $samplefile == *.fq ]]; then
    cp ${datadir}/${samplefile} ${workdir}/RawData/${samplename}_${n}.fastq
    else
    echo "$now |    The extension of file $samplefile is not recognised" >> ${workdir}/run_${timestamp}.log
    fi
    #Merge all files
    now=$(date +"%Y-%d-%m %H:%M:%S")
    echo "$now |    Merging $samplename files" >> ${workdir}/run_${timestamp}.log
    cat ${workdir}/RawData/${samplename}_* > ${workdir}/RawData/${samplename}.fastq
    rm ${workdir}/RawData/${samplename}_*
  done

  else
  #It is SR single file
    if [[ $sampleinfo == *.fastq.gz || $sampleinfo == *.fq.gz ]]; then
    cp ${datadir}/${samplefile} ${workdir}/RawData/${samplename}.fastq.gz
    pigz -d -p ${threads} ${workdir}/RawData/${samplename}.fastq.gz
    elif [[ $sampleinfo == *.fastq || $sampleinfo == *.fq ]]; then
    cp ${datadir}/${samplefile} ${workdir}/RawData/${samplename}.fastq
    else
    echo "$now |    ERROR: The extension of file $samplefile is not recognised" >> ${workdir}/run_${timestamp}.log
    fi
  fi
done < ${sampledatafile}

#Check if files were succesfully transferred
if [ -z "$(ls -A ${workdir})" ]; then
  now=$(date +"%Y-%d-%m %H:%M:%S")
  echo "$now |    ERROR: The data were not transferred"  >> ${workdir}/run_${timestamp}.log
  exit
else
  #Print stats
  filenumber=$(ls ${workdir}/RawData/| wc -l)
  now=$(date +"%Y-%d-%m %H:%M:%S")
  echo "$now |    $filenumber files belonging to $samplenumber samples were succesfully processed" >> ${workdir}/run_${timestamp}.log
fi

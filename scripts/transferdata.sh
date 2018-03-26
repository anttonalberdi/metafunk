source "$metafunkdirectory/settings.sh"

###### NOTE FOR MYSELF: SR multifile MIGHT NOT BE WORKING PROPERLY!! 26/03/2018

mkdir -p ${workingdirectory}/${project}/RawData

#Loop across samples specified in sample.data.txt
while read sample; do

  #Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  samplefile=$(echo $sample | cut -d ' ' -f2 )
  now=$(date +"%Y-%d-%m %H:%M:%S")

  if [[ $samplefile =~ "/" && ! $samplefile =~ ";" ]]; then
  #It is PE single file
    echo "$now |    Transferring PE sample $samplename" >> ${workingdirectory}/${project}/run.log
    #Get file names
    samplefile1=$(echo $samplefile | cut -d'/' -f1)
    samplefile2=$(echo $samplefile | cut -d'/' -f2)

    #Transfer both files
      #PE1
      if [[ $samplefile1 == *.fastq.gz || $samplefile1 == *.fq.gz ]]; then
      cp ${datadirectory}/${samplefile1} ${workingdirectory}/${project}/RawData/${samplename}_1.fastq.gz
      pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}_1.fastq.gz
      elif [[ $samplefile1 == *.fastq || $samplefile1 == *.fq ]]; then
      cp ${datadirectory}/${samplefile1} ${workingdirectory}/${project}/RawData/${samplename}_1.fastq.gz
      else
      echo "$now |    ERROR: The extension of file $samplefile1 is not recognised" >> ${workingdirectory}/${project}/run.log
      fi

      #PE2
      if [[ $samplefile2 == *.fastq.gz || $samplefile2 == *.fq.gz ]]; then
      cp ${datadirectory}/${samplefile2} ${workingdirectory}/${project}/RawData/${samplename}_2.fastq.gz
      pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}_2.fastq.gz
      elif [[ $samplefile1 == *.fastq || $samplefile1 == *.fq ]]; then
      cp ${datadirectory}/${samplefile1} ${workingdirectory}/${project}/RawData/${samplename}_2.fastq.gz
      else
      echo "$now |    ERROR: The extension of file $samplefile2 is not recognised" >> ${workingdirectory}/${project}/run.log
      fi

  elif [[ $samplefile =~ "/" && $samplefile =~ ";" ]]; then
  #It is PE multi-file
    echo "$now |    Transferring PE multifile sample $samplename" >> ${workingdirectory}/${project}/run.log
    #Get file names
    samplefile1=$(echo $samplefile | cut -d'/' -f1)
    samplefile2=$(echo $samplefile | cut -d'/' -f2)

    #PE1
    IFS='; ' read -r -a array <<< $samplefile1
    n=0
    for samplefile in "${array[@]}"; do
      n=$((n+1))
      if [[ $samplefile == *.fastq.gz || $samplefile == *.fq.gz ]]; then
      cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}_1_${n}.fastq.gz
      pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}_1_${n}.fastq.gz
      elif [[ $samplefile == *.fastq || $samplefile == *.fq ]]; then
      cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}_1_${n}.fastq
      else
      echo "$now |    ERROR: The extension of file $samplefile is not recognised" >> ${workingdirectory}/${project}/run.log
      fi
    done
    #Merge all files
    cat ${workingdirectory}/${project}/RawData/${samplename}_1_* > ${workingdirectory}/${project}/RawData/${samplename}_1.fastq
    rm ${workingdirectory}/${project}/RawData/${samplename}_1_*

    #PE2
    IFS='; ' read -r -a array <<< $samplefile2
    n=0
    for samplefile in "${array[@]}"; do
      n=$((n+1))
      if [[ $samplefile == *.fastq.gz || $samplefile == *.fq.gz ]]; then
      cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}_2_${n}.fastq.gz
      pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}_2_${n}.fastq.gz
      elif [[ $samplefile == *.fastq || $samplefile == *.fq ]]; then
      cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}_2_${n}.fastq
      else
      echo "$now |    ERROR: The extension of file $samplefile is not recognised" >> ${workingdirectory}/${project}/run.log
      fi
    done
    #Merge all files
    cat ${workingdirectory}/${project}/RawData/${samplename}_2_* > ${workingdirectory}/${project}/RawData/${samplename}_2.fastq
    rm ${workingdirectory}/${project}/RawData/${samplename}_2_*

  elif [[ ! $samplefile =~ "/" && $samplefile =~ ";" ]]; then
  #It is SR multifile
  echo "$now |    Transferring SR multifile sample $samplename" >> ${workingdirectory}/${project}/run.log
  #Get file names
  IFS='; ' read -r -a array <<< $samplefile
  n=0
  for samplefile in "${array[@]}"; do
    n=$((n+1))
    if [[ $samplefile == *.fastq.gz || $samplefile == *.fq.gz ]]; then
    cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}_${n}.fastq.gz
    pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}_${n}.fastq.gz
    elif [[ $samplefile == *.fastq || $samplefile == *.fq ]]; then
    cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}_${n}.fastq
    else
    echo "$now |    The extension of file $samplefile is not recognised" >> ${workingdirectory}/${project}/run.log
    fi
    #Merge all files
    cat ${workingdirectory}/${project}/RawData/${samplename}_* > ${workingdirectory}/${project}/RawData/${samplename}.fastq
    rm ${workingdirectory}/${project}/RawData/${samplename}_*
  done

  else
  #It is SR single file
  echo "$now |    Transferring SR single file sample $samplename" >> ${workingdirectory}/${project}/run.log
    if [[ $samplefile == *.fastq.gz || $samplefile == *.fq.gz ]]; then
    cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}.fastq.gz
    pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}.fastq.gz
    elif [[ $samplefile == *.fastq || $samplefile == *.fq ]]; then
    cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}.fastq
    else
    echo "$now |    ERROR: The extension of file $samplefile is not recognised" >> ${workingdirectory}/${project}/run.log
    fi
  fi
done < ${metafunkdirectory}/sample.data.txt

#Check if files were succesfully transferred
if [ -z "$(ls -A ${workingdirectory}/${project})" ]; then
  echo "ERROR: The data were not transferred"  >> ${workingdirectory}/${project}/run.log
  exit
else
  #Print stats
  filenumber=$(ls ${workingdirectory}/${project}/RawData/| wc -l)
  now=$(date +"%Y-%d-%m %H:%M:%S")
  echo "$now | $filenumber files were processed" >> ${workingdirectory}/${project}/run.log
fi

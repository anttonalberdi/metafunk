source "$metafunkdirectory/settings.sh"

mkdir -p ${workingdirectory}/${project}/RawData

#Loop across samples specified in sample.data.txt
while read sample; do

  #Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  samplefile=$(echo $sample | cut -d ' ' -f2 )

  echo $samplefile

  #Get samples if PE and/or multifile
  if [[ $samplefile =~ "/" && ! $samplefile =~ ";" ]]; then
    #It is PE
    echo "Transferring PE sample $samplename"
    #Get file names
    samplefile1=$(echo $samplefile | cut -d'/' -f1)
    samplefile2=$(echo $samplefile | cut -d'/' -f2)
    #Transfer both files
    cp ${datadirectory}/${samplefile1} ${workingdirectory}/${project}/RawData/${samplename}_1.fastq.gz
    pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}_1.fastq.gz
    cp ${datadirectory}/${samplefile2} ${workingdirectory}/${project}/RawData/${samplename}_2.fastq.gz
    pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}_2.fastq.gz
  elif [[ $samplefile =~ "/" && $samplefile =~ ";" ]]; then
    #It is PE multifile
    echo "Transferring PE multifile sample $samplename"
    #Get file names
    samplefile1=$(echo $samplefile | cut -d'/' -f1)
    samplefile2=$(echo $samplefile | cut -d'/' -f2)
    #Transfer 1st pairs
    IFS='; ' read -r -a array <<< $samplefile1
    n=0
    for samplefile in "${array[@]}"; do
    n=$((n+1))
    cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}_1_${n}.fastq.gz
    pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}_1_${n}.fastq.gz
    done
    #Merge all files
    cat ${workingdirectory}/${project}/RawData/${samplename}_1_* > ${workingdirectory}/${project}/RawData/${samplename}_1.fastq
    rm ${workingdirectory}/${project}/RawData/${samplename}_1_*
    #Transfer 2nd pairs
    IFS='; ' read -r -a array <<< $samplefile2
    n=0
    for samplefile in "${array[@]}"; do
    n=$((n+1))
    cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}_2_${n}.fastq
    pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplename}_2_${n}.fastq.gz
    done
    #Merge all files
    cat ${workingdirectory}/${project}/RawData/${samplename}_2_* > ${workingdirectory}/${project}/RawData/${samplename}_2.fastq
    rm ${workingdirectory}/${project}/RawData/${samplename}_2_*
  else
  #It is SR
  echo "Transferring SR sample $samplename"
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

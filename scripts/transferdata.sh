source "$metafunkdirectory/settings.sh"

mkdir -p ${workingdirectory}/${project}/RawData

while read sample; do
samplename=$(echo $sample | cut -d ' ' -f1 )
samplefile=$(echo $sample | cut -d ' ' -f2 )

#If uncompressed files
if [[ $comp == "no" ]]; then
  if [ ! -f "${datadirectory}/${samplefile}" ]; then
    echo "File ${samplefile} was not found. Check the settings are correct."
  fi
  cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/
  mv ${workingdirectory}/${project}/RawData/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}.fastq
fi

#If compressed files
if [[ $comp != "no" ]]; then
  if [ ! -f "${datadirectory}/${samplefile}" ]; then
    echo "File ${samplefile} was not found. Check the settings are correct."
  fi
cp ${datadirectory}/${samplefile} ${workingdirectory}/${project}/RawData/
pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplefile}
mv ${workingdirectory}/${project}/RawData/${samplefile} ${workingdirectory}/${project}/RawData/${samplename}.fastq
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

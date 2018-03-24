source "$metafunkdirectory/settings.sh"

mkdir -p ${workingdirectory}/${project}/RawData

while read samplefile; do

#If uncompressed files
if [[ $compression == "no" ]]; then
  cp ${datadirectory}/${samplefile}.${extension} ${workingdirectory}/${project}/RawData/
  mv ${workingdirectory}/${project}/RawData/${samplefile}.${extension} ${workingdirectory}/${project}/RawData/${samplefile}.fastq
fi

#If compressed files
if [[ ! $compression == "no" ]]; then
cp ${datadirectory}/${samplefile}.${extension}.${compression} ${workingdirectory}/${project}/RawData/
pigz -d -p ${threads} ${workingdirectory}/${project}/RawData/${samplefile}.${extension}.${compression}
mv ${workingdirectory}/${project}/RawData/${samplefile}.${extension} ${workingdirectory}/${project}/RawData/${samplefile}.fastq
fi

done < ${metafunkdirectory}/sample.data.txt

#Check if files were succesfully transferred
if [ -z "$(ls -A ${workingdirectory}/${project})" ]; then
  echo "ERROR: The data were not transderred."
  exit
else
  #Print stats
  filenumber=$(ls ${workingdirectory}/${project}/RawData/| wc -l)
  now=$(date +"%Y-%d-%m %H:%M:%S")
  echo "$now | $filenumber files were copied and uncompressed" >> ${workingdirectory}/${project}/run.log

fi

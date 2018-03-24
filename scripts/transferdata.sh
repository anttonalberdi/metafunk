source "$metafunkdirectory/settings.sh"

mkdir -p ${projectdirectory}/RawData

while read samplefile; do

#If uncompressed files
if [[ $compression == "no" ]]; then
  cp ${datadirectory}/${samplefile}.${extension} ${projectdirectory}/RawData/
  mv ${projectdirectory}/RawData/${samplefile}.${extension} ${projectdirectory}/RawData/${samplefile}.fastq
fi

#If compressed files
if [[ ! $compression == "no" ]]; then
cp ${datadirectory}/${samplefile}.${extension}.${compression} ${projectdirectory}/RawData/
pigz -d -p ${threads} ${projectdirectory}/RawData/${samplefile}.${samplefile}.${extension}.${compression}
mv ${projectdirectory}/RawData/${samplefile}.${extension} ${projectdirectory}/RawData/${samplefile}.fastq
fi

done < ${metafunkdirectory}/sample.data.txt
#Print stats
filenumber=$(ls ${projectdirectory}/RawData/| wc -l)
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | $filenumber files were copied and uncompressed" >> ${projectdirectory}/run.log

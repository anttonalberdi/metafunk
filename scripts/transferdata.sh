source "$metafunkdirectory/settings.sh"

mkdir -p ${projectdirectory}/RawData

while read samplefile; do
cp ${datadirectory}/${samplefile}.f*q.${compression} ${projectdirectory}/RawData/
#Uncompress files
if [[ ! $compression == "no" ]]; then
pigz -d -p ${threads} ${projectdirectory}/RawData/${samplefile}.${samplefile}.f*q.${compression}
fi
done < ${metafunkdirectory}/sample.data.txt
#Print stats
filenumber=$(ls ${projectdirectory}/RawData/| wc -l)
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | $filenumber files were copied and uncompressed" >> ${projectdirectory}/run.log

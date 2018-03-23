#Source dependencies
source "$metafunkdirectory/settings.sh"

mkdir -p ${workingdirectory}/${project}/QualityFiltered

while read samplefile; do
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Quality filtering sample $samplefile" >> ${workingdirectory}/${project}/run.log
AdapterRemoval --file1 ${workingdirectory}/${project}/RawData/${samplefile}.fastq --basename ${workingdirectory}/${project}/QualityFiltered/${samplefile} --minquality 30 --minlength 30 --trimqualities --trimns --maxns 5 --qualitymax 60 --threads ${threads}
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Sample $samplefile was succesfully filtered" >> ${workingdirectory}/${project}/run.log
done < ${metafunkdirectory}/sample.data.txt

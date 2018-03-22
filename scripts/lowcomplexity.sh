#Source dependencies
source "$metafunkdirectory/settings.sh"

#Create LowComplexFiltered directory
mkdir -p ${workingdirectory}/${project}/LowComplexFiltered

#Filter low complexity to all samples
while read samplefile; do 
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Filtering sample ${samplefile}" >> ${workingdirectory}/${project}/run.log
prinseq-lite.pl -lc_method "dust" -lc_threshold ${dustvalue} -fastq  ${workingdirectory}/${project}/RawData/${samplefile}.fastq -out_good ${workingdirectory}/${project}/LowComplexFiltered/${samplefile}.fastq -out_bad ${workingdirectory}/${project}/LowComplexFiltered/${samplefile}.lowQual
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		${samplefile} was successfully filtered" >> ${workingdirectory}/${project}/run.log
done < ${metafunkdirectory}/sample.data.txt

#Source dependencies
source "$metafunkdirectory/settings.sh"

mkdir -p ${workingdirectory}/${project}/RemoveDuplicates

#Iterate through samples
while read samplefile; do
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Removing duplicates from sample $samplefile" >> ${workingdirectory}/${project}/run.log
cat ${workingdirectory}/${project}/RawData/${samplefile}.fastq | seqkit rmdup -s -o ${workingdirectory}/${project}/RemoveDuplicates/${samplefile}.fastq
#Get statistics
before1=$(cat ${workingdirectory}/${project}/RawData/${samplefile}.fastq | wc -l)
before2=$((before1 / 4))
after1=$(cat ${workingdirectory}/${project}/RemoveDuplicates/${samplefile}.fastq | wc -l)
after2=$((after1 / 4))
percentage=$((after2 / before2 * 100))
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		$after2 duplicated reads ($percentage) were removed from sample $samplefile" >> ${workingdirectory}/${project}/run.log
done < ${metafunkdirectory}/sample.data.txt

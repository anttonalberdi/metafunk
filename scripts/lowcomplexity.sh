#Source dependencies
source "$metafunkdirectory/settings.sh"

#Create LowComplexFiltered directory
mkdir -p ${workingdirectory}/${project}/LowComplexFiltered

#Select source folder from which data will be retrieved (check if directories contain files)
if [[ "$(ls -A ${workingdirectory}/${project}/DuplicatesRemoved/)" ]]; then
sourcefolder="DuplicatesRemoved"
elif [[ "$(ls -A ${workingdirectory}/${project}/QualityFiltered/)" ]]; then
sourcefolder="QualityFiltered"
else
sourcefolder="RawData"
fi
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Removing low complexity reads from files in directory ${sourcefolder}" >> ${workingdirectory}/${project}/run.log

#Declare function
function lowcompjob() {

  sample=${1}
  metafunkdirectory=${2}
  sourcefolder=${3}

  source "$metafunkdirectory/settings.sh"

  #Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  samplefile=$(echo $sample | cut -d ' ' -f2 )

  now=$(date +"%Y-%d-%m %H:%M:%S")

  if [[ $samplefile =~ "/" ]]; then
    #It is PE
    echo "$now | 		Removing low complexity reads from PE sample ${samplename}" >> ${workingdirectory}/${project}/run.log
    prinseq-lite.pl -lc_method "dust" -lc_threshold ${dustvalue} -fastq  ${workingdirectory}/${project}/${sourcefolder}/${samplename}_1.fastq -out_good ${workingdirectory}/${project}/LowComplexFiltered/${samplename}_1
    prinseq-lite.pl -lc_method "dust" -lc_threshold ${dustvalue} -fastq  ${workingdirectory}/${project}/${sourcefolder}/${samplename}_2.fastq -out_good ${workingdirectory}/${project}/LowComplexFiltered/${samplename}_2
    #Compute statistics
    before1_1=$(cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}_1.fastq | wc -l)
    before1_2=$((before1_1 / 4))
    before2_1=$(cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}_2.fastq | wc -l)
    before2_2=$((before2_1 / 4))
    after1_1=$(cat ${workingdirectory}/${project}/LowComplexFiltered/${samplename}_1.fastq | wc -l)
    after1_2=$((after1_1 / 4))
    after2_1=$(cat ${workingdirectory}/${project}/LowComplexFiltered/${samplename}_2.fastq | wc -l)
    after2_2=$((after2_1 / 4))
    difference1=$((before1_2 - after1_2))
    difference2=$((before2_2 - after2_2))
    percentage1=$((100-(after1_2 * 100 / before1_2 )))
    percentage2=$((100-(after2_2 * 100 / before2_2 )))
    #Print statistics
    now=$(date +"%Y-%d-%m %H:%M:%S")
    echo "$now | 		From sample $samplename, $difference1 (PE1) and $difference2 (PE2) reads (${percentage1}% and ${percentage2}%) were removed due to low complexity" >> ${workingdirectory}/${project}/run.log
  else
    #It is SR
    echo "$now | 		Removing low complexity reads from SR sample ${samplename}" >> ${workingdirectory}/${project}/run.log
    prinseq-lite.pl -lc_method "dust" -lc_threshold ${dustvalue} -fastq  ${workingdirectory}/${project}/${sourcefolder}/${samplename}.fastq -out_good ${workingdirectory}/${project}/LowComplexFiltered/${samplename}
    #Compute statistics
    before1=$(cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}.fastq | wc -l)
    before2=$((before1 / 4))
    after1=$(cat ${workingdirectory}/${project}/LowComplexFiltered/${samplename}.fastq | wc -l)
    after2=$((after1 / 4))
    difference=$((before2 - after2))
    percentage=$((100-(after2 * 100 / before2 )))
    #Print statistics
    now=$(date +"%Y-%d-%m %H:%M:%S")
    echo "$now | 		From sample $samplename, $difference reads (${percentage}%) were removed due to low complexity" >> ${workingdirectory}/${project}/run.log
  fi
}

#Loop in parallel across samples specified in sample.data.txt
export -f lowcompjob
parallel -j ${threads} -k lowcompjob {} ${metafunkdirectory} ${sourcefolder} <${sampledatafile}

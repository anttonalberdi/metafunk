#Source settings file
source $settingsfile

#Create LowComplexFiltered directory
mkdir -p ${workdir}/LowComplexFiltered

#Select source folder from which data will be retrieved (check if directories contain files)
if [[ "$(ls -A ${workdir}/DuplicatesRemoved/)" ]]; then
sourcefolder="DuplicatesRemoved"
elif [[ "$(ls -A ${workdir}/QualityFiltered/)" ]]; then
sourcefolder="QualityFiltered"
else
sourcefolder="RawData"
fi
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Removing low complexity reads from files in directory ${sourcefolder}" >> ${workdir}/run_${timestamp}.log

#Declare function
function lowcompjob() {

  sample=${1}
  settingsfile=${2}
  sourcefolder=${3}

  source $settingsfile

  #Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  sampleinfo=$(echo $sample | cut -d ' ' -f2 )

  now=$(date +"%Y-%m-%d %H:%M:%S")

  if [[ $sampleinfo =~ "/" ]]; then
    #It is PE
    echo "$now | 		Removing low complexity reads from PE sample ${samplename}" >> ${workdir}/run_${timestamp}.log
    prinseq-lite.pl -lc_method "dust" -lc_threshold ${dustvalue} -fastq  ${workdir}/${sourcefolder}/${samplename}_1.fastq -out_good ${workdir}/LowComplexFiltered/${samplename}_1 -out_bad null
    prinseq-lite.pl -lc_method "dust" -lc_threshold ${dustvalue} -fastq  ${workdir}/${sourcefolder}/${samplename}_2.fastq -out_good ${workdir}/LowComplexFiltered/${samplename}_2 -out_bad null
    #Compute statistics
    before1_1=$(cat ${workdir}/${sourcefolder}/${samplename}_1.fastq | wc -l)
    before1_2=$((before1_1 / 4))
    before2_1=$(cat ${workdir}/${sourcefolder}/${samplename}_2.fastq | wc -l)
    before2_2=$((before2_1 / 4))
    after1_1=$(cat ${workdir}/LowComplexFiltered/${samplename}_1.fastq | wc -l)
    after1_2=$((after1_1 / 4))
    after2_1=$(cat ${workdir}/LowComplexFiltered/${samplename}_2.fastq | wc -l)
    after2_2=$((after2_1 / 4))
    difference1=$((before1_2 - after1_2))
    difference2=$((before2_2 - after2_2))
    percentage1=$((100-(after1_2 * 100 / before1_2 )))
    percentage2=$((100-(after2_2 * 100 / before2_2 )))
    #Print statistics
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now | 		From sample $samplename, $difference1 (PE1) and $difference2 (PE2) reads (${percentage1}% and ${percentage2}%) were removed due to low complexity" >> ${workdir}/run_${timestamp}.log
  	#Compress source files
		if [[ compress == TRUE ]]; then
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now | 		Compressing files ${sourcefolder}/${samplename}_1.fastq and ${sourcefolder}/${samplename}_2.fastq" >> ${workdir}/run_${timestamp}.log
    pigz -p ${threads} ${workdir}/${sourcefolder}/${samplename}_1.fastq
    pigz -p ${threads} ${workdir}/${sourcefolder}/${samplename}_2.fastq
    fi
  else
    #It is SR
    echo "$now | 		Removing low complexity reads from SR sample ${samplename}" >> ${workdir}/run_${timestamp}.log
    prinseq-lite.pl -lc_method "dust" -lc_threshold ${dustvalue} -fastq  ${workdir}/${sourcefolder}/${samplename}.fastq -out_good ${workdir}/LowComplexFiltered/${samplename} -out_bad null
    #Compute statistics
    before1=$(cat ${workdir}/${sourcefolder}/${samplename}.fastq | wc -l)
    before2=$((before1 / 4))
    after1=$(cat ${workdir}/LowComplexFiltered/${samplename}.fastq | wc -l)
    after2=$((after1 / 4))
    difference=$((before2 - after2))
    percentage=$((100-(after2 * 100 / before2 )))
    #Print statistics
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now | 		From sample $samplename, $difference reads (${percentage}%) were removed due to low complexity" >> ${workdir}/run_${timestamp}.log
    #Compress source file
    if [[ compress == TRUE ]]; then
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now | 		Compressing file ${sourcefolder}/${samplename}.fastq" >> ${workdir}/run_${timestamp}.log
    pigz -p ${threads} ${workdir}/${sourcefolder}/${samplename}.fastq
    fi
  fi
}

#Loop in parallel across samples specified in sample.data.txt
#Export function lowcompjob
export -f lowcompjob
parallel -j ${threads} -k lowcompjob {} ${settingsfile} ${sourcefolder} <${sampledatafile}

#Remove low complexity files
rm ${workdir}/LowComplexFiltered/*_prinseq_bad_*

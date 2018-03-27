#Source settings
source "$metafunkdirectory/settings.sh"

#Create QualityFiltered directory
mkdir -p ${workingdirectory}/${project}/QualityFiltered

#Loop across samples specified in sample.data.txt
while read sample; do

  #Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  samplefile=$(echo $sample | cut -d ' ' -f2 )
  now=$(date +"%Y-%d-%m %H:%M:%S")

  if [[ $samplefile =~ "/" ]]; then
    #It is PE
    echo "$now | 		Quality filtering sample $samplename" >> ${workingdirectory}/${project}/run.log
    #Perform quality filtering and rename output files
    AdapterRemoval --file1 ${workingdirectory}/${project}/RawData/${samplename}_1.fastq --file2 ${workingdirectory}/${project}/RawData/${samplename}_2.fastq --basename ${workingdirectory}/${project}/QualityFiltered/${samplename} --minquality ${minavgquality} --minlength ${minseqlength} --trimqualities --trimns --maxns 5 --qualitymax ${qualitymax} --threads ${threads}
    mv ${workingdirectory}/${project}/QualityFiltered/${samplename}.pair1.truncated ${workingdirectory}/${project}/QualityFiltered/${samplename}_1.fastq
    mv ${workingdirectory}/${project}/QualityFiltered/${samplename}.pair2.truncated ${workingdirectory}/${project}/QualityFiltered/${samplename}_2.fastq
    #Compute statistics
    before1_1=$(cat ${workingdirectory}/${project}/RawData/${samplename}_1.fastq | wc -l)
    before1_2=$((before1_1 / 4))
    before2_1=$(cat ${workingdirectory}/${project}/RawData/${samplename}_2.fastq | wc -l)
    before2_2=$((before2_1 / 4))
    after1_1=$(cat ${workingdirectory}/${project}/QualityFiltered/${samplename}_1.fastq | wc -l)
    after1_2=$((after1_1 / 4))
    after2_1=$(cat ${workingdirectory}/${project}/QualityFiltered/${samplename}_2.fastq | wc -l)
    after2_2=$((after2_1 / 4))
    difference1=$((before1_2 - after1_2))
    difference2=$((before2_2 - after2_2))
    percentage1=$((100-(after1_2 * 100 / before1_2 )))
    percentage2=$((100-(after2_2 * 100 / before2_2 )))
    #Print statistics
    now=$(date +"%Y-%d-%m %H:%M:%S")
    echo "$now | 		From sample $samplename, $difference1 (PE1) and $difference2 (PE2) reads (${percentage1}% and ${percentage2}%) were removed due to low quality" >> ${workingdirectory}/${project}/run.log
  else
    #It is SR
    echo "$now | 		Quality filtering sample $samplename" >> ${workingdirectory}/${project}/run.log
    #Perform quality filtering and rename output files
    AdapterRemoval --file1 ${workingdirectory}/${project}/RawData/${samplename}.fastq --basename ${workingdirectory}/${project}/QualityFiltered/${samplename} --minquality ${minavgquality} --minlength ${minseqlength} --trimqualities --trimns --maxns 5 --qualitymax ${qualitymax} --threads ${threads}
    mv ${workingdirectory}/${project}/QualityFiltered/${samplename}.truncated ${workingdirectory}/${project}/QualityFiltered/${samplename}.fastq
    #Compute statistics
    before1=$(cat ${workingdirectory}/${project}/RawData/${samplename}.fastq | wc -l)
    before2=$((before1 / 4))
    after1=$(cat ${workingdirectory}/${project}/DuplicatesRemoved/${samplename}.fastq | wc -l)
    after2=$((after1 / 4))
    difference=$((before2 - after2))
    percentage=$((100-(after2 * 100 / before2 )))
    #Print statistics
    now=$(date +"%Y-%d-%m %H:%M:%S")
    echo "$now | 		From sample $samplename, $difference reads (${percentage}%) were removed due to low quality" >> ${workingdirectory}/${project}/run.log
  fi

done < ${metafunkdirectory}/sample.data.txt

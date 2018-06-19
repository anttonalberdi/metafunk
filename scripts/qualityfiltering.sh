#Source settings
source $settingsfile

#Create QualityFiltered directory
mkdir -p ${workdir}/QualityFiltered

#Loop across samples specified in sample.data.txt
while read sample; do

  #Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  sampleinfo=$(echo $sample | cut -d ' ' -f2 )
  now=$(date +"%Y-%m-%d %H:%M:%S")

  if [[ $sampleinfo =~ "/" ]]; then
    #It is PE
    echo "$now |  Quality filtering sample $samplename" >> ${workdir}/run_${timestamp}.log
    #Perform quality filtering and rename output files
    AdapterRemoval --file1 ${workdir}/RawData/${samplename}_1.fastq --file2 ${workdir}/RawData/${samplename}_2.fastq --basename ${workdir}/QualityFiltered/${samplename} --minquality ${minavgquality} --minlength ${minseqlength} --trimqualities --trimns --maxns 5 --qualitymax ${qualitymax} --threads ${threads}
    mv ${workdir}/QualityFiltered/${samplename}.pair1.truncated ${workdir}/QualityFiltered/${samplename}_1.fastq
    mv ${workdir}/QualityFiltered/${samplename}.pair2.truncated ${workdir}/QualityFiltered/${samplename}_2.fastq
    #Compute statistics
    before1_1=$(cat ${workdir}/RawData/${samplename}_1.fastq | wc -l)
    before1_2=$((before1_1 / 4))
    before2_1=$(cat ${workdir}/RawData/${samplename}_2.fastq | wc -l)
    before2_2=$((before2_1 / 4))
    after1_1=$(cat ${workdir}/QualityFiltered/${samplename}_1.fastq | wc -l)
    after1_2=$((after1_1 / 4))
    after2_1=$(cat ${workdir}/QualityFiltered/${samplename}_2.fastq | wc -l)
    after2_2=$((after2_1 / 4))
    difference1=$((before1_2 - after1_2))
    difference2=$((before2_2 - after2_2))
    percentage1=$((100-(after1_2 * 100 / before1_2 )))
    percentage2=$((100-(after2_2 * 100 / before2_2 )))
    #Print statistics
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now | 		From sample $samplename, $difference1 (PE1) and $difference2 (PE2) reads (${percentage1}% and ${percentage2}%) were removed due to low quality" >> ${workdir}/run_${timestamp}.log
    #Compress source files
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now | 		Compressing files Rawdata/${samplename}_1.fastq and Rawdata/${samplename}_2.fastq" >> ${workdir}/run_${timestamp}.log
    pigz -p ${threads} ${workdir}/RawData/${samplename}_1.fastq
    pigz -p ${threads} ${workdir}/RawData/${samplename}_2.fastq
  else
    #It is SR
    echo "$now |  Quality filtering sample $samplename" >> ${workdir}/run_${timestamp}.log
    #Perform quality filtering and rename output files
    AdapterRemoval --file1 ${workdir}/RawData/${samplename}.fastq --basename ${workdir}/QualityFiltered/${samplename} --minquality ${minavgquality} --minlength ${minseqlength} --trimqualities --trimns --maxns 5 --qualitymax ${qualitymax} --threads ${threads}
    mv ${workdir}/QualityFiltered/${samplename}.truncated ${workdir}/QualityFiltered/${samplename}.fastq
    #Compute statistics
    before1=$(cat ${workdir}/RawData/${samplename}.fastq | wc -l)
    before2=$((before1 / 4))
    after1=$(cat ${workdir}/QualityFiltered/${samplename}.fastq | wc -l)
    after2=$((after1 / 4))
    difference=$((before2 - after2))
    percentage=$((100-(after2 * 100 / before2 )))
    #Print statistics
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now | 		From sample $samplename, $difference reads (${percentage}%) were removed due to low quality" >> ${workdir}/run_${timestamp}.log
    #Compress source file
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now | 		Compressing file Rawdata/${samplename}.fastq" >> ${workdir}/run_${timestamp}.log
    pigz -p ${threads} ${workdir}/RawData/${samplename}.fastq
  fi

done < ${sampledatafile}

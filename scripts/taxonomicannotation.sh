#Source settings file
source $settingsfile

#Create GeneMapping directory
mkdir -p ${workdir}/TaxonomicAnnotation

#Select source folder from which data will be retrieved
if [[ "$(ls -A ${workdir}/HumanDNARemoved/)" ]]; then
sourcefolder="HumanDNARemoved"
repair="no"
elif [[ "$(ls -A ${workdir}/HostDNARemoved/)" ]]; then
sourcefolder="HostDNARemoved"
repair="no"
elif [[ "$(ls -A ${workdir}/LowComplexFiltered/)" ]]; then
sourcefolder="LowComplexFiltered"
repair="yes"
elif [[ "$(ls -A ${workdir}/DuplicatesRemoved/)" ]]; then
sourcefolder="DuplicatesRemoved"
repair="yes"
elif [[ "$(ls -A ${workdir}/QualityFiltered/)" ]]; then
sourcefolder="QualityFiltered"
repair="yes"
else
sourcefolder="RawData"
repair="no"
fi

while read sample; do

		#Obtain data from sample.data.txt columns and get file name
		samplename=$(echo $sample | cut -d ' ' -f1)
		sampleinfo=$(echo $sample | cut -d ' ' -f2)
    
    if [[ $sampleinfo =~ "/" ]]; then
		#It is PE
      #Repair ends if necessary
      if [[ $repair == "yes" ]]; then
          now=$(date +"%Y-%m-%d %H:%M:%S")
          echo "$now | 		Repairing sample ${samplename}" >> ${workdir}/run_${timestamp}.log
          repair.sh in=${workdir}/${sourcefolder}/${samplename}_1.fastq in2=${workdir}/${sourcefolder}/${samplename}_2.fastq out=${workdir}/GeneMapping/${samplename}_1.fastq out2=${workdir}/GeneMapping/${samplename}_2.fastq
          sourcefolder="TaxonomicAnnotation"
      fi
 
    #Run metaphlan
    metaphlan2.py ${workdir}/${sourcefolder}/${query}_1.fastq,${workdir}/${sourcefolder}/${query}_2.fastq --input_type fastq --nproc ${threads}
		else
			#It is SR
    #Run metaphlan
    metaphlan2.py ${workdir}/${sourcefolder}/${query}.fastq --input_type fastq --nproc ${threads}

done < ${sampledatafile}


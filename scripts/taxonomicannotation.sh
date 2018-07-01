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
     	now=$(date +"%Y-%m-%d %H:%M:%S")
     	echo "$now | 		Assigning taxonomy to ${samplename}" >> ${workdir}/run_${timestamp}.log
    	metaphlan2.py ${workdir}/${sourcefolder}/${query}_1.fastq,${workdir}/${sourcefolder}/${query}_2.fastq --input_type fastq --nproc ${threads} --bowtie2out ${workdir}/TaxonomicAnnotation/${query}.bowtie2.bz2 -o ${workdir}/TaxonomicAnnotation/${query}.txt
		else
	#It is SR
    	#Run metaphlan
       	now=$(date +"%Y-%m-%d %H:%M:%S")
    	echo "$now | 		Assigning taxonomy to ${samplename}" >> ${workdir}/run_${timestamp}.log
	metaphlan2.py ${workdir}/${sourcefolder}/${query}.fastq --input_type fastq --nproc ${threads} -o ${workdir}/TaxonomicAnnotation/${query}.txt
	fi
done < ${sampledatafile}

#Merge all output files
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | 		Merging all taxonomic profiles" >> ${workdir}/run_${timestamp}.log
merge_metaphlan_tables.py ${workdir}/TaxonomicAnnotation/*.txt > ${workdir}/GeneTables/taxonomytable.txt

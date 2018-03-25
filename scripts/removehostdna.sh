#Source dependencies
source "$metafunkdirectory/settings.sh"

#Create LowComplexFiltered directory
mkdir -p ${workingdirectory}/${project}/HostDNARemoved
mkdir -p ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes

#Copy host reference genome to project directory
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Copying host genome(s)" >>  ${workingdirectory}/${project}/run.log

while read sample; do

	genomepath=$(echo $sample | cut -d ' ' -f4)
	genomefile=$(echo "${genomepath}"  | sed 's/.*\///')
	if [ ! -f ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile} ]; then
	cp ${genomepath}* ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes
	echo "$now |		Genome file $genomefile was copied to the project directory" >> ${workingdirectory}/${project}/run.log
	fi

done < ${metafunkdirectory}/sample.data.txt

#Index host reference genome
if [[ $indexedhostgenome == "yes" ]]; then

	while read sample; do

		genomepath=$(echo $sample | cut -d ' ' -f4)
		genomefile=$(echo "${genomepath}"  | sed 's/.*\///')
		now=$(date +"%Y-%d-%m %H:%M:%S")
		if [ ! -f ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile}.fai ]; then
		echo "$now | 		Indexing ${genomefile} genome" >> ${workingdirectory}/${project}/run.log
		samtools faidx ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile}
		bwa index ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile}
		now=$(date +"%Y-%d-%m %H:%M:%S")
		echo "$now | 		Genome ${genomefile} was succesfully indexed" >> ${workingdirectory}/${project}/run.log
		fi
	done < ${metafunkdirectory}/sample.data.txt

fi

#Select source folder from which data will be retrieved
if [[ "$(ls -A ${workingdirectory}/${project}/LowComplexFiltered/)" ]]; then
sourcefolder="LowComplexFiltered"
elif [[ "$(ls -A ${workingdirectory}/${project}/DuplicatesRemoved/)" ]]; then
sourcefolder="DuplicatesRemoved"
elif [[ "$(ls -A ${workingdirectory}/${project}/QualityFiltered/)" ]]; then
sourcefolder="QualityFiltered"
else
sourcefolder="RawData"
fi



now=$(date +"%Y-%d-%m %H:%M:%S")
#Remove host genome
if [[ $seqtype == "SR" ]]; then
echo "$now | 		Removing host DNA from SR data from directory ${sourcefolder}" >> ${workingdirectory}/${project}/run.log

		#Loop across samples specified in sample.data.txt
		while read sample; do

			#Obtain data from sample.data.txt columns and get file name
			samplename=$(echo $sample | cut -d ' ' -f1 )

			#Map reads against the reference genome and retrieve unmapped reads
			echo "				Removing host DNA from sample $samplename" >> ${workingdirectory}/${project}/run.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile} ${workingdirectory}/${project}/${sourcefolder}/${samplename}.fastq | samtools view -b -f4 - > ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam
			samtools fastq -0 ${workingdirectory}/${project}/HostDNARemoved/${samplename}.fastq ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam
			rm ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam

		done < ${metafunkdirectory}/sample.data.txt

elif [[ $seqtype == "PE" ]]; then
echo "$now | 		Removing host DNA from PE data from folder ${sourcefolder}" >> ${workingdirectory}/${project}/run.log

	#Loop across samples specified in sample.data.txt
	while read sample; do

		#Obtain data from sample.data.txt columns and get file name
		samplename=$(echo $sample | cut -d ' ' -f1 )
		genomepath=$(echo $sample | cut -d ' ' -f4)
		genomefile=$(echo "$genomepath"  | sed 's/.*\///')

		#Prevent repeating operation when looping through the 2nd pair
		if [[ ! -f ${workingdirectory}/${project}/HostDNARemoved/${samplename}_1.fastq ]]; then
		#Remove unpaired reads - Does not re-pair correctly!!! FIXX!!!
		now=$(date +"%Y-%d-%m %H:%M:%S")
		echo "$now | 			Repairing sample ${samplename}" >> ${workingdirectory}/${project}/run.log
		perl ${metafunkdirectory}/scripts/repairreads.pl -f1 ${workingdirectory}/${project}/${sourcefolder}/${samplename}_1.fastq -f2 ${workingdirectory}/${project}/${sourcefolder}/${samplename}_2.fastq -r '^@(\S+) [1|2](\S+)' -t -o ${workingdirectory}/${project}/HostDNARemoved/${samplename}
		#BBMap script below not repairing correctly
			#repair.sh in=${workingdirectory}/${project}/${sourcefolder}/${samplename}_1.fastq in2=${workingdirectory}/${project}/${sourcefolder}/${samplename}_2.fastq out=${workingdirectory}/${project}/HostDNARemoved/${samplename}_1.fastq out2=${workingdirectory}/${project}/HostDNARemoved/${samplename}_2.fastq
		#Map reads against the reference genome and retrieve unmapped reads
		now=$(date +"%Y-%d-%m %H:%M:%S")
		echo "$now | 			Removing host DNA from sample $samplename" >> ${workingdirectory}/${project}/run.log
		bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile} ${workingdirectory}/${project}/${sourcefolder}/${samplename}_1.fastq ${workingdirectory}/${project}/${sourcefolder}/${samplename}_2.fastq | samtools view -b -f12 - > ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam
		samtools fastq -1 ${workingdirectory}/${project}/HostDNARemoved/${samplename}_1.fastq -2 ${workingdirectory}/${project}/HostDNARemoved/${samplename}_2.fastq ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam
		rm ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam
		now=$(date +"%Y-%d-%m %H:%M:%S")
		echo "$now | 			Host DNA succesfully removed from sample $samplename" >> ${workingdirectory}/${project}/run.log
		fi
	done < ${metafunkdirectory}/sample.data.txt

else
echo "$now | 		ERROR: Sequencing read type has not been specified. It needs to be either SR or PE" >> ${workingdirectory}/${project}/run.log
exit
fi

#Source settings
source "$metafunkdirectory/settings.sh"

#Create LowComplexFiltered directory
mkdir -p ${workingdirectory}/${project}/HostDNARemoved
mkdir -p ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes

#Copy host reference genome to project directory
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Copying host genome(s)" >>  ${workingdirectory}/${project}/run.log

while read sample; do

	genomepath=$(echo $sample | cut -d ' ' -f3)
	genomefile=$(echo "${genomepath}"  | sed 's/.*\///')
	if [ ! -f ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile} ]; then
	echo "$now |		Copying genome file $genomefile to the project directory" >> ${workingdirectory}/${project}/run.log
	cp ${genomepath}* ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes
	now=$(date +"%Y-%d-%m %H:%M:%S")
	echo "$now |		Genome file $genomefile was copied to the project directory" >> ${workingdirectory}/${project}/run.log

done < ${sampledatafile}

#Index host reference genome
if [[ $indexhostgenome == "yes" ]]; then

	while read sample; do
		genomepath=$(echo $sample | cut -d ' ' -f3)
		genomefile=$(echo "${genomepath}"  | sed 's/.*\///')
		now=$(date +"%Y-%d-%m %H:%M:%S")
		if [ ! -f ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile}.fai ]; then
		echo "$now | 		Indexing ${genomefile} genome" >> ${workingdirectory}/${project}/run.log
		samtools faidx ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile}
		bwa index ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile}
		now=$(date +"%Y-%d-%m %H:%M:%S")
		echo "$now | 		Genome ${genomefile} was succesfully indexed" >> ${workingdirectory}/${project}/run.log
		fi
	done < ${sampledatafile}

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

#Map to host genome

now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Removing host DNA from files in directory ${sourcefolder}" >> ${workingdirectory}/${project}/run.log

while read sample; do

		#Obtain data from sample.data.txt columns and get file name
		samplename=$(echo $sample | cut -d ' ' -f1)
		samplefile=$(echo $sample | cut -d ' ' -f2)
		genomepath=$(echo $sample | cut -d ' ' -f3)
		genomefile=$(echo "$genomepath"  | sed 's/.*\///')

		if [[ $samplefile =~ "/" ]]; then
			#It is PE
			#Remove unpaired reads
			now=$(date +"%Y-%d-%m %H:%M:%S")
			echo "$now | 			Repairing sample ${samplename}" >> ${workingdirectory}/${project}/run.log
			repair.sh in=${workingdirectory}/${project}/${sourcefolder}/${samplename}_1.fastq in2=${workingdirectory}/${project}/${sourcefolder}/${samplename}_2.fastq out=${workingdirectory}/${project}/HostDNARemoved/${samplename}_1.fastq out2=${workingdirectory}/${project}/HostDNARemoved/${samplename}_2.fastq
			#Map reads against the reference genome and retrieve unmapped reads
			now=$(date +"%Y-%d-%m %H:%M:%S")
			echo "$now | 			Removing host DNA from sample $samplename" >> ${workingdirectory}/${project}/run.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile} ${workingdirectory}/${project}/HostDNARemoved/${samplename}_1.fastq ${workingdirectory}/${project}/HostDNARemoved/${samplename}_2.fastq | samtools view -b -f12 - > ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam
			now=$(date +"%Y-%d-%m %H:%M:%S")
			#Check if output file has been created; otherwise, print error message and kill the job
			if [[ ! -s ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam ]]; then
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workingdirectory}/${project}/run.log
				exit
			fi
			#Convert BAM file to FASTQ
			samtools fastq -1 ${workingdirectory}/${project}/HostDNARemoved/${samplename}_1.fastq -2 ${workingdirectory}/${project}/HostDNARemoved/${samplename}_2.fastq ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam
			rm ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam
			#Compute statistics
	    before1=$(cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}_1.fastq | wc -l)
	    before2=$((before1 / 4))
	    after1=$(cat ${workingdirectory}/${project}/DuplicatesRemoved/${samplename}_1.fastq | wc -l)
	    after2=$((after1 / 4))
	    difference=$((before2 - after2))
	    percentage=$((100-(after2 * 100 / before2 )))
			#Print statistics
	  	echo "$now | 		From sample $samplename, $difference PE reads (${percentage}%) were mapped to the host genome" >> ${workingdirectory}/${project}/run.log

		else

			#It is SR
			#Map reads against the reference genome and retrieve unmapped reads
			echo "				Removing host DNA from sample $samplename" >> ${workingdirectory}/${project}/run.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile} ${workingdirectory}/${project}/${sourcefolder}/${samplename}.fastq | samtools view -b -f4 - > ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam
			if [[ ! -s ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam ]]; then
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workingdirectory}/${project}/run.log
				exit
			fi
			#Convert BAM file to FASTQ
			samtools fastq -0 ${workingdirectory}/${project}/HostDNARemoved/${samplename}.fastq ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam
			rm ${workingdirectory}/${project}/HostDNARemoved/${samplename}.bam
			#Compute statistics
			before1=$(cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}.fastq | wc -l)
			before2=$((before1 / 4))
			after1=$(cat ${workingdirectory}/${project}/DuplicatesRemoved/${samplename}.fastq | wc -l)
			after2=$((after1 / 4))
			difference=$((before2 - after2))
			percentage=$((100-(after2 * 100 / before2 )))
			#Print statistics
			echo "$now | 		From sample $samplename, $difference reads (${percentage}%) were mapped to the host genome" >> ${workingdirectory}/${project}/run.log
		fi
done < ${sampledatafile}

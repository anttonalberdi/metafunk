#Source settings file
source $settingsfile

#Create LowComplexFiltered directory
mkdir -p ${workdir}/HostDNARemoved
mkdir -p ${workdir}/HostDNARemoved/ReferenceGenomes

#Copy host reference genome to project directory
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Copying host genome(s)" >>  ${workdir}/run.log

while read sample; do

	genomepath=$(echo $sample | cut -d ' ' -f3)
	genomefile=$(echo "${genomepath}"  | sed 's/.*\///')
	if [ ! -f ${workdir}/HostDNARemoved/ReferenceGenomes/${genomefile} ]; then
	echo "$now |		Copying genome file $genomefile to the project directory" >> ${workdir}/run.log
	cp ${genomepath}* ${workdir}/HostDNARemoved/ReferenceGenomes
	now=$(date +"%Y-%d-%m %H:%M:%S")
	echo "$now |		Genome file $genomefile was copied to the project directory" >> ${workdir}/run.log
	fi
done < ${sampledatafile}

#Index host reference genome
if [[ $indexhostgenome == "yes" ]]; then

	while read sample; do
		genomepath=$(echo $sample | cut -d ' ' -f3)
		genomefile=$(echo "${genomepath}"  | sed 's/.*\///')
		now=$(date +"%Y-%d-%m %H:%M:%S")
		if [ ! -f ${workdir}/HostDNARemoved/ReferenceGenomes/${genomefile}.fai ]; then
		echo "$now | 		Indexing ${genomefile} genome" >> ${workdir}/run.log
		samtools faidx ${workdir}/HostDNARemoved/ReferenceGenomes/${genomefile}
		bwa index ${workdir}/HostDNARemoved/ReferenceGenomes/${genomefile}
		now=$(date +"%Y-%d-%m %H:%M:%S")
		echo "$now | 		Genome ${genomefile} was succesfully indexed" >> ${workdir}/run.log
		fi
	done < ${sampledatafile}

fi

#Select source folder from which data will be retrieved
if [[ "$(ls -A ${workdir}/LowComplexFiltered/)" ]]; then
sourcefolder="LowComplexFiltered"
elif [[ "$(ls -A ${workdir}/DuplicatesRemoved/)" ]]; then
sourcefolder="DuplicatesRemoved"
elif [[ "$(ls -A ${workdir}/QualityFiltered/)" ]]; then
sourcefolder="QualityFiltered"
else
sourcefolder="RawData"
fi

#Map to host genome

now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Removing host DNA from files in directory ${sourcefolder}" >> ${workdir}/run.log

while read sample; do

		#Obtain data from sample.data.txt columns and get file name
		samplename=$(echo $sample | cut -d ' ' -f1)
		sampleinfo=$(echo $sample | cut -d ' ' -f2)
		genomepath=$(echo $sample | cut -d ' ' -f3)
		genomefile=$(echo "$genomepath"  | sed 's/.*\///')

		if [[ $sampleinfo =~ "/" ]]; then
			#It is PE
			#Remove unpaired reads
			now=$(date +"%Y-%d-%m %H:%M:%S")
			echo "$now | 			Repairing sample ${samplename}" >> ${workdir}/run.log
			repair.sh in=${workdir}/${sourcefolder}/${samplename}_1.fastq in2=${workdir}/${sourcefolder}/${samplename}_2.fastq out=${workdir}/HostDNARemoved/${samplename}_1.fastq out2=${workdir}/HostDNARemoved/${samplename}_2.fastq
			#Map reads against the reference genome and retrieve unmapped reads
			now=$(date +"%Y-%d-%m %H:%M:%S")
			echo "$now | 			Removing host DNA from sample $samplename" >> ${workdir}/run.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/HostDNARemoved/ReferenceGenomes/${genomefile} ${workdir}/HostDNARemoved/${samplename}_1.fastq ${workdir}/HostDNARemoved/${samplename}_2.fastq | samtools view -b -f12 - > ${workdir}/HostDNARemoved/${samplename}.bam
			now=$(date +"%Y-%d-%m %H:%M:%S")
			#Check if output file has been created; otherwise, print error message and kill the job
			if [[ ! -s ${workdir}/HostDNARemoved/${samplename}.bam ]]; then
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run.log
				exit
			fi
			#Convert BAM file to FASTQ
			samtools fastq -1 ${workdir}/HostDNARemoved/${samplename}_1.fastq -2 ${workdir}/HostDNARemoved/${samplename}_2.fastq ${workdir}/HostDNARemoved/${samplename}.bam
			rm ${workdir}/HostDNARemoved/${samplename}.bam
			#Compute statistics
	    before1=$(cat ${workdir}/${sourcefolder}/${samplename}_1.fastq | wc -l)
	    before2=$((before1 / 4))
	    after1=$(cat ${workdir}/DuplicatesRemoved/${samplename}_1.fastq | wc -l)
	    after2=$((after1 / 4))
	    difference=$((before2 - after2))
	    percentage=$((100-(after2 * 100 / before2 )))
			#Print statistics
	  	echo "$now | 		From sample $samplename, $difference PE reads (${percentage}%) were mapped to the host genome" >> ${workdir}/run.log

		else

			#It is SR
			#Map reads against the reference genome and retrieve unmapped reads
			echo "				Removing host DNA from sample $samplename" >> ${workdir}/run.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/HostDNARemoved/ReferenceGenomes/${genomefile} ${workdir}/${sourcefolder}/${samplename}.fastq | samtools view -b -f4 - > ${workdir}/HostDNARemoved/${samplename}.bam
			if [[ ! -s ${workdir}/HostDNARemoved/${samplename}.bam ]]; then
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run.log
				exit
			fi
			#Convert BAM file to FASTQ
			samtools fastq -0 ${workdir}/HostDNARemoved/${samplename}.fastq ${workdir}/HostDNARemoved/${samplename}.bam
			rm ${workdir}/HostDNARemoved/${samplename}.bam
			#Compute statistics
			before1=$(cat ${workdir}/${sourcefolder}/${samplename}.fastq | wc -l)
			before2=$((before1 / 4))
			after1=$(cat ${workdir}/DuplicatesRemoved/${samplename}.fastq | wc -l)
			after2=$((after1 / 4))
			difference=$((before2 - after2))
			percentage=$((100-(after2 * 100 / before2 )))
			#Print statistics
			echo "$now | 		From sample $samplename, $difference reads (${percentage}%) were mapped to the host genome" >> ${workdir}/run.log
		fi
done < ${sampledatafile}

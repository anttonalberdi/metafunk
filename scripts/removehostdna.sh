#Source settings file
source $settingsfile

#Create LowComplexFiltered directory
mkdir -p ${workdir}/HostDNARemoved
mkdir -p ${workdir}/HostDNARemoved/ReferenceGenomes

#Copy host reference genome to project directory
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | 		Copying host genome(s)" >>  ${workdir}/run_${timestamp}.log

while read sample; do

	genomepath=$(echo $sample | cut -d ' ' -f3)
	genomefile=$(echo "${genomepath}"  | sed 's/.*\///')
	if [ ! -f ${workdir}/HostDNARemoved/ReferenceGenomes/${genomefile} ]; then
	echo "$now |		Copying genome file $genomefile to the project directory" >> ${workdir}/run_${timestamp}.log
	cp ${genomepath}* ${workdir}/HostDNARemoved/ReferenceGenomes
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now |		Genome file $genomefile was copied to the project directory" >> ${workdir}/run_${timestamp}.log
	fi
done < ${sampledatafile}

#Index host reference genome
if [[ $indexhostgenome == "yes" ]]; then
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now | 		Indexing host genome(s)" >> ${workdir}/run_${timestamp}.log

	while read sample; do
		genomepath=$(echo $sample | cut -d ' ' -f3)
		genomefile=$(echo "${genomepath}"  | sed 's/.*\///')
		now=$(date +"%Y-%m-%d %H:%M:%S")
		if [ ! -f ${workdir}/HostDNARemoved/ReferenceGenomes/${genomefile}.fai ]; then
		echo "$now | 		Indexing ${genomefile} genome" >> ${workdir}/run_${timestamp}.log
		samtools faidx ${workdir}/HostDNARemoved/ReferenceGenomes/${genomefile}
		bwa index ${workdir}/HostDNARemoved/ReferenceGenomes/${genomefile}
		now=$(date +"%Y-%m-%d %H:%M:%S")
		echo "$now | 		Genome ${genomefile} was succesfully indexed" >> ${workdir}/run_${timestamp}.log
		fi
	done < ${sampledatafile}
else
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now | 		Host genome(s) will not be indexed" >> ${workdir}/run_${timestamp}.log

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

now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | 		Removing host DNA from files in directory ${sourcefolder}" >> ${workdir}/run_${timestamp}.log

while read sample; do

		#Obtain data from sample.data.txt columns and get file name
		samplename=$(echo $sample | cut -d ' ' -f1)
		sampleinfo=$(echo $sample | cut -d ' ' -f2)
		genomepath=$(echo $sample | cut -d ' ' -f3)
		genomefile=$(echo "$genomepath"  | sed 's/.*\///')

		if [[ $sampleinfo =~ "/" ]]; then
			#It is PE
			#Remove unpaired reads
			now=$(date +"%Y-%m-%d %H:%M:%S")
			echo "$now | 			Repairing sample ${samplename}" >> ${workdir}/run_${timestamp}.log
			repair.sh in=${workdir}/${sourcefolder}/${samplename}_1.fastq in2=${workdir}/${sourcefolder}/${samplename}_2.fastq out=${workdir}/HostDNARemoved/${samplename}_1.fastq out2=${workdir}/HostDNARemoved/${samplename}_2.fastq
			#Map reads against the reference genome and retrieve unmapped reads
			now=$(date +"%Y-%m-%d %H:%M:%S")
			echo "$now | 			Removing host DNA from sample $samplename" >> ${workdir}/run_${timestamp}.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/HostDNARemoved/ReferenceGenomes/${genomefile} ${workdir}/HostDNARemoved/${samplename}_1.fastq ${workdir}/HostDNARemoved/${samplename}_2.fastq | samtools view -b -f12 - > ${workdir}/HostDNARemoved/${samplename}.bam
			#Check if output file has been created; otherwise, print error message and kill the job
			if [[ ! -s ${workdir}/HostDNARemoved/${samplename}.bam ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			fi
			#Convert BAM file to FASTQ
			samtools fastq -1 ${workdir}/HostDNARemoved/${samplename}_1.fastq -2 ${workdir}/HostDNARemoved/${samplename}_2.fastq ${workdir}/HostDNARemoved/${samplename}.bam
			rm ${workdir}/HostDNARemoved/${samplename}.bam
			#Compute statistics
	    before1=$(cat ${workdir}/${sourcefolder}/${samplename}_1.fastq | wc -l)
	    before2=$((before1 / 4))
	    after1=$(cat ${workdir}/HostDNARemoved/${samplename}_1.fastq | wc -l)
	    after2=$((after1 / 4))
	    difference=$((before2 - after2))
	    percentage=$((100-(after2 * 100 / before2 )))
			#Print statistics
			now=$(date +"%Y-%m-%d %H:%M:%S")
	  	echo "$now | 			From sample $samplename, $difference PE reads (${percentage}%) were mapped to the host genome" >> ${workdir}/run_${timestamp}.log
			#Compress source files
		  #now=$(date +"%Y-%m-%d %H:%M:%S")
		  #echo "$now | 		Compressing files ${sourcefolder}/${samplename}_1.fastq and ${sourcefolder}/${samplename}_2.fastq" >> ${workdir}/run_${timestamp}.log
		  #pigz -p ${threads} ${workdir}/${sourcefolder}/${samplename}_1.fastq
		  #pigz -p ${threads} ${workdir}/${sourcefolder}/${samplename}_2.fastq
		else

			#It is SR
			#Map reads against the reference genome and retrieve unmapped reads
			now=$(date +"%Y-%m-%d %H:%M:%S")
			echo "				Removing host DNA from sample $samplename" >> ${workdir}/run_${timestamp}.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/HostDNARemoved/ReferenceGenomes/${genomefile} ${workdir}/${sourcefolder}/${samplename}.fastq | samtools view -b -f4 - > ${workdir}/HostDNARemoved/${samplename}.bam
			if [[ ! -s ${workdir}/HostDNARemoved/${samplename}.bam ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			fi
			#Convert BAM file to FASTQ
			samtools fastq -0 ${workdir}/HostDNARemoved/${samplename}.fastq ${workdir}/HostDNARemoved/${samplename}.bam
			rm ${workdir}/HostDNARemoved/${samplename}.bam
			#Compute statistics
			before1=$(cat ${workdir}/${sourcefolder}/${samplename}.fastq | wc -l)
			before2=$((before1 / 4))
			after1=$(cat ${workdir}/HostDNARemoved/${samplename}.fastq | wc -l)
			after2=$((after1 / 4))
			difference=$((before2 - after2))
			percentage=$((100-(after2 * 100 / before2 )))
			#Print statistics
			now=$(date +"%Y-%m-%d %H:%M:%S")
			echo "$now | 		From sample $samplename, $difference reads (${percentage}%) were mapped to the host genome" >> ${workdir}/run_${timestamp}.log
			#Compress source file
	    #now=$(date +"%Y-%m-%d %H:%M:%S")
	    #echo "$now | 		Compressing file ${sourcefolder}/${samplename}.fastq" >> ${workdir}/run_${timestamp}.log
	    #pigz -p ${threads} ${workdir}/${sourcefolder}/${samplename}.fastq
		fi
done < ${sampledatafile}

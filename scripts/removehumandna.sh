#Source settings file
source $settingsfile

#Create LowComplexFiltered directory
mkdir -p ${workdir}/HumanDNARemoved
mkdir -p ${workdir}/HumanDNARemoved/ReferenceGenome

#Copy host reference genome to project directory
humangenomefile=$(echo "${humangenomepath}" | sed 's/.*\///')
if [ ! -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile} ]; then
	now=$(date +"%Y-%d-%m %H:%M:%S")
	echo "$now | 		Copying human genome" >>  ${workdir}/run_${timestamp}.log
	cp ${humangenome}* ${workdir}/HumanDNARemoved/ReferenceGenome
	now=$(date +"%Y-%d-%m %H:%M:%S")
	echo "$now |		Human genome file $humangenomefile was copied to the project directory" >> ${workdir}/run_${timestamp}.log
else
	echo "$now | 		Human genome ${genomefile} already exists in the project directory" >> ${workdir}/run_${timestamp}.log
fi

#Index human reference genome
if [[ $indexhumangenome == "yes" ]]; then
	if [ ! -f ${workdir}/HumanDNARemoved/ReferenceGenome/${genomefile}.fai ]; then
		now=$(date +"%Y-%d-%m %H:%M:%S")
		echo "$now | 		Indexing human genome" >> ${workdir}/run_${timestamp}.log
		samtools faidx ${workdir}/HumanDNARemoved/ReferenceGenome/${genomefile}
		bwa index ${workdir}/HumanDNARemoved/ReferenceGenome/${genomefile}
		echo "$now | 		Human genome ${genomefile} was succesfully indexed" >> ${workdir}/run_${timestamp}.log
	else
		echo "$now | 		Human genome ${genomefile} is already indexed" >> ${workdir}/run_${timestamp}.log
	fi
fi

#Select source folder from which data will be retrieved
if [[ "$(ls -A ${workdir}/HostDNARemoved/)" ]]; then
sourcefolder="HostDNARemoved"
elif [[ "$(ls -A ${workdir}/LowComplexFiltered/)" ]]; then
sourcefolder="LowComplexFiltered"
elif [[ "$(ls -A ${workdir}/DuplicatesRemoved/)" ]]; then
sourcefolder="DuplicatesRemoved"
elif [[ "$(ls -A ${workdir}/QualityFiltered/)" ]]; then
sourcefolder="QualityFiltered"
else
sourcefolder="RawData"
fi

#Map to human genome

now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Removing human DNA from files in directory ${sourcefolder}" >> ${workdir}/run_${timestamp}.log

while read sample; do

		#Obtain data from sample.data.txt columns and get file name
		samplename=$(echo $sample | cut -d ' ' -f1)
		sampleinfo=$(echo $sample | cut -d ' ' -f2)
		humangenomefile=$(echo "${humangenomepath}" | sed 's/.*\///')

		if [[ $sampleinfo =~ "/" ]]; then
			#It is PE
			#Remove unpaired reads
			now=$(date +"%Y-%d-%m %H:%M:%S")
			echo "$now | 			Repairing sample ${samplename}" >> ${workdir}/run_${timestamp}.log
			repair.sh in=${workdir}/${sourcefolder}/${samplename}_1.fastq in2=${workdir}/${sourcefolder}/${samplename}_2.fastq out=${workdir}/HumanDNARemoved/${samplename}_1.fastq out2=${workdir}/HumanDNARemoved/${samplename}_2.fastq
			#Map reads against the reference genome and retrieve unmapped reads
			now=$(date +"%Y-%d-%m %H:%M:%S")
			echo "$now | 			Removing human DNA from sample $samplename" >> ${workdir}/run_${timestamp}.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/HumanDNARemoved/ReferenceGenome/${genomefile} ${workdir}/HumanDNARemoved/${samplename}_1.fastq ${workdir}/HumanDNARemoved/${samplename}_2.fastq | samtools view -b -f12 - > ${workdir}/HumanDNARemoved/${samplename}.bam
			now=$(date +"%Y-%d-%m %H:%M:%S")
			#Check if output file has been created; otherwise, print error message and kill the job
			if [[ ! -s ${workdir}/HumanDNARemoved/${samplename}.bam ]]; then
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			fi
			#Convert BAM file to FASTQ
			samtools fastq -1 ${workdir}/HumanDNARemoved/${samplename}_1.fastq -2 ${workdir}/HumanDNARemoved/${samplename}_2.fastq ${workdir}/HumanDNARemoved/${samplename}.bam
			rm ${workdir}/HumanDNARemoved/${samplename}.bam
			#Compute statistics
	    before1=$(cat ${workdir}/${sourcefolder}/${samplename}_1.fastq | wc -l)
	    before2=$((before1 / 4))
	    after1=$(cat ${workdir}/HumanDNARemoved/${samplename}_1.fastq | wc -l)
	    after2=$((after1 / 4))
	    difference=$((before2 - after2))
	    percentage=$((100-(after2 * 100 / before2 )))
			#Print statistics
	  	echo "$now | 		From sample $samplename, $difference PE reads (${percentage}%) were mapped to the human genome" >> ${workdir}/run_${timestamp}.log

		else

			#It is SR
			#Map reads against the reference genome and retrieve unmapped reads
			echo "				Removing human DNA from sample $samplename" >> ${workdir}/run_${timestamp}.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/HumanDNARemoved/ReferenceGenome/${genomefile} ${workdir}/${sourcefolder}/${samplename}.fastq | samtools view -b -f4 - > ${workdir}/HumanDNARemoved/${samplename}.bam
			if [[ ! -s ${workdir}/HumanDNARemoved/${samplename}.bam ]]; then
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			fi
			#Convert BAM file to FASTQ
			samtools fastq -0 ${workdir}/HumanDNARemoved/${samplename}.fastq ${workdir}/HumanDNARemoved/${samplename}.bam
			rm ${workdir}/HumanDNARemoved/${samplename}.bam
			#Compute statistics
			before1=$(cat ${workdir}/${sourcefolder}/${samplename}.fastq | wc -l)
			before2=$((before1 / 4))
			after1=$(cat ${workdir}/HumanDNARemoved/${samplename}.fastq | wc -l)
			after2=$((after1 / 4))
			difference=$((before2 - after2))
			percentage=$((100-(after2 * 100 / before2 )))
			#Print statistics
			echo "$now | 		From sample $samplename, $difference reads (${percentage}%) were mapped to the human genome" >> ${workdir}/run_${timestamp}.log
		fi
done < ${sampledatafile}

#Source settings
source "$metafunkdirectory/settings.sh"

#Create LowComplexFiltered directory
mkdir -p ${workingdirectory}/${project}/HumanDNARemoved
mkdir -p ${workingdirectory}/${project}/HumanDNARemoved/ReferenceGenome

#Copy host reference genome to project directory
humangenomefile=$(echo "${humangenomepath}" | sed 's/.*\///')
if [ ! -f ${workingdirectory}/${project}/HumanDNARemoved/ReferenceGenome/${humangenomefile} ]; then
	now=$(date +"%Y-%d-%m %H:%M:%S")
	echo "$now | 		Copying human genome" >>  ${workingdirectory}/${project}/run.log
	cp ${humangenome}* ${workingdirectory}/${project}/HumanDNARemoved/ReferenceGenome
	now=$(date +"%Y-%d-%m %H:%M:%S")
	echo "$now |		Human genome file $humangenomefile was copied to the project directory" >> ${workingdirectory}/${project}/run.log
else
	echo "$now | 		Human genome ${genomefile} already exists in the project directory" >> ${workingdirectory}/${project}/run.log
fi

#Index human reference genome
if [[ $indexhumangenome == "yes" ]]; then
	if [ ! -f ${workingdirectory}/${project}/HumanDNARemoved/ReferenceGenome/${genomefile}.fai ]; then
		now=$(date +"%Y-%d-%m %H:%M:%S")
		echo "$now | 		Indexing human genome" >> ${workingdirectory}/${project}/run.log
		samtools faidx ${workingdirectory}/${project}/HumanDNARemoved/ReferenceGenome/${genomefile}
		bwa index ${workingdirectory}/${project}/HumanDNARemoved/ReferenceGenome/${genomefile}
		echo "$now | 		Human genome ${genomefile} was succesfully indexed" >> ${workingdirectory}/${project}/run.log
	else
		echo "$now | 		Human genome ${genomefile} is already indexed" >> ${workingdirectory}/${project}/run.log
	fi
fi

#Select source folder from which data will be retrieved
if [[ "$(ls -A ${workingdirectory}/${project}/HostDNARemoved/)" ]]; then
sourcefolder="HostDNARemoved"
elif [[ "$(ls -A ${workingdirectory}/${project}/LowComplexFiltered/)" ]]; then
sourcefolder="LowComplexFiltered"
elif [[ "$(ls -A ${workingdirectory}/${project}/DuplicatesRemoved/)" ]]; then
sourcefolder="DuplicatesRemoved"
elif [[ "$(ls -A ${workingdirectory}/${project}/QualityFiltered/)" ]]; then
sourcefolder="QualityFiltered"
else
sourcefolder="RawData"
fi

#Map to human genome

now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Removing human DNA from files in directory ${sourcefolder}" >> ${workingdirectory}/${project}/run.log

while read sample; do

		#Obtain data from sample.data.txt columns and get file name
		samplename=$(echo $sample | cut -d ' ' -f1)
		samplefile=$(echo $sample | cut -d ' ' -f2)
		humangenomefile=$(echo "${humangenomepath}" | sed 's/.*\///')

		if [[ $samplefile =~ "/" ]]; then
			#It is PE
			#Remove unpaired reads
			now=$(date +"%Y-%d-%m %H:%M:%S")
			echo "$now | 			Repairing sample ${samplename}" >> ${workingdirectory}/${project}/run.log
			repair.sh in=${workingdirectory}/${project}/${sourcefolder}/${samplename}_1.fastq in2=${workingdirectory}/${project}/${sourcefolder}/${samplename}_2.fastq out=${workingdirectory}/${project}/HumanDNARemoved/${samplename}_1.fastq out2=${workingdirectory}/${project}/HumanDNARemoved/${samplename}_2.fastq
			#Map reads against the reference genome and retrieve unmapped reads
			now=$(date +"%Y-%d-%m %H:%M:%S")
			echo "$now | 			Removing human DNA from sample $samplename" >> ${workingdirectory}/${project}/run.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workingdirectory}/${project}/HumanDNARemoved/ReferenceGenome/${genomefile} ${workingdirectory}/${project}/HumanDNARemoved/${samplename}_1.fastq ${workingdirectory}/${project}/HumanDNARemoved/${samplename}_2.fastq | samtools view -b -f12 - > ${workingdirectory}/${project}/HumanDNARemoved/${samplename}.bam
			now=$(date +"%Y-%d-%m %H:%M:%S")
			#Check if output file has been created; otherwise, print error message and kill the job
			if [[ ! -s ${workingdirectory}/${project}/HumanDNARemoved/${samplename}.bam ]]; then
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workingdirectory}/${project}/run.log
				exit
			fi
			#Convert BAM file to FASTQ
			samtools fastq -1 ${workingdirectory}/${project}/HumanDNARemoved/${samplename}_1.fastq -2 ${workingdirectory}/${project}/HumanDNARemoved/${samplename}_2.fastq ${workingdirectory}/${project}/HumanDNARemoved/${samplename}.bam
			rm ${workingdirectory}/${project}/HumanDNARemoved/${samplename}.bam
			#Compute statistics
	    before1=$(cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}_1.fastq | wc -l)
	    before2=$((before1 / 4))
	    after1=$(cat ${workingdirectory}/${project}/HumanDNARemoved/${samplename}_1.fastq | wc -l)
	    after2=$((after1 / 4))
	    difference=$((before2 - after2))
	    percentage=$((100-(after2 * 100 / before2 )))
			#Print statistics
	  	echo "$now | 		From sample $samplename, $difference PE reads (${percentage}%) were mapped to the human genome" >> ${workingdirectory}/${project}/run.log

		else

			#It is SR
			#Map reads against the reference genome and retrieve unmapped reads
			echo "				Removing human DNA from sample $samplename" >> ${workingdirectory}/${project}/run.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workingdirectory}/${project}/HumanDNARemoved/ReferenceGenome/${genomefile} ${workingdirectory}/${project}/${sourcefolder}/${samplename}.fastq | samtools view -b -f4 - > ${workingdirectory}/${project}/HumanDNARemoved/${samplename}.bam
			if [[ ! -s ${workingdirectory}/${project}/HumanDNARemoved/${samplename}.bam ]]; then
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workingdirectory}/${project}/run.log
				exit
			fi
			#Convert BAM file to FASTQ
			samtools fastq -0 ${workingdirectory}/${project}/HumanDNARemoved/${samplename}.fastq ${workingdirectory}/${project}/HumanDNARemoved/${samplename}.bam
			rm ${workingdirectory}/${project}/HumanDNARemoved/${samplename}.bam
			#Compute statistics
			before1=$(cat ${workingdirectory}/${project}/${sourcefolder}/${samplename}.fastq | wc -l)
			before2=$((before1 / 4))
			after1=$(cat ${workingdirectory}/${project}/HumanDNARemoved/${samplename}.fastq | wc -l)
			after2=$((after1 / 4))
			difference=$((before2 - after2))
			percentage=$((100-(after2 * 100 / before2 )))
			#Print statistics
			echo "$now | 		From sample $samplename, $difference reads (${percentage}%) were mapped to the human genome" >> ${workingdirectory}/${project}/run.log
		fi
done < ${sampledatafile}

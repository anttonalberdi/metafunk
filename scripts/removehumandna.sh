#Source settings file
source $settingsfile

#Create LowComplexFiltered directory
mkdir -p ${workdir}/HumanDNARemoved
mkdir -p ${workdir}/HumanDNARemoved/ReferenceGenome

#Copy host reference genome to project directory
if [[ $humangenomepath == *.fasta.gz || $humangenomepath == *.fa.gz || $humangenomepath == *.fna.gz ]]; then
	humangenomefile=$(echo "${humangenomepath}"  | sed 's/.*\///' | sed 's/\.[^.]*$//')
elif [[ $humangenomepath == *.fasta || $humangenomepath == *.fa || $humangenomepath == *.fna ]]; then
	humangenomefile=$(echo "${humangenomepath}"  | sed 's/.*\///')
else
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now | 		ERROR! Genome ${humangenomefile} has an unsupported extension" >> ${workdir}/run_${timestamp}.log
fi

if [ ! -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile} ] && [ ! -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.gz ]; then
	if [[ $humangenomepath == *.fasta.gz || $humangenomepath == *.fa.gz || $humangenomepath == *.fna.gz ]]; then
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now | 		Copying human genome to project directory" >>  ${workdir}/run_${timestamp}.log
	cp ${humangenomepath}* ${workdir}/HumanDNARemoved/ReferenceGenome
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now |		Decompressing human genome" >> ${workdir}/run_${timestamp}.log
	gunzip ${workdir}/HumanDNARemoved/ReferenceGenomes/${humangenomefile}*.gz
	gunzip
	else
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now | 		Copying human genome to project directory" >>  ${workdir}/run_${timestamp}.log
	cp ${humangenomepath}* ${workdir}/HumanDNARemoved/ReferenceGenome
	fi
else
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now | 		Human genome ${genomefile} already exists in the project directory" >> ${workdir}/run_${timestamp}.log
fi

#Index human reference genome
	if [ ! -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.amb ] || [ ! -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.ann ] || [ ! -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.bwt ] || [ ! -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.fai ] || [ ! -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.pac ] || [ ! -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.sa ]; then
		now=$(date +"%Y-%m-%d %H:%M:%S")
		echo "$now | 		Indexing human genome" >> ${workdir}/run_${timestamp}.log
		samtools faidx ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}
		bwa index ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}

		now=$(date +"%Y-%m-%d %H:%M:%S")
		if [ -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.amb ] && [ -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.ann ] && [ -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.bwt ] && [ -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.fai ] && [ -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.pac ] && [ -f ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile}.sa ]; then
		echo "$now | 		Human genome ${humangenomefile} was succesfully indexed" >> ${workdir}/run_${timestamp}.log
		else
		echo "$now | 		There was an error when indexing the human genome ${humangenomefile}" >> ${workdir}/run_${timestamp}.log
		exit
		fi
	else
		echo "$now | 		Human genome ${humangenomefile} is already indexed" >> ${workdir}/run_${timestamp}.log
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

now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | 		Removing human DNA from files in directory ${sourcefolder}" >> ${workdir}/run_${timestamp}.log

while read sample; do

		#Obtain data from sample.data.txt columns and get file name
		samplename=$(echo $sample | cut -d ' ' -f1)
		sampleinfo=$(echo $sample | cut -d ' ' -f2)
		humangenomefile=$(echo "${humangenomepath}" | sed 's/.*\///')

		if [[ $sampleinfo =~ "/" ]]; then
			#It is PE
			#Decompress if needed
			if [[ -f ${workdir}/${sourcefolder}/${samplename}_1.fastq.gz ]]; then
			gunzip ${workdir}/${sourcefolder}/${samplename}_1.fastq.gz
			fi
			if [[ -f ${workdir}/${sourcefolder}/${samplename}_2.fastq.gz ]]; then
			gunzip ${workdir}/${sourcefolder}/${samplename}_2.fastq.gz
			fi
			#Repair paired-end reads using BBMap script repair.sh
			now=$(date +"%Y-%m-%d %H:%M:%S")
			echo "$now | 			Repairing sample ${samplename}" >> ${workdir}/run_${timestamp}.log
			repair.sh in=${workdir}/${sourcefolder}/${samplename}_1.fastq in2=${workdir}/${sourcefolder}/${samplename}_2.fastq out=${workdir}/HumanDNARemoved/${samplename}_1.fastq out2=${workdir}/HumanDNARemoved/${samplename}_2.fastq
			#Map reads against the reference genome and retrieve unmapped reads
			now=$(date +"%Y-%m-%d %H:%M:%S")
			echo "$now | 			Removing human DNA from sample $samplename" >> ${workdir}/run_${timestamp}.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile} ${workdir}/HumanDNARemoved/${samplename}_1.fastq ${workdir}/HumanDNARemoved/${samplename}_2.fastq | samtools view | grep -v -P '^@|NM:i:[0-2]\b' | samtools view -T ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile} -b - > ${workdir}/HumanDNARemoved/${samplename}.bam
			#OLD: bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile} ${workdir}/HumanDNARemoved/${samplename}_1.fastq ${workdir}/HumanDNARemoved/${samplename}_2.fastq | samtools view -b -f12 - > ${workdir}/HumanDNARemoved/${samplename}.bam
			#Check if output file has been created; otherwise, print error message and kill the job
			if [[ ! -s ${workdir}/HumanDNARemoved/${samplename}.bam ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			fi
			#Convert BAM file to FASTQ
			samtools fastq -s ${workdir}/HumanDNARemoved/${samplename}_singleton.fastq -1 ${workdir}/HumanDNARemoved/${samplename}_1.fastq -2 ${workdir}/HumanDNARemoved/${samplename}_2.fastq ${workdir}/HumanDNARemoved/${samplename}.bam
			rm ${workdir}/HumanDNARemoved/${samplename}.bam
			#Compute statistics
	    before1=$(cat ${workdir}/${sourcefolder}/${samplename}_1.fastq | wc -l)
	    before2=$((before1 / 4))
	    after1=$(cat ${workdir}/HumanDNARemoved/${samplename}_1.fastq | wc -l)
	    after2=$((after1 / 4))
	    difference=$((before2 - after2))
	    percentage=$((100-(after2 * 100 / before2 )))
			#Print statistics
			now=$(date +"%Y-%m-%d %H:%M:%S")
	  	echo "$now | 		From sample $samplename, $difference PE reads (${percentage}%) were mapped to the human genome" >> ${workdir}/run_${timestamp}.log
	#Compress source files
		if [[ compress == TRUE ]]; then
		  now=$(date +"%Y-%m-%d %H:%M:%S")
		  echo "$now | 		Compressing files ${sourcefolder}/${samplename}_1.fastq and ${sourcefolder}/${samplename}_2.fastq" >> ${workdir}/run_${timestamp}.log
		  pigz -p ${threads} ${workdir}/${sourcefolder}/${samplename}_1.fastq
		  pigz -p ${threads} ${workdir}/${sourcefolder}/${samplename}_2.fastq
		  fi
		else

		#It is SR
			#Decompress if needed
			if [[ -f ${workdir}/${sourcefolder}/${samplename}.fastq.gz ]]; then
			gunzip ${workdir}/${sourcefolder}/${samplename}.fastq.gz
			fi

			#Map reads against the reference genome and retrieve unmapped reads
			now=$(date +"%Y-%m-%d %H:%M:%S")
			echo "				Removing human DNA from sample $samplename" >> ${workdir}/run_${timestamp}.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile} ${workdir}/HumanDNARemoved/${samplename}.fastq | samtools view | grep -v -P '^@|NM:i:[0-2]\b' | samtools view -T ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile} -b - > ${workdir}/HumanDNARemoved/${samplename}.bam
			#OLD: bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/HumanDNARemoved/ReferenceGenome/${humangenomefile} ${workdir}/${sourcefolder}/${samplename}.fastq | samtools view -b -f4 - > ${workdir}/HumanDNARemoved/${samplename}.bam
			if [[ ! -s ${workdir}/HumanDNARemoved/${samplename}.bam ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
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
			now=$(date +"%Y-%m-%d %H:%M:%S")
			echo "$now | 		From sample $samplename, $difference reads (${percentage}%) were mapped to the human genome" >> ${workdir}/run_${timestamp}.log

	if [[ compress == TRUE ]]; then
	    now=$(date +"%Y-%m-%d %H:%M:%S")
	    echo "$now | 		Compressing file ${sourcefolder}/${samplename}.fastq" >> ${workdir}/run_${timestamp}.log
	    pigz -p ${threads} ${workdir}/${sourcefolder}/${samplename}.fastq
	fi
fi
done < ${sampledatafile}

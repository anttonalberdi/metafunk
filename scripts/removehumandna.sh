#Source dependencies
source "$metafunkdirectory/settings.sh"

#Create LowComplexFiltered directory
mkdir -p ${workingdirectory}/${project}/RemoveHumanDNA
mkdir -p ${workingdirectory}/${project}/RemoveHumanDNA/ReferenceGenome

#Copy host reference genome to project directory
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Copying human genome" >>  ${workingdirectory}/${project}/run.log
cp ${humangenome}* ${workingdirectory}/${project}/RemoveHumanDNA/ReferenceGenome
genomefile=$(echo "${humangenome}"  | sed 's/.*\///')
echo "				Genome file: $genomefile" >> ${workingdirectory}/${project}/run.log

#Index human reference genome
if [[ $indexedhhumangenome == "yes" ]]; then
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Indexing human genome" >> ${workingdirectory}/${project}/run.log
samtools faidx ${workingdirectory}/${project}/RemoveHumanDNA/ReferenceGenome/${genomefile}
bwa index ${workingdirectory}/${project}/RemoveHumanDNA/ReferenceGenome/${genomefile}
fi

#Remove human genome
if [[ $seqtype == "SR" ]]; then
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Removing human DNA from SR data" >> ${workingdirectory}/${project}/run.log
	#Iterate through samples
	while read samplefile; do
		#Select HostRemoved, LowComplex or Original data file
		if [[ $removehostdna == "yes" ]]; then
		inputfile=${workingdirectory}/${project}/RemoveHostDNA/${samplefile}.fastq
		elif [[ $removelowcomplexity == "yes" ]]; then
		inputfile=${workingdirectory}/${project}/LowComplexFiltered/${samplefile}.fastq
		else
		inputfile=${workingdirectory}/${project}/RawData/${samplefile}.fastq
		fi

		echo "				Input file: $inputfile" >> ${workingdirectory}/${project}/run.log

		#Map reads against the reference genome and retrieve unmapped reads
		bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workingdirectory}/${project}/RemoveHumanDNA/ReferenceGenome/${genomefile} ${inputfile} | samtools view -b -f4 - > ${workingdirectory}/${project}/RemoveHumanDNA/${samplefile}.bam
		samtools fastq -0 ${workingdirectory}/${project}/RemoveHumanDNA/${samplefile}.fastq ${workingdirectory}/${project}/RemoveHumanDNA/${samplefile}.bam
		rm ${workingdirectory}/${project}/RemoveHumanDNA/${samplefile}.bam
	done < ${metafunkdirectory}/sample.data.txt

elif [[ $seqtype == "PE" ]]; then
echo "Removing human DNA from PE data" >> ${workingdirectory}/${project}/run.log

else
echo "Sequencing read type has not been specified. It needs to be either SR or PE" >> ${workingdirectory}/${project}/run.log
fi

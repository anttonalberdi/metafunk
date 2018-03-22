#Source dependencies
source "$metafunkdirectory/settings.sh"

#Create LowComplexFiltered directory
mkdir -p ${workingdirectory}/${project}/RemoveHostDNA
mkdir -p ${workingdirectory}/${project}/RemoveHostDNA/referencegenome

#Copy host reference genome to project directory
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Copying host genome" >>  ${workingdirectory}/${project}/run.log
cp ${hostgenome}* ${workingdirectory}/${project}/RemoveHostDNA/referencegenome
genomefile=$(echo "${hostgenome}"  | sed 's/.*\///')
echo "			Genome file: $genomefile" >> ${workingdirectory}/${project}/run.log

#Index host reference genome
if [[ $indexedhostgenome == "yes" ]]; then  
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Indexing host genome" >> ${workingdirectory}/${project}/run.log
samtools faidx ${genomefile}
bwa index ${genomefile}
fi

#Remove host genome
if [[ $seqtype == "SR" ]]; then 
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Removing host DNA from SR data" >> ${workingdirectory}/${project}/run.log
	#Iterate through samples
	while read samplefile; do 
		#Select LowComplex or Original data file
		if [[ $removelowcomplexity == "yes" ]]; then  
		inputfile=${workingdirectory}/${project}/LowComplexFiltered/${samplefile}.fastq
		else
		inputfile=${workingdirectory}/${project}/RawData/${samplefile}.fastq
		fi

		echo "			Input file: $inputfile" >> ${workingdirectory}/${project}/run.log

		#Map reads against the reference genome and retrieve unmapped reads
		bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workingdirectory}/${project}/RemoveHostDNA/referencegenome/${genomefile} ${inputfile} | samtools view -b -f4 - > ${workingdirectory}/${project}/RemoveHostDNA/${samplefile}.bam
		samtools fastq -0 ${workingdirectory}/${project}/RemoveHostDNA/${samplefile}.fastq ${workingdirectory}/${project}/RemoveHostDNA/${samplefile}.bam
		rm ${workingdirectory}/${project}/RemoveHostDNA/${samplefile}.bam
	done < ${metafunkdirectory}/sample.data.txt

elif [[ $seqtype == "PE" ]]; then 
echo "Removing host DNA from PE data" >> ${workingdirectory}/${project}/run.log

else
echo "Sequencing read type has not been specified. It needs to be either SR or PE" >> run.log
fi



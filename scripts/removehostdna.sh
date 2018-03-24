#Source dependencies
source "$metafunkdirectory/settings.sh"

#Create LowComplexFiltered directory
mkdir -p ${workingdirectory}/${project}/HostDNARemoved
mkdir -p ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes

#Copy host reference genome to project directory
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | 		Copying host genome" >>  ${workingdirectory}/${project}/run.log

while read sample; do

	genomepath=$(echo $sample | cut -d ' ' -f4)
	genomefile=$(echo "${genomepath}"  | sed 's/.*\///')
	if [ ! -f ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile} ]; then
	cp ${genomepath}* ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes
	echo "				Genome file $genomefile was copied to the project directory" >> ${workingdirectory}/${project}/run.log
	fi

done < ${metafunkdirectory}/sample.data.txt

#Index host reference genome
if [[ $indexedhostgenome == "yes" ]]; then

	while read sample; do

		genomepath=$(echo $sample | cut -d ' ' -f4)
		genomefile=$(echo "${genomepath}"  | sed 's/.*\///')
		now=$(date +"%Y-%d-%m %H:%M:%S")
		echo "$now | 		Indexing ${genomefile} genome" >> ${workingdirectory}/${project}/run.log
		samtools faidx ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile}
		bwa index ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile}
		now=$(date +"%Y-%d-%m %H:%M:%S")
		echo "$now | 		Genome ${genomefile} was succesfully indexed" >> ${workingdirectory}/${project}/run.log

	done < ${metafunkdirectory}/sample.data.txt

fi

#Select source folder from which data will be retrieved
if [[ $removelowcomplexity == "yes" ]]; then
sourcefolder="LowComplexFiltered"
elif [[ $removeduplicates == "yes" ]]; then
sourcefolder="DuplicatesRemoved"
elif [[ $qualityfiltering == "yes" ]]; then
sourcefolder="QualityFiltered"
else
sourcefolder="RawData"
fi

now=$(date +"%Y-%d-%m %H:%M:%S")
#Remove host genome
if [[ $seqtype == "SR" ]]; then
echo "$now | 		Removing host DNA from SR data from folder ${sourcefolder}" >> ${workingdirectory}/${project}/run.log

		#Loop across samples specified in sample.data.txt
		while read sample; do

		  #Obtain data from sample.data.txt columns and get file name
		  samplename=$(echo $sample | cut -d ' ' -f1 )
		  sampleread=$(echo $sample | cut -d ' ' -f2 )
		  samplefile=${samplename}

			#Map reads against the reference genome and retrieve unmapped reads
			echo "				Removing host DNA from sample $sample" >> ${workingdirectory}/${project}/run.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile} ${workingdirectory}/${project}/${sourcefolder}/${samplefile}.fastq | samtools view -b -f4 - > ${workingdirectory}/${project}/HostDNARemoved/${samplefile}.bam
			samtools fastq -0 ${workingdirectory}/${project}/HostDNARemoved/${samplefile}.fastq ${workingdirectory}/${project}/HostDNARemoved/${samplefile}.bam
			rm ${workingdirectory}/${project}/HostDNARemoved/${samplefile}.bam

		done < ${metafunkdirectory}/sample.data.txt

elif [[ $seqtype == "PE" ]]; then
echo "$now | 		Removing host DNA from PE data from folder ${sourcefolder}" >> ${workingdirectory}/${project}/run.log

	#Loop across samples specified in sample.data.txt
	while read sample; do

		#Obtain data from sample.data.txt columns and get file name
		samplename=$(echo $sample | cut -d ' ' -f1 )
		sampleread=$(echo $sample | cut -d ' ' -f2 )
		samplefile=$(echo ${samplename}_${sampleread})

		#Map reads against the reference genome and retrieve unmapped reads
		echo "				Removing host DNA from sample $sample" >> ${workingdirectory}/${project}/run.log
		bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workingdirectory}/${project}/HostDNARemoved/ReferenceGenomes/${genomefile} ${workingdirectory}/${project}/${sourcefolder}/${samplefile}_1.fastq ${workingdirectory}/${project}/${sourcefolder}/${samplefile}_2.fastq | samtools view -b -f12 - > ${workingdirectory}/${project}/HostDNARemoved/${samplefile}.bam
		samtools fastq -1 ${workingdirectory}/${project}/HostDNARemoved/${samplefile}_1.fastq -2 ${workingdirectory}/${project}/HostDNARemoved/${samplefile}_1.fastq ${workingdirectory}/${project}/HostDNARemoved/${samplefile}.bam
		rm ${workingdirectory}/${project}/HostDNARemoved/${samplefile}.bam

else
echo "$now | 		ERROR: Sequencing read type has not been specified. It needs to be either SR or PE" >> ${workingdirectory}/${project}/run.log
exit
fi

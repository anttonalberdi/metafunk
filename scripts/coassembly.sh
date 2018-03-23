#Source dependencies
source "$metafunkdirectory/settings.sh"

#Create Coassembly directory
mkdir -p ${workingdirectory}/${project}/CoAssembly

#Convert to fasta  
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Converting FASTQ files to FASTA" >> ${workingdirectory}/${project}/run.log
	#Iterate through samples
	while read samplefile; do 
		#Select LowComplex or Original data file
		now=$(date +"%Y-%d-%m %H:%M:%S")
		if [[ $removehostdna == "yes" ]]; then 
		echo "$now | 		Converting $samplefile from host removed files" >> ${workingdirectory}/${project}/run.log
		inputfile=${workingdirectory}/${project}/RemoveHostDNA/${samplefile}.fastq
		elif [[ $removelowcomplexity == "yes" ]]; then  
		echo "$now | 		Converting $samplefile from low complexity removed files" >> ${workingdirectory}/${project}/run.log
		inputfile=${workingdirectory}/${project}/LowComplexFiltered/${samplefile}.fastq
		else
		echo "$now | 		Converting $samplefile from original raw files" >> ${workingdirectory}/${project}/run.log
		inputfile=${workingdirectory}/${project}/RawData/${samplefile}.fastq
		fi
	fastq_to_fasta -i ${inputfile} -o ${workingdirectory}/${project}/CoAssembly/${samplefile}.fasta
	done < ${metafunkdirectory}/sample.data.txt

#Concatenate
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Concatenating all samples" >> ${workingdirectory}/${project}/run.log
cat ${workingdirectory}/${project}/CoAssembly/*.fasta > ${workingdirectory}/${project}/CoAssembly/Allsamples.fasta

#Remove co-assembly directory if exists
if [[ $overridecoassembly == "yes" ]]; then 
rm -r ${workingdirectory}/${project}/CoAssembly/Megahit
fi
    
#Co-assembly
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Co-assembling all reads" >> ${workingdirectory}/${project}/run.log
megahit -t ${threads} -r ${workingdirectory}/${project}/CoAssembly/Allsamples.fasta -o ${workingdirectory}/${project}/CoAssembly/Megahit
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Co-assembly succesfully finished" >> ${workingdirectory}/${project}/run.log

#Index the assembly
#echo "Indexing the co-assembly" >> ${workingdirectory}/${project}/run.log
#samtools faidx ${workingdirectory}/${project}/CoAssembly/Megahit/final.contigs.fa
#bwa index ${workingdirectory}/${project}/CoAssembly/Megahit/final.contigs.fa



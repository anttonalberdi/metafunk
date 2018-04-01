#Source settings
source "$metafunkdirectory/settings.sh"

#Create Coassembly directory
mkdir -p ${workingdirectory}/${project}/CoAssembly

#Select source folder from which data will be retrieved
if [[ "$(ls -A ${workingdirectory}/${project}/HumanDNARemoved/)" ]]; then
sourcefolder="HumanDNARemoved"
elif [[ "$(ls -A ${workingdirectory}/${project}/HostDNARemoved/)" ]]; then
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

#Convert to fasta
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Converting FASTQ files to FASTA" >> ${workingdirectory}/${project}/run.log

while read samplefile; do

	#Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  samplefile=$(echo $sample | cut -d ' ' -f2 )
  now=$(date +"%Y-%m-%d %H:%M:%S")

  if [[ $samplefile =~ "/" ]]; then
    #It is PE
		fastq_to_fasta -i ${workingdirectory}/${project}/${sourcefolder}/${samplename}_1.fastq -o ${workingdirectory}/${project}/CoAssembly/${samplename}_1.fasta
		fastq_to_fasta -i ${workingdirectory}/${project}/${sourcefolder}/${samplename}_2.fastq -o ${workingdirectory}/${project}/CoAssembly/${samplename}_2.fasta
		now=$(date +"%Y-%m-%d %H:%M:%S")
		echo "$now | 	${samplename} succesfully converted to fasta" >> ${workingdirectory}/${project}/run.log
	else
		#It is SR
		fastq_to_fasta -i ${workingdirectory}/${project}/${sourcefolder}/${samplename}.fastq -o ${workingdirectory}/${project}/CoAssembly/${samplename}.fasta
		now=$(date +"%Y-%m-%d %H:%M:%S")
		echo "$now | 	Sample ${samplename} succesfully converted to fasta" >> ${workingdirectory}/${project}/run.log
	fi
done < ${sampledatafile}


#Interleave samples if PE  - NOT NECESSARY - TWEAK MEGAHIT CODE
if [[ $samplefile =~ "/" ]]; then
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Interleaving samples" >> ${workingdirectory}/${project}/run.log

while read samplefile; do
	seqtk mergepe ${workingdirectory}/${project}/CoAssembly/${samplename}_1.fasta ${workingdirectory}/${project}/CoAssembly/${samplename}_2.fasta > ${workingdirectory}/${project}/CoAssembly/${samplename}.fasta
	rm ${workingdirectory}/${project}/CoAssembly/${samplename}_[1-2].fasta
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now | 	Sample ${samplename} succesfully interleaved" >> ${workingdirectory}/${project}/run.log
done < ${sampledatafile}
fi

#Concatenate
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Concatenating all samples" >> ${workingdirectory}/${project}/run.log
cat ${workingdirectory}/${project}/CoAssembly/*.fasta > ${workingdirectory}/${project}/CoAssembly/Allsamples.fasta

#Remove co-assembly directory if exists
if [[ $overridecoassembly == "yes" ]]; then
rm -r ${workingdirectory}/${project}/CoAssembly/Megahit
fi

#Kill process if concatenated file does not exists
if [[ ! -f ${workingdirectory}/${project}/CoAssembly/Allsamples.fasta ]]; then
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | ERROR: The concatenated file required for running the co-assembly does not exist" >> ${workingdirectory}/${project}/run.log
exit
fi

#Co-assembly
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Co-assembling all reads" >> ${workingdirectory}/${project}/run.log
if [[ $samplefile =~ "/" ]]; then
	#It is PE
	megahit -t ${threads} --12 ${workingdirectory}/${project}/CoAssembly/Allsamples.fasta -o ${workingdirectory}/${project}/CoAssembly/Megahit
else
	#It is SR
	megahit -t ${threads} -r ${workingdirectory}/${project}/CoAssembly/Allsamples.fasta -o ${workingdirectory}/${project}/CoAssembly/Megahit
fi
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Co-assembly succesfully finished" >> ${workingdirectory}/${project}/run.log

#Index the assembly
if [[ $indexassembly == "yes" ]]; then
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now | Indexing the co-assembly" >> ${workingdirectory}/${project}/run.log
	samtools faidx ${workingdirectory}/${project}/CoAssembly/Megahit/final.contigs.fa
	bwa index ${workingdirectory}/${project}/CoAssembly/Megahit/final.contigs.fa
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now | Co-assembly succesfully indexed" >> ${workingdirectory}/${project}/run.log
fi

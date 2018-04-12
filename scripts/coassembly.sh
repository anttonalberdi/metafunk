#Source settings file
source $settingsfile

#Create Coassembly directory
mkdir -p ${workdir}/CoAssembly

#Select source folder from which data will be retrieved
if [[ "$(ls -A ${workdir}/HumanDNARemoved/)" ]]; then
sourcefolder="HumanDNARemoved"
elif [[ "$(ls -A ${workdir}/HostDNARemoved/)" ]]; then
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

#Convert to fasta
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Converting FASTQ files to FASTA" >> ${workdir}/run_${timestamp}.log

#Declare fastq-to-fasta function
function fastqtofasta() {

sample=${1}
settingsfile=${2}
sourcefolder=${3}

source $settingsfile

	#Obtain data from sample.data.txt columns
  samplename=$(echo $sample | cut -d ' ' -f1 )
  samplefile=$(echo $sample | cut -d ' ' -f2 )
  now=$(date +"%Y-%m-%d %H:%M:%S")

  if [[ $samplefile =~ "/" ]]; then
    #It is PE
    if [ ! -f ${workdir}/CoAssembly/${samplename}_1.fasta ]; then
  		fastq_to_fasta -i ${workdir}/${sourcefolder}/${samplename}_1.fastq -o ${workdir}/CoAssembly/${samplename}_1.fasta
  		fastq_to_fasta -i ${workdir}/${sourcefolder}/${samplename}_2.fastq -o ${workdir}/CoAssembly/${samplename}_2.fasta
      now=$(date +"%Y-%m-%d %H:%M:%S")
  		echo "$now | 	${samplename} succesfully converted to fasta" >> ${workdir}/run_${timestamp}.log
    else
      now=$(date +"%Y-%m-%d %H:%M:%S")
      echo "$now | 	${samplename} fasta files already exists" >> ${workdir}/run_${timestamp}.log
    fi
	else
		#It is SR
    if [ ! -f ${workdir}/CoAssembly/${samplename}.fasta ]; then
  		fastq_to_fasta -i ${workdir}/${sourcefolder}/${samplename}.fastq -o ${workdir}/CoAssembly/${samplename}.fasta
  		now=$(date +"%Y-%m-%d %H:%M:%S")
  		echo "$now | 	Sample ${samplename} succesfully converted to fasta" >> ${workdir}/run_${timestamp}.log
    else
      now=$(date +"%Y-%m-%d %H:%M:%S")
      echo "$now | 	${samplename} fasta file already exists" >> ${workdir}/run_${timestamp}.log
    fi
	fi
}

#Loop in parallel across samples specified in sample.data.txt
export -f fastqtofasta
parallel -j ${threads} -k fastqtofasta {} ${settingsfile} ${sourcefolder} <${sampledatafile}

#Remove co-assembly directory if exists
if [[ $overridecoassembly == "yes" ]]; then
rm -r ${workdir}/CoAssembly/Megahit
fi

#Co-assembly
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Co-assembling all reads" >> ${workdir}/run_${timestamp}.log
if [[ $samplefile =~ "/" ]]; then
	#It is PE
  PE1=$(ls -p ${workdir}/CoAssembly/*_1.fasta | tr '\n' ',' | sed 's/.$//')
  PE2=$(ls -p ${workdir}/CoAssembly/*_2.fasta | tr '\n' ',' | sed 's/.$//')
	megahit -t ${threads} -1 ${PE1} -2 ${PE2} -o ${workdir}/CoAssembly/Megahit
else
	#It is SR
  SR=$(ls -p ${workdir}/CoAssembly/*.fasta | tr '\n' ',' | sed 's/.$//')
	megahit -t ${threads} -r ${SR} -o ${workdir}/CoAssembly/Megahit
fi
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Co-assembly succesfully finished" >> ${workdir}/run_${timestamp}.log

#Index the assembly
if [[ $indexassembly == "yes" ]]; then
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now | Indexing the co-assembly" >> ${workdir}/run_${timestamp}.log
	samtools faidx ${workdir}/CoAssembly/Megahit/final.contigs.fa
	bwa index ${workdir}/CoAssembly/Megahit/final.contigs.fa
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now | Co-assembly succesfully indexed" >> ${workdir}/run_${timestamp}.log
fi

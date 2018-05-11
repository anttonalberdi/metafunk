#Source settings file
source $settingsfile


#Create GeneMapping directory
mkdir -p ${workdir}/GeneMapping

#Index genes
if [ ! -f ${workdir}/GenePrediction/assembly.genes.fna.fai ]; then
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now |        Indexing gene catalogue" >> ${workdir}/run_${timestamp}.log
	samtools faidx ${workdir}/GenePrediction/assembly.genes.fna
	bwa index ${workdir}/GenePrediction/assembly.genes.fna
	now=$(date +"%Y-%m-%d %H:%M:%S")
	if [ -f ${workdir}/GenePrediction/assembly.genes.fna.fai ]; then
	echo "$now |        Gene catalogue succesfully indexed" >> ${workdir}/run_${timestamp}.log
	else
	echo "$now |        There was an error while indexing the gene catalogue" >> ${workdir}/run_${timestamp}.log
	fi
else
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now |        Gene catalogue is already indexed" >> ${workdir}/run_${timestamp}.log
fi

#Get gene lengths
if [ ! -f ${workdir}/GenePrediction/assembly.genes.lengths ]; then
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now |        Calculating gene lengths" >> ${workdir}/run_${timestamp}.log
python ${metafunkdirectory}/scripts/contiglengths.py -i ${workdir}/GenePrediction/assembly.genes.fna > ${workdir}/GenePrediction/assembly.genes.lengths
filesize=$(ls -l ${workdir}/GenePrediction/assembly.genes.lengths | awk '{print $5}')
now=$(date +"%Y-%m-%d %H:%M:%S")
	if [[ $filesize > 0 ]]; then
	echo "$now |        Gene length file was successfully created" >> ${workdir}/run_${timestamp}.log
	else
	echo "$now |        There was an error while generating the gene length file" >> ${workdir}/run_${timestamp}.log
	fi
else
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now |        Gene lengths are already computed" >> ${workdir}/run_${timestamp}.log

fi

#Select source folder from which data will be retrieved
if [[ "$(ls -A ${workdir}/HumanDNARemoved/)" ]]; then
sourcefolder="HumanDNARemoved"
repair="no"
elif [[ "$(ls -A ${workdir}/HostDNARemoved/)" ]]; then
sourcefolder="HostDNARemoved"
repair="no"
elif [[ "$(ls -A ${workdir}/LowComplexFiltered/)" ]]; then
sourcefolder="LowComplexFiltered"
repair="yes"
elif [[ "$(ls -A ${workdir}/DuplicatesRemoved/)" ]]; then
sourcefolder="DuplicatesRemoved"
repair="yes"
elif [[ "$(ls -A ${workdir}/QualityFiltered/)" ]]; then
sourcefolder="QualityFiltered"
repair="yes"
else
sourcefolder="RawData"
repair="no"
fi


#Map reads back to genes
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Mapping reads from directory $sourcefolder to gene catalogue" >> ${workdir}/run_${timestamp}.log

while read sample; do

		#Obtain data from sample.data.txt columns and get file name
		samplename=$(echo $sample | cut -d ' ' -f1)
		sampleinfo=$(echo $sample | cut -d ' ' -f2)

		if [[ $sampleinfo =~ "/" ]]; then
			#It is PE
			#Repair reads
			if [[ $repair == "yes" ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
				echo "$now | 		Repairing sample ${samplename}" >> ${workdir}/run_${timestamp}.log
				repair.sh in=${workdir}/${sourcefolder}/${samplename}_1.fastq in2=${workdir}/${sourcefolder}/${samplename}_2.fastq out=${workdir}/GeneMapping/${samplename}_1.fastq out2=${workdir}/GeneMapping/${samplename}_2.fastq
				sourcefolder="GeneMapping"
			fi
			#Map reads against the gene catalogue
			now=$(date +"%Y-%m-%d %H:%M:%S")
			echo "$now | 		Mapping $samplename reads" >> ${workdir}/run_${timestamp}.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/GenePrediction/assembly.genes.fna ${workdir}/${sourcefolder}/${samplename}_1.fastq ${workdir}/${sourcefolder}/${samplename}_2.fastq | samtools view -b -F0x4 - > ${workdir}/GeneMapping/${samplename}.bam
			#Check if output file has been created; otherwise, print error message and kill the job
			if [[ ! -s ${workdir}/GeneMapping/${samplename}.bam ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
				echo "$now | 		ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			fi
			#Generate coverage table
			samtools sort -T ${workdir}/GeneMapping/${samplename}.tmp.bam -o ${workdir}/GeneMapping/${samplename}.sorted.bam ${workdir}/GeneMapping/${samplename}.bam
			samtools flagstat ${workdir}/GeneMapping/${samplename}.bam > ${workdir}/GeneMapping/${samplename}.flagstat
				echo "$now | 		Calculating coverage for $samplename" >> ${workdir}/run_${timestamp}.log
			bedtools genomecov -ibam ${workdir}/GeneMapping/${samplename}.sorted.bam -g ${workdir}/GenePrediction/assembly.genes.lengths > ${workdir}/GeneMapping/${samplename}.cov
			awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' ${workdir}/GeneMapping/${samplename}.cov > ${workdir}/GeneMapping/${samplename}.cov.csv
			now=$(date +"%Y-%m-%d %H:%M:%S")
			if [[ ! -s ${workdir}/GeneMapping/${samplename}.cov.csv ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
				echo "$now | 		ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			else
				echo "$now | 		Coverage of $samplename successfully calculated" >> ${workdir}/run_${timestamp}.log
			fi
		else
			#It is SR
			#Map reads against the gene catalogue
			now=$(date +"%Y-%m-%d %H:%M:%S")
			echo "$now | 			Mapping $samplename reads" >> ${workdir}/run_${timestamp}.log
			bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/GenePrediction/assembly.genes.fna ${workdir}/${sourcefolder}/${samplename}.fastq | samtools view -b -F0x4 - > ${workdir}/GeneMapping/${samplename}.bam
			#Check if output file has been created; otherwise, print error message and kill the job
			if [[ ! -s ${workdir}/GeneMapping/${samplename}.bam ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			fi
			#Generate coverage table
			samtools sort -T ${workdir}/GeneMapping/${samplename}.tmp.bam -o ${workdir}/GeneMapping/${samplename}.sorted.bam ${workdir}/GeneMapping/${samplename}.bam
			samtools flagstat ${workdir}/GeneMapping/${samplename}.bam > ${workdir}/GeneMapping/${samplename}.flagstat
				echo "$now |      Calculating coverage for $samplename" >> ${workdir}/run_${timestamp}.log
			bedtools genomecov -ibam ${workdir}/GeneMapping/${samplename}.sorted.bam -g ${workdir}/GenePrediction/assembly.genes.lengths > ${workdir}/GeneMapping/${samplename}.cov
			awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' ${workdir}/GeneMapping/${samplename}.cov > ${workdir}/GeneMapping/${samplename}.cov.csv
			now=$(date +"%Y-%m-%d %H:%M:%S")
			if [[ ! -s ${workdir}/GeneMapping/${samplename}.cov.csv ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			else
				echo "$now |      Coverage of $samplename successfully calculated" >> ${workdir}/run_${timestamp}.log
			fi
		fi
done < ${sampledatafile}

#Collate coverages
now=$(date +"%Y-%m-%d %H:%M:%S")
mkdir -p ${workdir}/GeneTables
echo "$now |      Collating coverage files" >> ${workdir}/run_${timestamp}.log
perl ${metafunkdirectory}/scripts/collatecoverages.pl ${workdir}/GeneMapping/ > ${workdir}/GeneTables/GeneCoverageTable.csv

#Generate hit table
echo "$now | Generating hit table" >> ${workdir}/run_${timestamp}.log
export WORKDIR="${workdir}"
Rscript ${metafunkdirectory}/scripts/createhittable.r --no-save
filesize=$(ls -l ${workdir}/GeneTables/HitTable.csv | awk '{print $5}')
now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $filesize > 0 ]]; then
echo "$now |        Hit table was successfully created" >> ${workdir}/run_${timestamp}.log
else
echo "$now |        There was an error while generating the hit table" >> ${workdir}/run_${timestamp}.log
fi

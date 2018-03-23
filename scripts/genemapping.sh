#Source dependencies
source "$metafunkdirectory/settings.sh"

mkdir -p ${workingdirectory}/${project}/GeneMapping

#Index genes
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now |        Indexing gene catalogue" >> ${workingdirectory}/${project}/run.log
samtools faidx ${workingdirectory}/${project}/GenePrediction/assembly.genes.fna
bwa index ${workingdirectory}/${project}/GenePrediction/assembly.genes.fna
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now |        Gene catalogue indexed" >> ${workingdirectory}/${project}/run.log

#Get gene lengths
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now |        Calculating gene lengths" >> ${workingdirectory}/${project}/run.log
python ${metafunkdirectory}/scripts/contiglengths.py -i ${workingdirectory}/${project}/GenePrediction/assembly.genes.fna > ${workingdirectory}/${project}/GenePrediction/assembly.genes.lengths
filesize=$(ls -l ${workingdirectory}/${project}/GenePrediction/assembly.genes.lengths | awk '{print $5}')
now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $filesize > 0 ]]; then
echo "$now |        Gene length file was successfully created" >> ${workingdirectory}/${project}/run.log
else
echo "$now |        There was an error while generating the gene length file" >> ${workingdirectory}/${project}/run.log
fi

mkdir -p ${workingdirectory}/${project}/GeneMapping

#Map reads back to genes
if [[ $seqtype == "SR" ]]; then
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Mapping SR reads to the gene catalogue" >> ${workingdirectory}/${project}/run.log
	#Iterate through samples
	while read samplefile; do
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now |      Mapping sample $samplefile" >> ${workingdirectory}/${project}/run.log
		#Select LowComplex or Original data file
		if [[ $removelowcomplexity == "yes" ]]; then
		inputfile=${workingdirectory}/${project}/LowComplexFiltered/${samplefile}.fastq
		else
		inputfile=${workingdirectory}/${project}/RawData/${samplefile}.fastq
		fi
		bwa mem -t ${threads} -R '@RG\tID:Project\tCN:User\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workingdirectory}/${project}/GenePrediction/assembly.genes.fna ${inputfile} | samtools view -b -F0x4 - > ${workingdirectory}/${project}/GeneMapping/${samplefile}.bam
		samtools sort -T ${workingdirectory}/${project}/GeneMapping/${samplefile}.tmp.bam -o ${workingdirectory}/${project}/GeneMapping/${samplefile}.sorted.bam ${workingdirectory}/${project}/GeneMapping/${samplefile}.bam
		samtools flagstat ${workingdirectory}/${project}/GeneMapping/${samplefile}.bam > ${workingdirectory}/${project}/GeneMapping/${samplefile}.flagstat
		bedtools genomecov -ibam ${workingdirectory}/${project}/GeneMapping/${samplefile}.sorted.bam -g ${1}/9.0-GenePrediction/assembly.lengths > ${workingdirectory}/${project}/GeneMapping/${samplefile}.cov
		awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' ${workingdirectory}/${project}/GeneMapping/${samplefile}.cov > ${workingdirectory}/${project}/GeneMapping/${samplefile}.cov.csv
echo "$now |      Sample $samplefile successfully mapped" >> ${workingdirectory}/${project}/run.log
	done < ${metafunkdirectory}/sample.data.txt

elif [[ $seqtype == "PE" ]]; then
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Removing host DNA from PE data" >> ${workingdirectory}/${project}/run.log

else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Sequencing read type has not been specified. It needs to be either SR or PE" >> run.log
fi

#Collate coverages
perl ${metafunkdirectory}/scripts/collatecoverages.pl ${workingdirectory}/${project}/GeneMapping/ > ${workingdirectory}/${project}/GeneCoverageTable.csv

#Generate hit table
echo "$now | Generating hit table" >> ${workingdirectory}/${project}/run.log
export WORKDIR="${workingdirectory}/${project}"
Rscript ${workdir}/scripts/createhittable.r --no-save
filesize=$(ls -l ${workingdirectory}/${project}/HitTable.csv | awk '{print $5}')
now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $filesize > 0 ]]; then
echo "$now |        Hit table was successfully created" >> ${workingdirectory}/${project}/run.log
else
echo "$now |        There was an error while generating the hit table" >> ${workingdirectory}/${project}/run.log
fi

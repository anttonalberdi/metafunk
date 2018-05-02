#Source settings file
source $settingsfile


#Create ContigMapping directory
mkdir -p ${workdir}/ContigMapping

#Cut-up assembly
if [[ $cutcontigs == "yes" ]]; then
  python ${metafunkdirectory}/scripts/cutcontigs.py -c ${chunksize} -o 0 -m ${workdir}/CoAssembly/Megahit/final.contigs.fa > ${workdir}/CoAssembly/Megahit/final.contigs.cut.fa
  python ${metafunkdirectory}/scripts/filtercontiglength.py -m ${mincontigsize} ${workdir}/CoAssembly/Megahit/final.contigs.cut.fa > ${workdir}/CoAssembly/Megahit/final.contigs.cut.filt.fa
fi

#Index contigs
if [[ $cutcontigs == "yes" ]]; then
  if [ ! -f ${workdir}/CoAssembly/Megahit/final.contigs.cut.filt.fa.fai ]; then
    echo "$now |        Indexing chunked assembly" >> ${workdir}/run_${timestamp}.log
    samtools faidx ${workdir}/CoAssembly/Megahit/final.contigs.cut.filt.fa
  	bwa index ${workdir}/CoAssembly/Megahit/final.contigs.cut.filt.fa
  	now=$(date +"%Y-%m-%d %H:%M:%S")
  	echo "$now |        Chunked assembly succesfully indexed" >> ${workdir}/run_${timestamp}.log
  else
    echo "$now |        Chunked assembly is already indexed" >> ${workdir}/run_${timestamp}.log
  fi
else
  if [ ! -f ${workdir}/CoAssembly/Megahit/final.contigs.fa.fai ]; then
    echo "$now |        Indexing assembly" >> ${workdir}/run_${timestamp}.log
    samtools faidx ${workdir}/CoAssembly/Megahit/final.contigs.fa
    bwa index ${workdir}/CoAssembly/Megahit/final.contigs.fa
    now=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$now |        Assembly succesfully indexed" >> ${workdir}/run_${timestamp}.log
  else
    echo "$now |        Assembly is already indexed" >> ${workdir}/run_${timestamp}.log
  fi
fi

#Get contig lengths
if [ ! -f ${workdir}/ContigMapping/contig.lengths ]; then
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now |        Calculating contig lengths" >> ${workdir}/run_${timestamp}.log
if [[ $cutcontigs == "yes" ]]; then
  python ${metafunkdirectory}/scripts/contiglengths.py -i ${workdir}/CoAssembly/Megahit/final.contigs.cut.filt.fa > ${workdir}/ContigMapping/contig.lengths
  else
  python ${metafunkdirectory}/scripts/contiglengths.py -i ${workdir}/CoAssembly/Megahit/final.contigs.fa > ${workdir}/ContigMapping/contig.lengths
fi
filesize=$(ls -l ${workdir}/ContigMapping/contig.lengths | awk '{print $5}')
now=$(date +"%Y-%m-%d %H:%M:%S")
	if [[ $filesize > 0 ]]; then
	echo "$now |        Contig length file was successfully created" >> ${workdir}/run_${timestamp}.log
	else
	echo "$now |        There was an error while generating the contig length file" >> ${workdir}/run_${timestamp}.log
	fi
else
	now=$(date +"%Y-%m-%d %H:%M:%S")
	echo "$now |        Contig lengths are already computed" >> ${workdir}/run_${timestamp}.log
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

#Map reads back to contigs
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Mapping reads from directory $sourcefolder to contigs" >> ${workdir}/run_${timestamp}.log

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
				repair.sh in=${workdir}/${sourcefolder}/${samplename}_1.fastq in2=${workdir}/${sourcefolder}/${samplename}_2.fastq out=${workdir}/ContigMapping/${samplename}_1.fastq out2=${workdir}/ContigMapping/${samplename}_2.fastq
				sourcefolder="ContigMapping"
			fi
			#Map reads against the gene catalogue
			now=$(date +"%Y-%m-%d %H:%M:%S")
			echo "$now | 		Mapping $samplename reads" >> ${workdir}/run_${timestamp}.log
      if [[ $cutcontigs == "yes" ]]; then
			  bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/CoAssembly/Megahit/final.contigs.cut.filt.fa ${workdir}/${sourcefolder}/${samplename}_1.fastq ${workdir}/${sourcefolder}/${samplename}_2.fastq | samtools view -b -F0x4 - > ${workdir}/ContigMapping/${samplename}.bam
        else
        bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/CoAssembly/Megahit/final.contigs.fa ${workdir}/${sourcefolder}/${samplename}_1.fastq ${workdir}/${sourcefolder}/${samplename}_2.fastq | samtools view -b -F0x4 - > ${workdir}/ContigMapping/${samplename}.bam
      fi
      #Check if output file has been created; otherwise, print error message and kill the job
			if [[ ! -s ${workdir}/ContigMapping/${samplename}.bam ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
				echo "$now | 		ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			fi
			#Generate coverage table
			samtools sort -T ${workdir}/ContigMapping/${samplename}.tmp.bam -o ${workdir}/ContigMapping/${samplename}.sorted.bam ${workdir}/ContigMapping/${samplename}.bam
			samtools flagstat ${workdir}/ContigMapping/${samplename}.bam > ${workdir}/ContigMapping/${samplename}.flagstat
				echo "$now | 		Calculating coverage for $samplename" >> ${workdir}/run_${timestamp}.log
			bedtools genomecov -ibam ${workdir}/ContigMapping/${samplename}.sorted.bam -g ${workdir}/ContigMapping/contig.lengths > ${workdir}/ContigMapping/${samplename}.cov
			awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' ${workdir}/ContigMapping/${samplename}.cov > ${workdir}/ContigMapping/${samplename}.cov.csv
			now=$(date +"%Y-%m-%d %H:%M:%S")
			if [[ ! -s ${workdir}/ContigMapping/${samplename}.cov.csv ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
				echo "$now | 		ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			else
				echo "$now | 		Coverage of $samplefile successfully calculated" >> ${workdir}/run_${timestamp}.log
			fi
		else
			#It is SR
			#Map reads against the gene catalogue
			now=$(date +"%Y-%m-%d %H:%M:%S")
			echo "$now | 			Mapping $samplename reads" >> ${workdir}/run_${timestamp}.log
      if [[ $cutcontigs == "yes" ]]; then
        bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/CoAssembly/Megahit/final.contigs.cut.filt.fa ${workdir}/${sourcefolder}/${samplename}.fastq | samtools view -b -F0x4 - > ${workdir}/ContigMapping/${samplename}.bam
        else
        bwa mem -t ${threads} -R '@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample' ${workdir}/CoAssembly/Megahit/final.contigs.fa ${workdir}/${sourcefolder}/${samplename}.fastq | samtools view -b -F0x4 - > ${workdir}/ContigMapping/${samplename}.bam
      fi
			#Check if output file has been created; otherwise, print error message and kill the job
			if [[ ! -s ${workdir}/ContigMapping/${samplename}.bam ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			fi
			#Generate coverage table
			samtools sort -T ${workdir}/ContigMapping/${samplename}.tmp.bam -o ${workdir}/ContigMapping/${samplename}.sorted.bam ${workdir}/ContigMapping/${samplename}.bam
			samtools flagstat ${workdir}/ContigMapping/${samplename}.bam > ${workdir}/ContigMapping/${samplename}.flagstat
				echo "$now |      Calculating coverage for $samplename" >> ${workdir}/run_${timestamp}.log
			bedtools genomecov -ibam ${workdir}/ContigMapping/${samplename}.sorted.bam -g ${workdir}/ContigMapping/contig.lengths > ${workdir}/ContigMapping/${samplename}.cov
			awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' ${workdir}/ContigMapping/${samplename}.cov > ${workdir}/ContigMapping/${samplename}.cov.csv
			now=$(date +"%Y-%m-%d %H:%M:%S")
			if [[ ! -s ${workdir}/ContigMapping/${samplename}.cov.csv ]]; then
				now=$(date +"%Y-%m-%d %H:%M:%S")
				echo "$now | 			ERROR: There was an error when mapping sample $samplename" >> ${workdir}/run_${timestamp}.log
				exit
			else
				echo "$now |      Coverage of $samplefile successfully calculated" >> ${workdir}/run_${timestamp}.log
			fi
		fi
done < ${sampledatafile}

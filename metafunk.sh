#Source dependencies
source "${1}/settings.sh"

#########
# Create and set working directory 
#########

mkdir -p ${workingdirectory}
cd ${workingdirectory}
mkdir ${project}
cd ${project}
projectdirectory=${workingdirectory}/${project}

#########
# Print settings to log file
#########

echo "##### SETTINGS ####" > ${projectdirectory}/run.log
echo "Number of threads: $threads" >> ${projectdirectory}/run.log
echo "Sequencing read type: $seqtype" >> ${projectdirectory}/run.log
echo "Sequencing read length: $readlength" >> ${projectdirectory}/run.log
echo '###################/n' >> ${projectdirectory}/run.log

#########
# Copy data to project directory
#########

if [[ $copydata == "yes" ]]; then  
mkdir ${projectdirectory}/RawData
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Copying and uncompressing data files" >> ${projectdirectory}/run.log
while read samplefile; do 
cp ${datadirectory}/${samplefile}.fastq.gz ${projectdirectory}/RawData/
#Uncompress files
pigz -d -p ${threads} ${projectdirectory}/RawData/${samplefile}.fastq.gz
done < ${metafunkdirectory}/sample.data.txt
#Print stats
filenumber=$(ls ${projectdirectory}/RawData/| wc -l)
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | $filenumber files were copied and uncompressed" >> ${projectdirectory}/run.log
else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Data will not be copied and uncompressed" >> ${projectdirectory}/run.log
fi

#########
# Remove low complexity reads
#########

if [[ $removelowcomplexity == "yes" ]]; then  
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Running low complexity filter" >> ${projectdirectory}/run.log
export metafunkdirectory
sh ${metafunkdirectory}/scripts/lowcomplexity.sh
else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Low complexity will not be filtered" >> ${projectdirectory}/run.log
fi

#########
# Remove host DNA
#########

if [[ $removehostdna == "yes" ]]; then  
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Removing host DNA" >> ${projectdirectory}/run.log
export metafunkdirectory
sh ${metafunkdirectory}/scripts/removehostdna.sh
else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Host DNA will not be removed" >> ${projectdirectory}/run.log
fi

#########
# Perform co-assembly
#########

if [[ $coassembly == "yes" ]]; then  
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Co-assembling reads" >> ${projectdirectory}/run.log
export metafunkdirectory
sh ${metafunkdirectory}/scripts/coassembly.sh
else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Reads will not be co-assembled" >> ${projectdirectory}/run.log
fi

#########
# Predict genes
#########

if [[ $geneprediction == "yes" ]]; then  
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Predicting genes" >> ${projectdirectory}/run.log
export metafunkdirectory
sh ${metafunkdirectory}/scripts/geneprediction.sh
else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Gene prediction will not be performed" >> ${projectdirectory}/run.log
fi

#########
# Map reads back to the genes
#########

if [[ $genemapping == "yes" ]]; then  
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Mapping reads back to genes" >> ${projectdirectory}/run.log
export metafunkdirectory
sh ${metafunkdirectory}/scripts/genemapping.sh
else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Gene mapping will not be performed" >> ${projectdirectory}/run.log
fi


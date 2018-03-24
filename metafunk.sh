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
echo "###################" >> ${projectdirectory}/run.log
echo "" >> ${projectdirectory}/run.log

#########
# Check dependencies
#########

export metafunkdirectory
sh ${metafunkdirectory}/scripts/checkdependencies.sh

#########
# Copy data to project directory
#########
now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $copydata == "yes" ]]; then
echo "$now | Copying and uncompressing data files" >> ${projectdirectory}/run.log
export metafunkdirectory
sh ${metafunkdirectory}/scripts/transferdata.sh
else
echo "$now | Data will not be copied and uncompressed" >> ${projectdirectory}/run.log
fi

#########
# Quality filtering
#########

now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $qualityfiltering == "yes" ]]; then
echo "$now | Performing quality filtering" >> ${projectdirectory}/run.log
export metafunkdirectory
sh ${metafunkdirectory}/scripts/qualityfiltering.sh
else
echo "$now | Quality filtering will not be performed" >> ${projectdirectory}/run.log
fi


#########
# Remove duplicates
#########

if [[ $removeduplicates == "yes" ]]; then
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Removing duplicates" >> ${projectdirectory}/run.log
export metafunkdirectory
sh ${metafunkdirectory}/scripts/removeduplicates.sh
else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Duplicates will not be removed" >> ${projectdirectory}/run.log
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
# Map reads back to the genes and generate Coverage and Hit tables
#########
now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $genemapping == "yes" ]]; then
echo "$now | Mapping reads back to genes" >> ${projectdirectory}/run.log
export metafunkdirectory
sh ${metafunkdirectory}/scripts/genemapping.sh
else
echo "$now | Gene mapping will not be performed" >> ${projectdirectory}/run.log
fi

#########
# Normalise Coverage and Hit tables
#########
now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $tss == "yes" || $css == "yes" ]]; then
  echo "$now | Normalising hit and coverage tables" >> ${projectdirectory}/run.log
  export metafunkdirectory
  sh ${metafunkdirectory}/scripts/normalisetables.sh
  else
  echo "$now | Hit and coverage tables will not be normalised" >> ${projectdirectory}/run.log
fi

#########
# Perform functional annotation
#########
now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $kegg == "yes" ]]; then
  echo "$now | Starting KEGG functional annotation" >> ${projectdirectory}/run.log
  export metafunkdirectory
  sh ${metafunkdirectory}/scripts/functionalannotation_kegg.sh
  else
  echo "$now | KEGG functional annotation will not be performed" >> ${projectdirectory}/run.log
fi

now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $eggnog == "yes" ]]; then
  echo "$now | Starting EggNog functional annotation" >> ${projectdirectory}/run.log
  export metafunkdirectory
  sh ${metafunkdirectory}/scripts/functionalannotation_eggnog.sh
  else
  echo "$now | EggNog functional annotation will not be performed" >> ${projectdirectory}/run.log
fi

#MetaFunk version
version="0.1"

#########
# Get options
#########

while getopts w:d:s:f:t: option; do
  case "${option}"
  in
  w) workdir=${OPTARG};;
  d) sampledatafile=${OPTARG};;
  s) settingsfile=${OPTARG};;
  f) datadir=${OPTARG};;
  t) threads=${OPTARG};;
  esac
done

#Get timestamp

timestamp=$(date +"%Y%d%m_%H%M%S")

#Get metafunk directory
metafunkdirectory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#Source dependencies
source "${settingsfile}"

#########
# Create and set working directory
#########

mkdir -p ${workdir}
cd ${workdir}

#########
# Print settings and sample information to log file
#########
echo " " > ${workdir}/run_${timestamp}.log
echo "Running MetaFunk v$version pipeline with the following settings and samples:" > ${workdir}/run_${timestamp}.log
echo " " >> ${workdir}/run_${timestamp}.log
echo "##### SETTINGS ####" >> ${workdir}/run_${timestamp}.log
echo "Metafunk directory: $metafunkdirectory" >> ${workdir}/run_${timestamp}.log
echo "Working directory: $workdir" >> ${workdir}/run_${timestamp}.log
echo "Sample data file: $sampledatafile" >> ${workdir}/run_${timestamp}.log
echo "Settings file: $settingsfile" >> ${workdir}/run_${timestamp}.log
echo "Data directory: $datadir" >> ${workdir}/run_${timestamp}.log
echo "Number of threads: $threads" >> ${workdir}/run_${timestamp}.log
echo "Sequencing platform: $platform" >> ${workdir}/run_${timestamp}.log
echo " " >> ${workdir}/run_${timestamp}.log
echo "##### SAMPLES ####" >> ${workdir}/run_${timestamp}.log
while read sample; do
  samplename=$(echo $sample | cut -d ' ' -f1 )
  sampleinfo=$(echo $sample | cut -d ' ' -f2 )
  if [[ $sampleinfo =~ "/" && ! $sampleinfo =~ ";" ]]; then
  echo "  $samplename PE SF" >> ${workdir}/run_${timestamp}.log
  elif [[ $sampleinfo =~ "/" && $sampleinfo =~ ";" ]]; then
  echo "  $samplename PE MF" >> ${workdir}/run_${timestamp}.log
  elif [[ ! $sampleinfo =~ "/" && $sampleinfo =~ ";" ]]; then
  echo "  $samplename SR MF" >> ${workdir}/run_${timestamp}.log
  else
  echo "  $samplename SR SF" >> ${workdir}/run_${timestamp}.log
  fi
done < ${sampledatafile}
echo "###################" >> ${workdir}/run_${timestamp}.log
echo "" >> ${workdir}/run_${timestamp}.log

#########
# Check sample.data file
#########
# Check if the number of columns is as expected in all rows
cat ${sampledatafile} | awk -F ' ' -v NCOLS=3 'NF!=NCOLS{printf "Wrong number of columns at line %d\n", NR; exit}'

#########
# Check dependencies
#########

export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/checkdependencies.sh

#########
# Transfer data to project directory
#########
now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $copydata == "yes" ]]; then
samplenumber=$(cat sample.data.txt | wc -l)
echo "$now | Copying and uncompressing data files of $samplenumber samples" >> ${workdir}/run_${timestamp}.log
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/transferdata.sh
else
echo "$now | Data will not be copied and uncompressed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Quality filtering
#########

now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $qualityfiltering == "yes" ]]; then
echo "$now | Performing quality filtering" >> ${workdir}/run_${timestamp}.log
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/qualityfiltering.sh
else
echo "$now | Quality filtering will not be performed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Remove duplicates
#########

if [[ $removeduplicates == "yes" ]]; then
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/removeduplicates.sh
else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Duplicates will not be removed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Remove low complexity reads
#########

if [[ $removelowcomplexity == "yes" ]]; then
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/lowcomplexity.sh
else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Low complexity will not be filtered" >> ${workdir}/run_${timestamp}.log
fi

#########
# Remove host DNA
#########

if [[ $removehostdna == "yes" ]]; then
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Removing host DNA" >> ${workdir}/run_${timestamp}.log
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/removehostdna.sh
else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Host DNA will not be removed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Perform co-assembly
#########

if [[ $coassembly == "yes" ]]; then
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Co-assembling reads" >> ${workdir}/run_${timestamp}.log
export metafunkdirectory
sh ${metafunkdirectory}/scripts/coassembly.sh
else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Reads will not be co-assembled" >> ${workdir}/run_${timestamp}.log
fi

#########
# Predict genes
#########

if [[ $geneprediction == "yes" ]]; then
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Predicting genes" >> ${workdir}/run_${timestamp}.log
export metafunkdirectory
sh ${metafunkdirectory}/scripts/geneprediction.sh
else
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Gene prediction will not be performed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Map reads back to the genes and generate Coverage and Hit tables
#########
now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $genemapping == "yes" ]]; then
echo "$now | Mapping reads back to genes" >> ${workdir}/run_${timestamp}.log
export metafunkdirectory
sh ${metafunkdirectory}/scripts/genemapping.sh
else
echo "$now | Gene mapping will not be performed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Normalise Coverage and Hit tables
#########
now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $tss == "yes" || $css == "yes" ]]; then
  echo "$now | Normalising hit and coverage tables" >> ${workdir}/run_${timestamp}.log
  export metafunkdirectory
  sh ${metafunkdirectory}/scripts/normalisetables.sh
  else
  echo "$now | Hit and coverage tables will not be normalised" >> ${workdir}/run_${timestamp}.log
fi

#########
# Perform functional annotation
#########
now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $kegg == "yes" ]]; then
  echo "$now | Starting KEGG functional annotation" >> ${workdir}/run_${timestamp}.log
  export metafunkdirectory
  sh ${metafunkdirectory}/scripts/functionalannotation_kegg.sh
  else
  echo "$now | KEGG functional annotation will not be performed" >> ${workdir}/run_${timestamp}.log
fi

now=$(date +"%Y-%d-%m %H:%M:%S")
if [[ $eggnog == "yes" ]]; then
  echo "$now | Starting EggNog functional annotation" >> ${workdir}/run_${timestamp}.log
  export metafunkdirectory
  sh ${metafunkdirectory}/scripts/functionalannotation_eggnog.sh
  else
  echo "$now | EggNog functional annotation will not be performed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Final message
#########
now=$(date +"%Y-%d-%m %H:%M:%S")
echo "$now | Congratulations, the pipeline was succesfully run" >> ${workdir}/run_${timestamp}.log

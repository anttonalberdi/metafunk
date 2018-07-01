#MetaFunk version
version="0.1"

#########
# Get options
#########

while getopts w:d:s:f:t:m: option; do
  case "${option}"
  in
  w) workdir=${OPTARG};;
  d) sampledatafile=${OPTARG};;
  s) settingsfile=${OPTARG};;
  f) datadir=${OPTARG};;
  t) threads=${OPTARG};;
  m) modules=${OPTARG};;
  \?) echo "Invalid option: -$OPTARG"
  exit ;;
  :) echo "Option -$OPTARG requires an argument."
  exit ;;
  esac
done

#Get timestamp
timestamp=$(date +"%Y%m%d_%H%M%S")

#Get metafunk directory
metafunkdirectory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#Source dependencies
source "${settingsfile}"

#########
# Process module data
#########

if [[ $modules == "1" || $modules =~ ",1," || $modules == 1,* || $modules == *,1 ]];
  then copydata="yes"
  else copydata="no"
fi
if [[ $modules == "2" || $modules =~ ",2," || $modules == 2,* || $modules == *,2 ]];
  then qualityfiltering="yes"
  else qualityfiltering="no"
fi
if [[ $modules == "3" || $modules =~ ",3," || $modules == 3,* || $modules == *,3 ]];
  then removeduplicates="yes"
  else removeduplicates="no"
fi
if [[ $modules == "4" || $modules =~ ",4," || $modules == 4,* || $modules == *,4 ]];
  then removelowcomplexity="yes"
  else removelowcomplexity="no"
fi
if [[ $modules == "5" || $modules =~ ",5," || $modules == 5,* || $modules == *,5 ]];
  then removehostdna="yes"
  else removehostdna="no"
fi
if [[ $modules == "6" || $modules =~ ",6," || $modules == 6,* || $modules == *,6 ]];
  then removehumandna="yes"
  else removehumandna="no"
fi
if [[ $modules == "7" || $modules =~ ",7," || $modules == 7,* || $modules == *,7 ]];
  then coassembly="yes"
  else coassembly="no"
fi
if [[ $modules == "8" || $modules =~ ",8," || $modules == 8,* || $modules == *,8 ]];
  then geneprediction="yes"
  else geneprediction="no"
fi
if [[ $modules == "9" || $modules =~ ",9," || $modules == 9,* || $modules == *,9 ]];
  then genemapping="yes"
  else genemapping="no"
fi
if [[ $modules == "10" || $modules =~ ",10," || $modules == 10,* || $modules == *,10 ]];
  then functional="yes"
  else functional="no"
fi
if [[ $modules == "11" || $modules =~ ",11," || $modules == 11,* || $modules == *,11 ]];
  then taxonomic="yes"
  else taxonomic="no"
fi
if [[ $modules == "12" || $modules =~ ",12," || $modules == 12,* || $modules == *,12 ]];
  then contigmapping="yes"
  else contigmapping="no"
fi

if [ -z "$modules" ]; then
  echo "No modules were specified."
  exit
fi

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
cat ${sampledatafile} | awk -F ' ' -v NCOLS=4 'NF!=NCOLS{printf "Wrong number of columns at line %d\n", NR; exit}'

#########
# Check dependencies
#########

export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/checkdependencies.sh

#########
# Transfer data to project directory
#########
now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $copydata == "yes" ]]; then
samplenumber=$(cat sample.data.txt | wc -l)
echo "$now | DATA TRANSFER" >> ${workdir}/run_${timestamp}.log
echo "$now | Copying and uncompressing data files of $samplenumber samples" >> ${workdir}/run_${timestamp}.log
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/transferdata.sh
else
echo "$now | Data will not be copied and uncompressed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Quality filtering
#########

now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $qualityfiltering == "yes" ]]; then
echo "$now | QUALITY FILTERING" >> ${workdir}/run_${timestamp}.log
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/qualityfiltering.sh
else
echo "$now | Quality filtering will not be performed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Remove duplicates
#########

now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $removeduplicates == "yes" ]]; then
echo "$now | DUPLICATE REMOVAL" >> ${workdir}/run_${timestamp}.log
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/removeduplicates.sh
else
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Duplicates will not be removed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Remove low complexity reads
#########

now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $removelowcomplexity == "yes" ]]; then
echo "$now | LOW COMPLEXITY FILTERING" >> ${workdir}/run_${timestamp}.log
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/lowcomplexity.sh
else
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Low complexity will not be filtered" >> ${workdir}/run_${timestamp}.log
fi

#########
# Remove host DNA
#########

if [[ $removehostdna == "yes" ]]; then
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | HOST DNA REMOVAL" >> ${workdir}/run_${timestamp}.log
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/removehostdna.sh
else
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Host DNA will not be removed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Remove human DNA
#########

if [[ $removehumandna == "yes" ]]; then
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | HUMAN DNA REMOVAL" >> ${workdir}/run_${timestamp}.log
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/removehumandna.sh
else
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Human DNA will not be removed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Perform co-assembly
#########

if [[ $coassembly == "yes" ]]; then
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | COASSEMBLY" >> ${workdir}/run_${timestamp}.log
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/coassembly.sh
else
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Reads will not be co-assembled" >> ${workdir}/run_${timestamp}.log
fi

#########
# Predict genes
#########

if [[ $geneprediction == "yes" ]]; then
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | GENE PREDICTION" >> ${workdir}/run_${timestamp}.log
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/geneprediction.sh
else
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Gene prediction will not be performed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Map reads back to the genes and generate Coverage and Hit tables
#########
if [[ $genemapping == "yes" ]]; then
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | GENE MAPPING" >> ${workdir}/run_${timestamp}.log
export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
sh ${metafunkdirectory}/scripts/genemapping.sh
else
echo "$now | Gene mapping will not be performed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Normalise Coverage and Hit tables
#########
if [[ $tss == "yes" || $css == "yes" ]]; then
  export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/normalisetables.sh
  else
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | Hit and coverage tables will not be normalised" >> ${workdir}/run_${timestamp}.log
fi

#########
# Perform functional annotation
#########
now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $functional == "yes" ]]; then
  if [[ $kegg == "yes" ]]; then
    echo "$now | Starting KEGG functional annotation" >> ${workdir}/run_${timestamp}.log
    export metafunkdirectory; export timestamp; export settingsfile
    sh ${metafunkdirectory}/scripts/functionalannotation_kegg.sh
    else
    echo "$now | KEGG functional annotation will not be performed" >> ${workdir}/run_${timestamp}.log
  fi

  now=$(date +"%Y-%m-%d %H:%M:%S")
  if [[ $eggnog == "yes" ]]; then
    echo "$now | Starting EggNog functional annotation" >> ${workdir}/run_${timestamp}.log
    export metafunkdirectory; export timestamp; export settingsfile
    sh ${metafunkdirectory}/scripts/functionalannotation_eggnog.sh
    else
    echo "$now | EggNog functional annotation will not be performed" >> ${workdir}/run_${timestamp}.log
  fi
fi

#########
# Contig Mapping
#########
if [[ $contigmapping == "yes" ]]; then
  export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/contigmapping.sh
  else
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | Reads will not be mapped to contigs" >> ${workdir}/run_${timestamp}.log
fi

#########
# Estimate number of genomes
#########
now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $genomeestimation == "yes" ]]; then
  echo "$now | Starting genome estimation" >> ${workdir}/run_${timestamp}.log
  export WORKDIR="${workdir}"; export METAFUNKDIR="${metafunkdirectory}"
  Rscript ${metafunkdirectory}/scripts/genomestimation.r
fi

#########
# Perform taxonomic profiling
#########
now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $taxonomic == "yes" ]]; then
  echo "$now | Starting taxonomic profiling" >> ${workdir}/run_${timestamp}.log
  export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/taxonomicannotation.sh
  else
  echo "$now | Taxonomic profiling will not be conducted" >> ${workdir}/run_${timestamp}.log
fi

#########
# Final message
#########
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Congratulations, the pipeline was succesfully run" >> ${workdir}/run_${timestamp}.log

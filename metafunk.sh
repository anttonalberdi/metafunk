#MetaFunk version
version="0.1"

#########
# Get options
#########

while getopts w:d:s:f:t:m:ch option; do
  case "${option}"
  in
  w) workdir=${OPTARG};;
  d) sampledatafile=${OPTARG};;
  s) settingsfile=${OPTARG};;
  f) datadir=${OPTARG};;
  t) threads=${OPTARG};;
  m) modules=${OPTARG};;
  c) compress="TRUE";;
  k) keep="TRUE";;
  h) echo ""
     echo "Metafunk pipeline (alpha version), by Antton Alberdi."
     echo ""
     echo "MANDATORY ARGUMENTS:"
     echo "    -w   Working directory. Absolute path to the directory where the Metafunk pipeline will create and edit files (e.g. /home/antton/wolf_project/)"
     echo "    -d   Sample data file. Absolute path to the file containing information about samples and files (e.g. /home/antton/wolf_project/sample.data.txt)"
     echo "    -s   Settings file Absolute path to the file containing additional settings (e.g. /home/antton/wolf_project/settings.sh)"
     echo "    -f   Raw data directory. Absolute path to the directory where the raw data files are located (e.g. /home/antton/raw_data/wolfs/)"
     echo "    -t   Number of threads [int] (e.g. 8)"
     echo "    -m   Module(s) to run separated with commas. See Metafunk wiki for details: https://github.com/anttonalberdi/metafunk/wiki (e.g. 1,2,3,4)"
     echo ""
     echo "OPTIONAL ARGUMENTS:"
     echo "    -c   Compress files. If used (add "-c" to the command line), source data files will be compressed after accomplishing each step (module)"
     echo "    -k   Keep intermediate files. If used (add "-k" to the command line), intermediate files not necessary for the following steps will NOT be removed"
     echo ""
     exit 0
     ;;
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
if [[ $modules == "15" || $modules =~ ",15," || $modules == 15,* || $modules == *,15 ]];
  then normalisation="yes"
  else normalisation="no"
fi
if [[ $modules == "16" || $modules =~ ",16," || $modules == 16,* || $modules == *,16 ]];
  then functionalcontrasting="yes"
  else functionalcontrasting="no"
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
now=$(date +"%Y-%m-%d %H:%M:%S")

#Check sample name uniqueness
allsamples=$(cut -d' ' -f1 ${sampledatafile} | wc -l)
uniquesamples=$(cut -d' ' -f1 ${sampledatafile} | uniq | wc -l)
if [[ ${allsamples} != ${uniquesamples} ]]; then
  echo "$now | ERROR! Sample names are duplicated" >> ${workdir}/run_${timestamp}.log
  exit
fi

#Check sample file uniqueness
alldata=$(cut -d' ' -f2 ${sampledatafile} | wc -l)
uniquedata=$(cut -d' ' -f2 ${sampledatafile} | uniq | wc -l)
if [[ ${alldata} != ${uniquedata} ]]; then
  echo "$now | ERROR! Sample files are duplicated" >> ${workdir}/run_${timestamp}.log
  exit
fi

#Check genome files


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
  samplenumber=$(cat ${sampledatafile} | wc -l)
  echo "$now | DATA TRANSFER" >> ${workdir}/run_${timestamp}.log
  echo "$now | Copying and uncompressing data files of $samplenumber samples" >> ${workdir}/run_${timestamp}.log
  #Load necessary modules and check they work
  module load ${soft_pigz}
  module load ${soft_parallel}
  dependencylist="pigz,parallel"
  export workdir; export dependencylist; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/checkdependencies.sh
  #Launch script
  export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/transferdata.sh
  #Unload necessary modules
  module unload ${soft_pigz}
  module unload ${soft_parallel}
else
  echo "$now | Data will not be copied and uncompressed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Quality filtering
#########

now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $qualityfiltering == "yes" ]]; then
  echo "$now | QUALITY FILTERING" >> ${workdir}/run_${timestamp}.log
  #Load necessary modules
  module load ${soft_pigz}
  module load ${soft_adapterremoval}
  dependencylist="pigz,adapterremoval"
  export workdir; export dependencylist; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/checkdependencies.sh
  #Launch script
  export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp; export compress
  sh ${metafunkdirectory}/scripts/qualityfiltering.sh
  #Unload necessary modules
  module unload ${soft_pigz}
  module unload ${soft_adapterremoval}
else
  echo "$now | Quality filtering will not be performed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Remove duplicates
#########

now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $removeduplicates == "yes" ]]; then
  echo "$now | DUPLICATE REMOVAL" >> ${workdir}/run_${timestamp}.log
  #Load necessary modules
  module load ${soft_pigz}
  module load ${soft_parallel}
  module load ${soft_seqkit}
  dependencylist="pigz,parallel,seqkit"
  export workdir; export dependencylist; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/checkdependencies.sh
  #Launch script
  export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp; export compress
  sh ${metafunkdirectory}/scripts/removeduplicates.sh
  #Unload necessary modules
  module unload ${soft_pigz}
  module unload ${soft_parallel}
  module unload ${soft_seqkit}
else
  echo "$now | Duplicates will not be removed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Remove low complexity reads
#########

now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $removelowcomplexity == "yes" ]]; then
  echo "$now | LOW COMPLEXITY FILTERING" >> ${workdir}/run_${timestamp}.log
  #Load necessary modules
  module load ${soft_pigz}
  module load ${soft_parallel}
  module load ${soft_perl}
  module load ${soft_prinseq}
  dependencylist="pigz,parallel,perl,prinseq"
  export workdir; export dependencylist; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/checkdependencies.sh
  #Launch script
  export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp; export compress
  sh ${metafunkdirectory}/scripts/lowcomplexity.sh
  #Unload necessary modules
  module unload ${soft_pigz}
  module unload ${soft_parallel}
  module unload ${soft_perl}
  module unload ${soft_prinseq}
else
  echo "$now | Low complexity will not be filtered" >> ${workdir}/run_${timestamp}.log
fi

#########
# Remove host DNA
#########

now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $removehostdna == "yes" ]]; then
  echo "$now | HOST DNA REMOVAL" >> ${workdir}/run_${timestamp}.log
  #Load necessary modules
  module load ${soft_pigz}
  module load ${soft_parallel}
  module load ${soft_openssl}
  module load ${soft_samtools}
  module load ${soft_bwa}
  module load ${soft_jre}
  module load ${soft_bbmap}
  dependencylist="pigz,parallel,openssl,samtools,bwa,jre,bbmap"
  export workdir; export dependencylist; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/checkdependencies.sh
  #Launch script
  export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp; export compress
  sh ${metafunkdirectory}/scripts/removehostdna.sh
  #Unload necessary modules
  module unload ${soft_pigz}
  module unload ${soft_parallel}
  module unload ${soft_openssl}
  module unload ${soft_samtools}
  module unload ${soft_bwa}
  module unload ${soft_jre}
  module unload ${soft_bbmap}
else
  echo "$now | Host DNA will not be removed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Remove human DNA
#########

if [[ $removehumandna == "yes" ]]; then
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | HUMAN DNA REMOVAL" >> ${workdir}/run_${timestamp}.log
  #Load necessary modules
  module load ${soft_pigz}
  module load ${soft_parallel}
  module load ${soft_openssl}
  module load ${soft_samtools}
  module load ${soft_bwa}
  module load ${soft_jre}
  module load ${soft_bbmap}
  dependencylist="pigz,parallel,openssl,samtools,bwa,jre,bbmap"
  export workdir; export dependencylist; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/checkdependencies.sh
  #Launch script
  export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp; export compress
  sh ${metafunkdirectory}/scripts/removehumandna.sh
  #Unload necessary modules
  module unload ${soft_pigz}
  module unload ${soft_parallel}
  module unload ${soft_openssl}
  module unload ${soft_samtools}
  module unload ${soft_bwa}
  module unload ${soft_jre}
  module unload ${soft_bbmap}
else
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | Human DNA will not be removed" >> ${workdir}/run_${timestamp}.log
fi

#########
# Perform co-assembly
#########

now=$(date +"%Y-%m-%d %H:%M:%S")
if [[ $coassembly == "yes" ]]; then
  echo "$now | COASSEMBLY" >> ${workdir}/run_${timestamp}.log
  #Load necessary modules
  module load ${soft_pigz}
  module load ${soft_parallel}
  module load ${soft_fastx}
  module load ${soft_megahit}
  module load ${soft_openssl}
  module load ${soft_samtools}
  module load ${soft_bwa}
  dependencylist="pigz,parallel,fastx,megahit,openssl,samtools,bwa"
  export workdir; export dependencylist; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/checkdependencies.sh
  #Launch script
  export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp; export compress
  sh ${metafunkdirectory}/scripts/coassembly.sh
  #Unload necessary modules
  module unload ${soft_pigz}
  module unload ${soft_parallel}
  module unload ${soft_fastx}
  module unload ${soft_megahit}
  module unload ${soft_openssl}
  module unload ${soft_samtools}
  module unload ${soft_bwa}
else
  echo "$now | Reads will not be co-assembled" >> ${workdir}/run_${timestamp}.log
fi

#########
# Predict genes
#########

now=$(date +"%Y-%m-%d %H:%M:%S")
  if [[ $geneprediction == "yes" ]]; then
  echo "$now | GENE PREDICTION" >> ${workdir}/run_${timestamp}.log
  #Load necessary modules
  module load ${soft_pigz}
  module load ${soft_parallel}
  module load ${soft_prodigal}
  dependencylist="pigz,parallel,prodigal"
  export workdir; export dependencylist; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/checkdependencies.sh
  #Launch script
  export workdir; export sampledatafile; export settingsfile; export datadir; export threads; export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/geneprediction.sh
  #Unload necessary modules
  module unload ${soft_pigz}
  module unload ${soft_parallel}
  module unload ${soft_prodigal}
else
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
# Normalise Coverage and Hit tables
#########
if [[ $normalisation == "yes" ]]; then
  export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/normalisation.sh
  else
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | Hit and coverage tables will not be normalised" >> ${workdir}/run_${timestamp}.log
fi

#########
# Functional contrasting
#########
if [[ $functionalcontrasting == "yes" ]]; then
  export metafunkdirectory; export timestamp
  sh ${metafunkdirectory}/scripts/functionalcontrasting.sh
  else
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | Functional contrasting will not be conducted" >> ${workdir}/run_${timestamp}.log
fi

#########
# Final message
#########
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Congratulations, the pipeline was succesfully run" >> ${workdir}/run_${timestamp}.log

#Source settings file
source $settingsfile

now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Contrasting functional profiles" >> ${workdir}/run_${timestamp}.log
export WORKDIR="${workdir}"
export SAMPLEDATAFILE="${sampledatafile}"
export NORMALISATIONMETHOD="${normalisationmethod}"
Rscript ${metafunkdirectory}/scripts/functionalcontrasting_KEGG.r --no-save

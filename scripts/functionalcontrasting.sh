#Source settings file
source $settingsfile

now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Contrasting KEGG functional profiles" >> ${workdir}/run_${timestamp}.log
export WORKDIR="${workdir}"
export SAMPLEDATAFILE="${sampledatafile}"
export NORMALISATIONMETHOD="${normalisationmethod}"
export ANNOTATIONTHRESHOLD="${annotationthreshold}"
Rscript ${metafunkdirectory}/scripts/functionalcontrasting_KEGG.r --no-save

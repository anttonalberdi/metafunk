#Source settings file
source $settingsfile

#Perform TSS normalisation
if [[ $tss == "yes" ]]; then
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | TSS-normalising hit and coverage tables" >> ${workdir}/run_${timestamp}.log
  export WORKDIR="${workdir}"
  export NORMALISATIONSCALE="${normalisationscale}"
  export NORMALISATIONDECIMALS="${normalisationdecimals}"
  Rscript ${metafunkdirectory}/scripts/normalisation_tss.r --no-save
fi
#Perform CSS normalisation
if [[ $css == "yes" ]]; then
  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | CSS-normalising hit and coverage tables" >> ${workdir}/run_${timestamp}.log
  export WORKDIR="${workdir}"
  export NORMALISATIONSCALE="${normalisationscale}"
  export NORMALISATIONDECIMALS="${normalisationdecimals}"
  Rscript ${metafunkdirectory}/scripts/normalisation_css.r --no-save
fi

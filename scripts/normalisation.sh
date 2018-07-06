#Source settings file
source $settingsfile

  now=$(date +"%Y-%m-%d %H:%M:%S")
  echo "$now | Normalising hit and coverage tables" >> ${workdir}/run_${timestamp}.log
  export WORKDIR="${workdir}"
  export NORMALISATIONMETHOD="${normalisationmethod}"
  export NORMALISATIONSCALE="${normalisationscale}"
  export NORMALISATIONDECIMALS="${normalisationdecimals}"
  Rscript ${metafunkdirectory}/scripts/normalisation.r --no-save


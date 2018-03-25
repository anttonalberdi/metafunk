#Source dependencies
source "$metafunkdirectory/settings.sh"

#Perform TSS normalisation
if [[ $tss == "yes" ]]; then
  export WORKDIR="${workingdirectory}/${project}"
  export NORMALISATIONSCALE="${normalisationscale}"
  export NORMALISATIONDECIMALS="${normalisationdecimals}"
  Rscript scripts/createhittable.r --no-save

fi

#Perform CSS normalisation
if [[ $css == "yes" ]]; then
  echo "To be written"
fi

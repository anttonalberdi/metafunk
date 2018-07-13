#Source settings file
source $settingsfile

now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Contrasting KEGG functional profiles" >> ${workdir}/run_${timestamp}.log

#Create directory
mkdir -p ${workdir}/FunctionalStats

#Apply annotation threshold
awk -v threshold="${keggthreshold}" '$5 <= threshold {print $0}' ${workdir}/GeneAnnotationKEGG/assembly.genes.KEGG.annotated.txt > ${workdir}/GeneAnnotationKEGG/assembly.genes.KEGG.annotated.${keggthreshold}.txt

#Run analysis in R
export WORKDIR="${workdir}"
export METAFUNKDIR="${metafunkdirectory}"
export SAMPLEDATAFILE="${sampledatafile}"
export NORMALISATIONMETHOD="${normalisationmethod}"
export KEGGTHRESHOLD="${keggthreshold}"
export TIMESTAMP="${timestamp}"
Rscript ${metafunkdirectory}/scripts/functionalcontrasting.r --no-save

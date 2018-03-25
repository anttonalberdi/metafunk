#Source dependencies
source "$metafunkdirectory/settings.sh"

echo "  Checking dependencies:" >> ${workingdirectory}/${project}/run.log

#Check software needed for transferring and uncompressing data
if [[ $copydata == "yes" ]]; then
test=$(command -v pigz)
    if [[ ${#test} > 0 ]]; then
    echo "    Pigz is installed" >> ${workingdirectory}/${project}/run.log
    else
    echo "    ERROR: Pigz is NOT installed" >> ${workingdirectory}/${project}/run.log
    fi
fi

#Check software needed for removing duplicate sequences
if [[ $removeduplicates == "yes" ]]; then
test=$(command -v seqkit)
    if [[ ${#test} > 0 ]]; then
    echo "    Seqkit is installed" >> ${workingdirectory}/${project}/run.log
    else
    echo "    ERROR: Seqkit is NOT installed" >> ${workingdirectory}/${project}/run.log
    fi
fi

#Check software needed for removing low complexity reads
if [[ $removelowcomplexity == "yes" ]]; then
test=$(command -v prinseq-lite.pl)
    if [[ ${#test} > 0 ]]; then
    echo "    Prinseq is installed" >> ${workingdirectory}/${project}/run.log
    else
    echo "    ERROR: Prinseq is NOT installed" >> ${workingdirectory}/${project}/run.log
    fi
fi

#Check software needed for removing host DNA
if [[ $removehostdna == "yes" || $removehumandna == "yes" || $genemapping == "yes" ]]; then
  test=$(command -v bbmap.sh)
      if [[ ${#test} > 0 ]]; then
      echo "    BBMap is installed" >> ${workingdirectory}/${project}/run.log
      else
      echo "    ERROR: BBMap is NOT installed" >> ${workingdirectory}/${project}/run.log
      fi
test=$(command -v samtools)
    if [[ ${#test} > 0 ]]; then
    echo "    Samtools is installed" >> ${workingdirectory}/${project}/run.log
    else
    echo "    ERROR: Samtools is NOT installed" >> ${workingdirectory}/${project}/run.log
    fi
test=$(command -v bwa)
    if [[ ${#test} > 0 ]]; then
    echo "    Bwa is installed" >> ${workingdirectory}/${project}/run.log
    else
    echo "    ERROR: Bwa is NOT installed" >> ${workingdirectory}/${project}/run.log
    fi
fi

#Check software needed for performing co-assembly
if [[ $coassembly == "yes" ]]; then
  if [[ $seqtype == "PE" ]]; then
    test=$(command -v seqtk)
        if [[ ${#test} > 0 ]]; then
        echo "    Seqtk is installed" >> ${workingdirectory}/${project}/run.log
        else
        echo "    ERROR: Seqtk is NOT installed" >> ${workingdirectory}/${project}/run.log
        fi
  fi
test=$(command -v fastx_toolkit)
    if [[ ${#test} > 0 ]]; then
    echo "    Fastx-toolkit is installed" >> ${workingdirectory}/${project}/run.log
    else
    echo "    ERROR: Fastx-toolkit is NOT installed" >> ${workingdirectory}/${project}/run.log
    fi
test=$(command -v megahit)
    if [[ ${#test} > 0 ]]; then
    echo "    Megahit is installed" >> ${workingdirectory}/${project}/run.log
    else
    echo "    ERROR: Megahit is NOT installed" >> ${workingdirectory}/${project}/run.log
    fi

fi

#Check software needed for gene prediction
if [[ $geneprediction == "yes" ]]; then
test=$(command -v prodigal)
    if [[ ${#test} > 0 ]]; then
    echo "    Prodigal is installed" >> ${workingdirectory}/${project}/run.log
    else
    echo "    ERROR: Prodigal is NOT installed" >> ${workingdirectory}/${project}/run.log
    fi
fi

#Check software needed for gene mapping
if [[ $genemapping == "yes" ]]; then
    test=$(command -v bedtools)
    if [[ ${#test} > 0 ]]; then
    echo "    Bedtools is installed" >> ${workingdirectory}/${project}/run.log
    else
    echo "    ERROR: Bedtools is NOT installed" >> ${workingdirectory}/${project}/run.log
    fi

    test=$(command -v anaconda)
    if [[ ${#test} > 0 ]]; then
    echo "    Anaconda is installed" >> ${workingdirectory}/${project}/run.log
    else
    echo "    ERROR: Anaconda is NOT installed" >> ${workingdirectory}/${project}/run.log
    fi
fi

#Check software needed for normalisation
if [[ $normalisation == "yes" ]]; then
test=$(command -v R)
    if [[ ${#test} > 0 ]]; then
    echo "    R is installed" >> ${workingdirectory}/${project}/run.log
    else
    echo "    ERROR: R is NOT installed" >> ${workingdirectory}/${project}/run.log
    fi
fi

#Check software needed for functional annotation
if [[ $kegg == "yes" || $eggnog == "yes" ]]; then
test=$(command -v diamond)
    if [[ ${#test} > 0 ]]; then
    echo "    Diamond is installed" >> ${workingdirectory}/${project}/run.log
    else
    echo "    ERROR: Diamond is NOT installed" >> ${workingdirectory}/${project}/run.log
    fi
fi

echo "" >> ${workingdirectory}/${project}/run.log

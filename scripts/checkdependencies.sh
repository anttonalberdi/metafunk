#Source dependencies
source "${settingsfile}"

echo "##### DEPENDENCIES ####" >> ${workdir}/run_${timestamp}.log

#Check software needed for transferring and uncompressing data
if [[ $copydata == "yes" ]]; then
test=$(command -v pigz)
    if [[ ${#test} > 0 ]]; then
    echo "  Pigz is installed" >> ${workdir}/run_${timestamp}.log
    else
    echo "  ERROR: The required software Pigz is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

#Check software needed for removing duplicate sequences
if [[ $removeduplicates == "yes" ]]; then
test=$(command -v seqkit)
    if [[ ${#test} > 0 ]]; then
    echo "  Seqkit is installed" >> ${workdir}/run_${timestamp}.log
    else
    echo "  ERROR: The required software Seqkit is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

#Check software needed for removing low complexity reads
if [[ $removelowcomplexity == "yes" ]]; then
test=$(command -v prinseq-lite.pl)
    if [[ ${#test} > 0 ]]; then
    echo "  Prinseq is installed" >> ${workdir}/run_${timestamp}.log
    else
    echo "  ERROR: The required software Prinseq is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

#Check software needed for removing host DNA
if [[ $removehostdna == "yes" || $removehumandna == "yes" || $genemapping == "yes" ]]; then
  test=$(command -v bbmap.sh)
      if [[ ${#test} > 0 ]]; then
      echo "  BBMap is installed" >> ${workdir}/run_${timestamp}.log
      else
      echo "  ERROR: The required software BBMap is NOT installed" >> ${workdir}/run_${timestamp}.log
      fi
test=$(command -v samtools)
    if [[ ${#test} > 0 ]]; then
    echo "  Samtools is installed" >> ${workdir}/run_${timestamp}.log
    else
    echo "  ERROR: The required software Samtools is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
test=$(command -v bwa)
    if [[ ${#test} > 0 ]]; then
    echo "  Bwa is installed" >> ${workdir}/run_${timestamp}.log
    else
    echo "  ERROR: The required software Bwa is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

#Check software needed for performing co-assembly
if [[ $coassembly == "yes" ]]; then
  if [[ $seqtype == "PE" ]]; then
    test=$(command -v seqtk)
        if [[ ${#test} > 0 ]]; then
        echo "  Seqtk is installed" >> ${workdir}/run_${timestamp}.log
        else
        echo "  ERROR: The required software Seqtk is NOT installed" >> ${workdir}/run_${timestamp}.log
        exit
        fi
  fi
test=$(command -v fastq_to_fasta)
    if [[ ${#test} > 0 ]]; then
    echo "  Fastx-toolkit is installed" >> ${workdir}/run_${timestamp}.log
    else
    echo "  ERROR: The required software Fastx-toolkit is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
test=$(command -v megahit)
    if [[ ${#test} > 0 ]]; then
    echo "  Megahit is installed" >> ${workdir}/run_${timestamp}.log
    else
    echo "  ERROR: The required software Megahit is NOT installed" >> ${workdir}/run_${timestamp}.log
    fi

fi

#Check software needed for gene prediction
if [[ $geneprediction == "yes" ]]; then
test=$(command -v prodigal)
    if [[ ${#test} > 0 ]]; then
    echo "  Prodigal is installed" >> ${workdir}/run_${timestamp}.log
    else
    echo "  ERROR: The required software Prodigal is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

#Check software needed for gene mapping
if [[ $genemapping == "yes" ]]; then
    test=$(command -v bedtools)
    if [[ ${#test} > 0 ]]; then
    echo "  Bedtools is installed" >> ${workdir}/run_${timestamp}.log
    else
    echo "  ERROR: The required software Bedtools is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi

    test=$(command -v anaconda)
    if [[ ${#test} > 0 ]]; then
    echo "  Anaconda is installed" >> ${workdir}/run_${timestamp}.log
    else
    echo "  ERROR: The required software Anaconda is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

#Check software needed for normalisation
if [[ $normalisation == "yes" ]]; then
test=$(command -v R)
    if [[ ${#test} > 0 ]]; then
    echo "  R is installed" >> ${workdir}/run_${timestamp}.log
    else
    echo "  ERROR: The required software R is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

#Check software needed for functional annotation
if [[ $kegg == "yes" || $eggnog == "yes" ]]; then
test=$(command -v diamond)
    if [[ ${#test} > 0 ]]; then
    echo "  Diamond is installed" >> ${workdir}/run_${timestamp}.log
    else
    echo "  ERROR: The required software Diamond is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

echo "" >> ${workdir}/run_${timestamp}.log

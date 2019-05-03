#Source dependencies
source "${settingsfile}"

echo "##### DEPENDENCIES ####" >> ${workdir}/run_${timestamp}.log

if [[ $dependencylist =~ "pigz" ]]; then
test=$(command -v pigz)
    if [[ ${#test} = 0 ]]; then
    echo "  ERROR: The required software Pigz is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

if [[ $dependencylist =~ "seqkit" ]]; then
test=$(command -v seqkit)
    if [[ ${#test} = 0 ]]; then
    echo "  ERROR: The required software Seqkit is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

if [[ $dependencylist =~ "prinseq" ]]; then
test=$(command -v prinseq-lite.pl)
    if [[ ${#test} = 0 ]]; then
    echo "  ERROR: The required software Prinseq is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

if [[ $dependencylist =~ "bbmap" ]]; then
  test=$(command -v bbmap.sh)
  if [[ ${#test} = 0 ]]; then
  echo "  ERROR: The required software BBMap is NOT installed" >> ${workdir}/run_${timestamp}.log
  stop
  fi
fi

if [[ $dependencylist =~ "samtools" ]]; then
test=$(command -v samtools)
  if [[ ${#test} = 0 ]]; then
  echo "  ERROR: The required software Samtools is NOT installed" >> ${workdir}/run_${timestamp}.log
  exit
  fi
fi

if [[ $dependencylist =~ "bwa" ]]; then
test=$(command -v bwa)
    if [[ ${#test} = 0 ]]; then
    echo "  ERROR: The required software Bwa is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

if [[ $dependencylist =~ "fastx" ]]; then
test=$(command -v fastq_to_fasta)
    if [[ ${#test} = 0 ]]; then
    echo "  ERROR: The required software Fastx-toolkit is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

if [[ $dependencylist =~ "megahit" ]]; then
test=$(command -v megahit)
    if [[ ${#test} = 0 ]]; then
    echo "  ERROR: The required software Megahit is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

if [[ $dependencylist =~ "megahit" ]]; then
test=$(command -v prodigal)
    if [[ ${#test} = 0 ]]; the
    echo "  ERROR: The required software Prodigal is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

if [[ $dependencylist =~ "bedtools" ]]; then
test=$(command -v bedtools)
    if [[ ${#test} = 0 ]]; then
    echo "  ERROR: The required software Bedtools is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

if [[ $dependencylist =~ "anaconda" ]]; then
test=$(command -v anaconda)
    if [[ ${#test} > 0 ]]; then
    echo "  ERROR: The required software Anaconda is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

if [[ $dependencylist =~ "r" ]]; then
test=$(command -v R)
    if [[ ${#test} > 0 ]]; then
    echo "  ERROR: The required software R is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

if [[ $dependencylist =~ "diamond" ]]; then
test=$(command -v diamond)
    if [[ ${#test} > 0 ]]; then
    echo "  ERROR: The required software Diamond is NOT installed" >> ${workdir}/run_${timestamp}.log
    exit
    fi
fi

echo "" >> ${workdir}/run_${timestamp}.log

echo "Settings file loaded"

project="test1"
datadirectory="/home/projects/pr_46704/people/antalb/metafunk_v0.1/testdata"
workingdirectory="/home/projects/pr_46704/people/antalb/metafunk_v0.1/test_20180322"
metafunkdirectory="/home/projects/pr_46704/people/antalb/metafunk_v0.1"

#### Software dependencies
module load prinseq-lite/0.20.4
module load samtools/1.4
module load bwa/0.7.15
module load pigz/2.3.4
module load fastx_toolkit/0.0.14
module load megahit/1.1.1
module load prodigal/2.6.3
module load bedtools/2.26.0
module load anaconda2/4.0.0

threads="16"

#Data type
readlength="80"
seqtype="SR" #either SR or PE

##### Data transfer settings
copydata="yes"

##### Low complexity settings
removelowcomplexity="yes"
dustvalue="7"

##### Host DNA settings
removehostdna="yes"
hostgenome="/home/projects/pr_46704/people/antalb/databases/L.Dalen_14_wolf.scf.noHets.fasta" #absolute path to fasta file
indexedhostgenome="yes"

##### Coassembly settings
coassembly="yes"

##### Gene prediction settings
geneprediction="yes"

##### Gene mapping settings
genemapping="yes"


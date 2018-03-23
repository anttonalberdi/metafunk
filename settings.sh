echo "Settings file successfully loaded"

#######################################
######### Project information #########
#######################################

#You can use the same working directory for multiple projects (e.g. different settings)

project="test1"
datadirectory="/home/projects/pr_46704/people/antalb/metafunk_v0.1/testdata"
workingdirectory="/home/projects/pr_46704/people/antalb/metafunk_v0.1/test_20180323"
metafunkdirectory="/home/projects/pr_46704/people/antalb/metafunk"

#######################################
######## Software dependencies ########
#######################################

module load prinseq-lite/0.20.4
module load samtools/1.4
module load bwa/0.7.15
module load pigz/2.3.4
module load fastx_toolkit/0.0.14
module load megahit/1.1.1
module load prodigal/2.6.3
module load bedtools/2.26.0
module load anaconda2/4.0.0
module load R/3.2.1

threads="16"

#Data type
readlength="80"
seqtype="SR" #either SR or PE

##### Data transfer settings
copydata="no"

##### Low complexity settings
removelowcomplexity="no"
dustvalue="7"

##### Host DNA settings
removehostdna="yes"
hostgenome="/home/projects/pr_46704/people/antalb/databases/L.Dalen_14_wolf.scf.noHets.fasta" #absolute path to fasta file
indexedhostgenome="no"

##### Coassembly settings
coassembly="yes"
overridecoassembly="yes" #override if the coassembly directory exists

##### Gene prediction settings
geneprediction="yes"

##### Gene mapping settings
genemapping="yes"

##### Coverage/Hit table normalisation
tss="yes"
css="yes"

echo "Settings file successfully loaded"

#######################################
######### Project information #########
#######################################

#You can use the same working directory for multiple projects (e.g. different settings)
project="test1"
datadirectory="/home/projects/pr_46704/people/antalb/metafunk_v0.1/testdata"
workingdirectory="/home/projects/pr_46704/people/antalb/metafunk_v0.1/test_20180325"
metafunkdirectory="/home/projects/pr_46704/people/antalb/metafunk"

#Data type - Platform (Illumina, BGI), Read length (e.g. 80, 100, 150...) and read type (SR or PE)
platform="Illumina"
readlength="150"
seqtype="PE"

#Computational information
threads="16"

#######################################
######## Software dependencies ########
#######################################

#Pigz (https://github.com/madler/pigz)
module load pigz/2.3.4
#AdapterRemoval (https://github.com/MikkelSchubert/adapterremoval)
module load adapterremoval/2.2.2
#SeqKit (https://github.com/shenwei356/seqkit)
module load seqkit/0.7.1
#Prinseq (http://prinseq.sourceforge.net/)
module load prinseq-lite/0.20.4
#Samtools (http://www.htslib.org/doc/samtools.html)
module load samtools/1.4
#Bwa (http://bio-bwa.sourceforge.net/)
module load bwa/0.7.15
#SeqTK (https://github.com/lh3/seqtk)
module load jre/1.8.0
module load bbmap/36.49
module load seqtk/1.0-r82-dirty
#FastX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
module load fastx_toolkit/0.0.14
#Megahit (https://github.com/voutcn/megahit)
module load megahit/1.1.1
module load prodigal/2.6.3
module load bedtools/2.26.0
module load anaconda2/4.0.0
module load R/3.2.1
module load diamond/0.9.13

#######################################
########## Module switcher ############
#######################################

#Select (yes / no) which modules you want to be ran.
#Read the wiki to learn more about what each module does: https://github.com/anttonalberdi/metafunk/wiki

copydata="no"
qualityfiltering="no"
removeduplicates="no"
removelowcomplexity="no"
removehostdna="yes"
removehumandna="no"
coassembly="no"
geneprediction="no"
genemapping="no"

#######################################
############## Settings ###############
#######################################

##### Data transfer settings
compressed="yes"

##### Quality filtering settings
minavgquality="30"
minseqlength="30"
qualitymax="60"

##### Remove duplicate sequences

##### Low complexity settings
dustvalue="7"

##### Host DNA removal settings
indexedhostgenome="yes"
repairingregex="^@(\S+) [1|2](\S+)"

##### Human DNA removal settings
humangenome="/home/projects/pr_46704/people/antalb/databases/Homo_sapiens.fasta" #absolute path to fasta file
indexedhumangenome="no"

##### Coassembly settings
overridecoassembly="yes" #override if the coassembly directory exists

##### Gene prediction settings

##### Gene mapping settings

##### Coverage/Hit table normalisation
tss="no"
css="no"
normalisationscale="1000000" # range of values in nucleotide hits (e.g. 0-1, 0-100, 0-10000, etc.)
normalisationdecimals="2" #increase for smaller scales (e.g. 0-1)

##### Functional annotation KEGG
kegg="no"
keggdatabase="/home/projects/pr_46704/people/antalb/databases/KEGG_species_prokaryotes.dmnd" #absolute path to diamond (.dmnd) file
keggevalue="0.0000000001"

##### Functional annotation EggNog
eggnog="no"

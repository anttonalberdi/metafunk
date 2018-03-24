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

#Pigz (https://github.com/madler/pigz)
module load pigz/2.3.4
#AdapterRemoval (https://github.com/MikkelSchubert/adapterremoval)
module load adapterremoval/2.2.2
#SeqKit (https://github.com/shenwei356/seqkit)
module load seqkit/0.7.1
module load prinseq-lite/0.20.4
module load samtools/1.4
module load bwa/0.7.15
module load seqtk/1.0-r82-dirty
module load fastx_toolkit/0.0.14
module load megahit/1.1.1
module load prodigal/2.6.3
module load bedtools/2.26.0
module load anaconda2/4.0.0
module load R/3.2.1
module load diamond/0.9.13

#######################################
############## Settings ###############
#######################################

threads="16"

#Data type
readlength="80"
seqtype="SR" #either SR or PE

##### Data transfer settings
copydata="yes"

##### Quality filtering
qualityfiltering="no"

##### Remove duplicate sequences
removeduplicates="no"

##### Low complexity settings
removelowcomplexity="no"
dustvalue="7"

##### Host DNA removal settings
removehostdna="no"
hostgenome="/home/projects/pr_46704/people/antalb/databases/L.Dalen_14_wolf.scf.noHets.fasta" #absolute path to fasta file
indexedhostgenome="no"

##### Human DNA removal settings
removehumandna="no"
humangenome="/home/projects/pr_46704/people/antalb/databases/Homo_sapiens.fasta" #absolute path to fasta file
indexedhumangenome="no"

##### Coassembly settings
coassembly="no"
overridecoassembly="yes" #override if the coassembly directory exists

##### Gene prediction settings
geneprediction="no"

##### Gene mapping settings
genemapping="no"

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

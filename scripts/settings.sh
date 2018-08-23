#######################################
######### Project information #########
#######################################

#Data type - Platform (Illumina, BGI)
platform="Illumina"

#######################################
######## Software dependencies ########
#######################################
# Software loading commands

#Parallel (https://www.gnu.org/software/parallel/)
module load parallel/20170822
#Pigz (https://github.com/madler/pigz)
module load pigz/2.3.4
#AdapterRemoval (https://github.com/MikkelSchubert/adapterremoval)
module load adapterremoval/2.2.2
#SeqKit (https://github.com/shenwei356/seqkit)
module load seqkit/0.7.1
#Prinseq (http://prinseq.sourceforge.net/)
module load prinseq-lite/0.20.4
#Samtools (http://www.htslib.org/doc/samtools.html)
module load samtools/1.8
#Bwa (http://bio-bwa.sourceforge.net/)
module load bwa/0.7.15
#SeqTK (https://github.com/lh3/seqtk)
module load jre/1.8.0
#BBMap (https://sourceforge.net/projects/bbmap/)
module load bbmap/36.49
#SeqTK (https://github.com/lh3/seqtk)
module load seqtk/1.0-r82-dirty
#FastX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
module load fastx_toolkit/0.0.14
#Megahit (https://github.com/voutcn/megahit)
module load megahit/1.1.1
#Prodigal (https://github.com/hyattpd/Prodigal)
module load prodigal/2.6.3
#BedTools (http://bedtools.readthedocs.io/en/latest/)
module load bedtools/2.26.0
#Anaconda (https://conda.io/docs/)
module load anaconda2/4.0.0
#R (https://www.r-project.org/)
module load R/3.2.1
#Diamond (https://github.com/bbuchfink/diamond)
module load diamond/0.9.13
#Metaphlan
module load metaphlan/2.6.0
#Bowtie
module load bowtie2/2.3.2

#######################################
############## Settings ###############
#######################################

##### Quality filtering settings
minavgquality="30"
minseqlength="30"
qualitymax="60"

##### Low complexity settings
dustvalue="7"

##### Host DNA removal settings
indexhostgenome="yes"

##### Human DNA removal settings
humangenomepath="/home/projects/pr_46704/people/antalb/databases/Homo_sapiens.fasta" #absolute path to fasta file
indexhumangenome="no"

##### Coassembly settings
overridecoassembly="yes" #override if the coassembly directory exists
indexassembly="no"

##### Gene prediction settings

##### Gene mapping settings

##### Coverage/Hit table normalisation
tss="yes"
css="no"
normalisationscale="1000000" # range of values in nucleotide hits (e.g. 0-1, 0-100, 0-1000000, etc.)
normalisationdecimals="2" #increase for smaller scales (e.g. 0-1)

##### Functional annotation KEGG
kegg="yes"
keggdatabase="/home/projects/pr_46704/people/antalb/databases/KEGG_species_prokaryotes_metafunk.pep" #absolute path to fasta (.fasta, .fa, .pep) or diamond (.dmnd) file
genes_ko="/home/projects/pr_46704/people/antalb/databases/keggs_database/genes_ko.list"
genes_path="/home/projects/pr_46704/people/antalb/databases/keggs_database/genes_pathway.list"
keggevalue="0.0000000001"

##### Functional annotation EggNog
eggnog="no"
eggnogdatabase="/home/projects/pr_46704/people/antalb/databases/eggnog/bactNOG.raw_algs.tar.gz"
eggnogannotations="/home/projects/pr_46704/people/antalb/databases/eggnog/NOG.annotations.tsv"
eggnogmembers="/home/projects/pr_46704/people/antalb/databases/eggnog/NOG.members.tsv"

##### Binning
cutcontigs="yes"
chunksize="10000"
mincontigsize="1000"

##### Source load message
echo "Settings file successfully loaded"

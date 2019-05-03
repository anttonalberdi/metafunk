#######################################
######## Software dependencies ########
#######################################

#Parallel (https://www.gnu.org/software/parallel/)
soft_parallel=parallel/20170822
#Pigz (https://github.com/madler/pigz)
soft_pigz=pigz/2.3.4
#AdapterRemoval (https://github.com/MikkelSchubert/adapterremoval)
soft_adapterremoval=adapterremoval/2.2.2
#SeqKit (https://github.com/shenwei356/seqkit)
soft_seqkit=seqkit/0.7.1
#Perl
soft_perl=perl/5.24.0
#libgtextutils
soft_libgtext=libgtextutils/0.7
#Prinseq (http://prinseq.sourceforge.net/)
soft_prinseq=prinseq-lite/0.20.4
#OpenSSL
soft_openssl=openssl/1.0.0
#Samtools (http://www.htslib.org/doc/samtools.html)
soft_samtools=samtools/1.4
#Bwa (http://bio-bwa.sourceforge.net/)
soft_bwa=bwa/0.7.15
#SeqTK (https://github.com/lh3/seqtk)
soft_jre=jre/1.8.0
#BBMap (https://sourceforge.net/projects/bbmap/)
soft_bbmap=bbmap/36.49
#SeqTK (https://github.com/lh3/seqtk)
soft_seqtk=seqtk/1.0-r82-dirty
#FastX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
soft_fastx=fastx_toolkit/0.0.14
#Megahit (https://github.com/voutcn/megahit)
soft_megahit=megahit/1.1.1
#Prodigal (https://github.com/hyattpd/Prodigal)
soft_prodigal=prodigal/2.6.3
#BedTools (http://bedtools.readthedocs.io/en/latest/)
soft_bedtools=bedtools/2.26.0
#Anaconda (https://conda.io/docs/)
soft_anaconda=anaconda2/4.0.0
#R (https://www.r-project.org/)
soft_r=R/3.2.1
#Diamond (https://github.com/bbuchfink/diamond)
soft_diamond=diamond/0.9.13

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
humangenomepath="/home/projects/ku-cbd/people/antalb/cat_ext/GCF_000001405.38_GRCh38.p12_genomic.fna" #absolute path to fasta file
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

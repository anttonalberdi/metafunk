# MetaFunk (alpha)
Pipeline for functional metagenomics profiling. Beware this is still an alpha developmental version.

MetaFunk is suitable for Illumina and BGISeq single-read (SR) and paired-end (PE) sequencing data, and by March 2018, it incorporates the following steps:

1. Data transference
2. Quality filtering
3. Duplicate removal
4. Low complexity read removal
5. Host DNA removal (for non-human samples)
6. Human DNA removal (contamination for non-human samples)
7. Co-assembly
8. Gene prediction
9. Gene mapping
10. Functional profiling (KEGG and EggNog databases)
11. Taxonomic profiling

## Quick start
1. Clone or download metafunk to your working environment

`cd /home/software/`

`git clone https://github.com/anttonalberdi/metafunk.git`

2. Edit the sample data file (sample.data.txt) and save it in another directory*

`nano /home/software/metafunk/sample.data.txt`

`mv /home/software/metafunk/sample.data.txt /home/projects/wolfproject/sample.data.txt`

3. Modify the settings.sh file and save it in another directory*

`nano /home/software/metafunk/settings.sh`

`mv /home/software/metafunk/settings.sh /home/projects/wolfproject/settins.sh`

4. Run the following script:

`sh /[path]/metafunk.sh -w [working directory] -d [/path/sample.data.txt] -s [/path/settings.sh] -f [data file directory] -t [number of threads] -m [module list]`

example:

`sh /home/software/metafunk/metafunk.sh -w /home/projects/wolfproject -d /home/projects/wolfproject/sample.data.txt -s /home/projects/wolfproject/settings.sh -f /home/rawdata/wolfdata  -t 24 -m 1,2,3,4,5,6,7,8,9`


*read wiki for further details

## Change log
2018/03/30 | The pipeline is ready until the co-assembly step for SR and PE, SF and MF datasets.

2018/03/22 | The adventure has started. This is still a protoalpha version of the metafunk pipeline.

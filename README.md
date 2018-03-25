# MetaFunk (alpha)
Pipeline for functional metagenomics profiling. Beware this is still an alpha developmental version.

MetaFunk is suitable for Illumina and BGISeq single-read (SR) and paired-end (PE) sequencing data, and by March 2018, it incorporates the following steps:

- Quality filtering
- Duplicate removal
- Low complexity read removal
- Host DNA removal (for non-human samples)
- Human DNA removal (contamination for non-human samples)
- Co-assembly
- Gene prediction
- Gene mapping
- Functional profiling (KEGG and EggNog databases)
- Taxonomic profiling

----

To make it work, clone or download metafunk to your working environment, modify the settings.sh file and run the following script:

cd [metafunkdirectory]
sh metafunk.sh pwd

2018/03/22 | The adventure has started. This is still a protoalpha version of the metafunk pipeline.

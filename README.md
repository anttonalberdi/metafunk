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

## Quick start
1. Clone or download metafunk to your working environment
2. Prepare the sample data file (sample.data.txt)*
3. Modify the settings.sh file*
4. Run the following script:

`metafunkdir="/[absolutepath]/metafunk"

sh $metafunkdir/metafunk.sh $metafunkdir`

*read wiki for further details

## Change log
2018/03/30 | The pipeline is ready until the co-assembly step for SR and PE, SF and MF datasets.
2018/03/22 | The adventure has started. This is still a protoalpha version of the metafunk pipeline.

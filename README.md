# Wagner et al 2020

This repository contains the workflow for analyzing the RNA-seq data in 

* Wagner _et al_ (2020): Mass Cytometric and Transcriptomic Profiling of Epithelial-Mesenchymal Transitions in Human Mammary Cell Lines

The workflow structure is based on [`ARMOR`](https://github.com/csoneson/armor) ([Orjuela, Huang, Hembach _et al_, 2019](https://www.g3journal.org/content/9/7/2089.long)). To run the workflow and regenerate the results, follow the instructions below.

### Preparation

#### Download the data

The FASTQ files have been uploaded to ArrayExpress, with accession number E-MTAB-9365. Download the 32 FASTQ files, and place them in a subfolder named `FASTQ`.

#### Download the reference files

The analysis was performed using the Gencode v34 reference. Download the following files, unzip them, and place them in a subfolder named `reference_files`:

* [genome](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz)
* [transcriptome](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz)
* [gtf file](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz)

#### Run the workflow

To run the workflow, first set up the conda environments:

```
snakemake --use-conda setup --cores 1
```

Check that all inputs are correctly specified:

```
snakemake --use-conda checkinputs --cores 1
```

Then you can run the workflow as follows:

```
snakemake --use-conda --cores 16
```

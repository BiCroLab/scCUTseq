# Workflow for (sc)CUTseq preprocessing and copy number calling

## Preprocessing
### Dependencies
Python:
* argparse
* gzip
* re
* os
* multiprocessing
* functools
* collections
* itertools
* operator
* Biopython
* distance
* glob
* pandas

Other software:
* snakemake
* alfred
* umi_tools

### Preprocessing workflow
This workflow expects standard illumina demultiplexed `fastq` files as input (can be either multiple `fastq` files for each lane or merged).
To process the a single multiplex library edit the `config.yaml` file in the `preprocessing` folder with the appropriate values.

A reference genome, the link to the `alfred` binary, a link the the `processBarcode.py` and a list of barcodes are also required. The list of barcodes should be a standard text file with each row containing one barcode. 

The `barcode pattern` depends on the sequencing method used. For standard CUTseq (96 linkers) this is typically `UUUUUUUUBBBBBBBBDDDD` and for scCUTseq (384 linkers) this is typically `UUUUUUUUBBBBBBBBBBBDDDD`.

After successfully editing the `config.yaml` you can run the entire snakemake file from the command line. If you are in the `preprocessing` directory this can be done in the following way: `snakemake -j {number of threads to use}`

## CNA calling
### Dependencies
Python:
* snakemake
* glob
* os

R:
* data.table
* argparser
* DNAcopy
* ParDNAcopy
* pbapply
* gtools
* RColorBrewer
* scales
* ggdendro
* cowplot
* ggplot2
* tidyverse

### CNA calling workflow
After finishing the preprocessing workflow you should have aligned bamfiles. One for each of the barcodes provided in the barcode file.  
Just like for the preprocessing workflow, here you also just have to change the `config.yaml` and then execute the snakemake. There are, just like in the preprocessing workflow, some file requirements. The main files you need is a blacklist (provided in this repository under `/files/hg19/hg19-blacklist.v2_adjusted.bed` and a couple of files regarding the variable widths of the bins. Opposed to fixed-width binning, the variable-width bins take into account the mappability of the genome (with the length of the reads that you use). Within the directory `/files/hg19/` there are included some of the main files. Specifically, there are files for different binsizes, `1mb`, `500kb`, `250kb`, `175kb`, `100kb` and `50kb`. Also, for each of these binsizes there are files with variable widths generated for the following readlengths: `48`, `76` and `150`.  
If you want to generate your own files I suggest you check out the following script: https://github.com/EngeLab/DNTRseq/blob/master/scripts/build_genome_files_bybin.sh

If you have all the required file and specified all the paths and parameters correctly, you can run snakemake like before `snakemake -j {number of threads to use}`. Note: you can specify multiple binsizes, i.e. `binsize: [100000, 50000]`, and the analysis will run for all the specified binsizes.

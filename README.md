![PRONAME_logo](./images/PRONAME_logo.jpg?raw=true "PRONAME logo")

# PRONAME: PROcessing NAnopore MEtabarcoding data

PRONAME is an open-source bioinformatics pipeline that allows processing Nanopore metabarcoding sequencing data. The pipeline is written mainly in bash and is precompiled in a conda environment package which simply needs to be decompressed and activated to be ready to use. The PRONAME package includes all developed scripts, dependencies and precompiled reference databases.

The pipeline is divided into four steps: (i) Nanopore sequencing data is first imported into PRONAME to trim adapter and primer sequences (optional) and to visualize raw read length and quality (`proname_import`). (ii) One of the main advantages of the second script of the pipeline (`proname_filter`) is that it allows diffentiating simplex from duplex reads and, thus, take advantage of higher-accuracy duplex reads introduced with the V14 sequencing chemistry. Reads that do not meet length and quality criteria are then filtered out. (iii) The next script of the pipeline (`proname_refine`) performs a read clustering, uses Medaka, i.e. a Nanopore data-dedicated tool, to correct sequencing errors by polishing, and discards chimera sequences. (iv) The last script (`proname_taxonomy`) allows performing the taxonomic analysis of the generated high-accuracy consensus sequences. The pipeline offers the possibility to import the generated files into QIIME2 for further analyses (diversity, abundance, etc.), if desired.

# Installation

If you don't have [miniconda3](https://docs.anaconda.com/free/miniconda/) installed, please execute the following commands to install it:

~~~
# Installing miniconda3
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~~~


The PRONAME environment package is available here: FIGSHARE LINK

You can then run these commands to install and activate the PRONAME environment:

~~~
# Creating the environment directory
mkdir ~/miniconda3/envs/proname

# Extracting the environment
tar -xzf proname.tar.gz -C ~/miniconda3/envs/proname

# Activating the environment
source ~/miniconda3/envs/proname/bin/activate

# Removing prefixes
conda-unpack
~~~


And that's it! You are now ready to analyze your nanopore metabarcoding data with PRONAME. There are four scripts constituting the pipeline:

* proname_import
* proname_filter
* proname_refine
* proname_taxonomy

These scripts must be run in this order, with their required arguments. 
The best way to go is to type the name of each script followed by "--help" (e.g. `proname_import --help`) to get the list of all arguments and a usage example. A tutorial detailing the whole workflow is also presented below.

# Tutorial

## 0. Before PRONAME

The Nanopore sequencing data to import into PRONAME must be fastq files, i.e. basecalled reads. If you have raw-signal data (fast5 or pod5 files), you should first basecall them preferably with [Dorado](https://github.com/nanoporetech/dorado). For this tutorial, all fastq files (one file per sample) have been placed in the `RawData` directory.

## 1. proname_import

The first step is to import sequencing data into PRONAME. Since adapter and primer sequences have not been removed yet, the `trimadapters` and `trimprimers` arguments are set to "yes" and the primer sequences are provided (5'-3'). The V14 sequencing chemistry was used, so the `duplex` argument is set to "yes", so that seperate length-vs-quality scatterplots are generated for simplex, duplex and simplex+duplex reads.

~~~
proname_import \
  --inputpath RawData \
  --threads 48 \
  --duplex yes \
  --trimadapters yes \
  --trimprimers yes \
  --fwdprimer AGRGTTYGATYMTGGCTCAG \
  --revprimer CGACATCGAGGTGCCAAAC
~~~

Here is the complete list of available arguments for `proname_import`:

| Command | Arguments | Description |
| ------- | --------- | ----------- |
| proname_import | --inputpath | Path to the folder containing raw fastq files. |
|  | --threads | Number of threads to use for the Guppy adapter-trimming step and/or the Cutadapt primmer-trimming step. |
|  | --duplex | Indicate whether your sequencing data include duplex reads or not. Duplex reads are high-quality reads that were introduced with the kit 14 chemistry. [Option: "yes" or "no"] |
|  | --trimadapters | Indicate whether your sequencing data contain adapters that should be trimmed. [Option: "yes" or "no"] |
|  | --trimprimers | Indicate whether your sequencing data contain primers that should be trimmed. [Option: "yes" or "no"] |
|  | --fwdprimer | The sequence of the forward primer used during PCR to amplify DNA. If barcoded primers were used to multiplex samples, please provide here only the target-specific part of the primer in 5'->3' orientation. This argument is required if --trimprimers is set to "yes". |
|  | --revprimer | The sequence of the reverse primer used during PCR to amplify DNA. If barcoded primers were used to multiplex samples, please provide here only the target-specific part of the primer in 5'->3' orientation. This argument is required if --trimprimers is set to "yes". |
|  | --version | Print the version of the pipeline. |
|  | --help | Print the help menu. |






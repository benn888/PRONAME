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

The first step is to import sequencing data into PRONAME. Since adapter and primer sequences have not been removed yet, the `--trimadapters` and `--trimprimers` arguments are set to "yes" and the primer sequences are provided (5'-3'). Given that the V14 sequencing chemistry was used, the `--duplex` argument is set to "yes", so that seperate length-vs-quality scatterplots are generated for simplex, duplex and simplex+duplex reads.

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

| Command | Arguments | Description | Mandatory arguments |
| ------- | --------- | ----------- | ------------------- |
| proname_import | --inputpath | Path to the folder containing raw fastq files. | X |
|  | --threads | Number of threads to use for the Guppy adapter-trimming step and/or the Cutadapt primmer-trimming step. You can know the number of available threads on your computer by running the command 'nproc --all' [Default: 2] |  |
|  | --duplex | Indicate whether your sequencing data include duplex reads or not. Duplex reads are high-quality reads that were introduced with the kit 14 chemistry. [Option: "yes" or "no"] | X |
|  | --trimadapters | Indicate whether your sequencing data contain adapters that should be trimmed. [Option: "yes" or "no"] | X |
|  | --trimprimers | Indicate whether your sequencing data contain primers that should be trimmed. [Option: "yes" or "no"] | X |
|  | --fwdprimer | The sequence of the forward primer used during PCR to amplify DNA. If barcoded primers were used to multiplex samples, please provide here only the target-specific part of the primer in 5'->3' orientation. This argument is required if --trimprimers is set to "yes". | ~ |
|  | --revprimer | The sequence of the reverse primer used during PCR to amplify DNA. If barcoded primers were used to multiplex samples, please provide here only the target-specific part of the primer in 5'->3' orientation. This argument is required if --trimprimers is set to "yes". | ~ |
|  | --version | Print the version of the pipeline. |  |
|  | --help | Print the help menu. |  |

The analysis of the `simplex_duplex_read_distribution.tsv` generated file shows that enough duplex reads were sequenced: 

| Sample_name | Simplex_reads | Duplex_reads |
| ----------- | ------------- | ------------ |
| sample10 | 405356 | 38219 |
| sample1 | 254762 | 58214 |
| sample2 | 473202 | 33534 |
| sample3 | 250737 | 48763 |
| sample4 | 214401 | 14949 |
| sample5 | 260615 | 20368 |
| sample6 | 359033 | 53987 |
| sample7 | 366997 | 47476 |
| sample8 | 476416 | 101848 |
| sample9 | 711157 | 55596 |

We can thus work only with duplex reads and discard simplex reads. Analyzing the `LengthvsQualityScatterPlot_duplex.png` file helps determine which length and quality thresholds to choose at the following step:

![Duplex_scatterplot](./images/LengthvsQualityScatterPlot_duplex.png?raw=true "Duplex_scatterplot")

## 2. proname_filter

Based on the simplex/duplex read distribution and the length-vs-quality scatterplot generated at the previous step, we decide to carry on the analysis only with duplex reads, and to discard those with a Q score lower than 15 and those shorter than 3500 bp or longer than 5000 bp:

~~~
proname_filter \
  --datatype duplex \
  --filtminlen 3500 \
  --filtmaxlen 5000 \
  --filtminqual 15 \
  --threads 48 \
  --inputpath RawData
~~~

Here is the complete list of available arguments for `proname_filter`:

| Command | Arguments | Description | Mandatory arguments |
| ------- | --------- | ----------- | ------------------- |
| proname_filter | --datatype | Indicate whether you want to work with simplex reads, duplex reads or both. [Option: "simplex", "duplex" or "both"] | X |
|  | --filtminlen | Reads with a length below this threshold will be discarded during quality filtering. [Option: integer] | X |
|  | --filtmaxlen | Reads with a length above this threshold will be discarded during quality filtering. [Option: integer] | X |
|  | --filtminqual | Reads with a quality score below this threshold will be discarded during quality filtering. [Option: integer] | X |
|  | --threads | Number of threads to use. You can know the number of available threads on your computer by running the command 'nproc --all' [Default: 2] |  |
|  | --inputpath | Path to the folder containing raw fastq files. This must be the same path than the one provided while running proname_import. | X |
|  | --deletefiles | Delete all non-essential files, i.e. files generated with proname_import that are no more needed for the rest of the analysis through PRONAME. [Option: "yes" or "no", Default: no] |  |
|  | --version | Print the version of the pipeline. |  |
|  | --help | Print the help menu. |  |

The `HQ_duplex_read_distribution.tsv` file indicates how many high-quality duplex reads remained after filtering:

| Sample_name | HQ_duplex_reads |
| ----------- | --------------- |
| sample10 | 21326 |
| sample1 | 49029 |
| sample2 | 24610 |
| sample3 | 41600 |
| sample4 | 11212 |
| sample5 | 17422 |
| sample6 | 38546 |
| sample7 | 32267 |
| sample8 | 73018 |
| sample9 | 31843 |


The new scatterplot `LengthvsQualityScatterPlot_HQ_duplex.png` allows visualizing the impact of these filterings on the read quality and length distributions:

![HQ_Duplex_scatterplot](./images/LengthvsQualityScatterPlot_HQ_duplex.png?raw=true "HQ_Duplex_scatterplot")

## 3. proname_refine

The HQ duplex reads will now undergo a serie of processing steps wrapped in the proname_refine script:

* They are clustered according to a similarity threshold defined by the `--clusterid` argument;
* The centroid sequence is extracted from each cluster;
* 300 other reads, defined as 'sub-reads', are randomly extracted from each cluster;
* Each centroid sequence is polished using its associated sub-reads to generate an error-corrected consensus sequence;
* Chimera sequences are removed;
* The generated files can then be imported into QIIME2 if desired, by setting the `--qiime2import` argument to 'yes'.

~~~
proname_refine \
  --clusterid 0.90 \
  --inputpath RawData \
  --medakamodel r1041_e82_400bps_sup_v4.2.0 \
  --chimeradb /home/utilisateur/miniconda3/envs/proname/db/rEGEN-B/regenB_sequences.fasta \
  --qiime2import yes
~~~

Here is the complete list of available arguments for `proname_refine`:

| Command | Arguments | Description | Mandatory arguments |
| ------- | --------- | ----------- | ------------------- |
| proname_refine | --clusterid | The percentage of identity at which clustering should be performed. [Option: decimal between 0 and 1] | X |
|  | --clusterthreads | Number of threads to use for the clustering step. You can know the number of available threads on your computer by running the command 'nproc --all' [Default: 2] |  |
|  | --inputpath | Path to the folder containing raw fastq files. This must be the same path than the one provided while running proname_import and proname_filter. | X |
|  | --subsampledreads | Number of subsampled reads that will be aligned against the centroid sequence during polishing. [Default: 300] |  |
|  | --medakabatchsize | Controls memory use. Medaka developers set the default value to 100 but it was reduced to 20 in PRONAME because it is one of the main reasons why medaka may crash due to insufficient available memory. Feel free to increase this value if your working machine has enough memory. [Default: 20] |  |
|  | --medakathreads | Number of threads to use for the polishing step. The default value has been set to 1 because the lack of memory is one of the main reasons why medaka may crash. Feel free to increase this value if your working machine has enough memory. You can know the number of available threads on your computer by running the command 'nproc --all' [Default: 1] |  |
|  | --medakamodel | Basecalling model used to generate raw fastq files. This model will be used by medaka to polish data. The list of available models can be found by running 'medaka tools list\_models' | X |
|  | --chimeradb | Path to the reference database to use for the chimera detection. | X |
|  | --qiime2import | Indicate whether the generated representative sequences and table must be imported into QIIME2. [Option: "yes" or "no"] | X |
|  | --deletefiles | Delete all non-essential files, i.e. files generated with proname_filter that are no more needed for the rest of the analysis through PRONAME. [Option: "yes" or "no", Default: no] |  |
|  | --version | Print the version of the pipeline. |  |
|  | --help | Print the help menu. |  |

## 4. proname_taxonomy

The files generated at the previous step gathering all the consensus sequences (`rep-seqs.qza`) and associated frequency table (`rep-table.qza`) are used to perform the taxonomic analysis and produce a taxonomy file and a taxa barplot:

~~~
proname_taxonomy \
  --qseqs rep_seqs.qza \
  --qtable rep_table.qza \
  --db /home/utilisateur/miniconda3/envs/proname/db/rEGEN-B/regenB_sequences.fasta \
  --reftax /home/utilisateur/miniconda3/envs/proname/db/rEGEN-B/regenB_taxonomy.tsv \
  --assay assay1
~~~

Here is the complete list of available arguments for `proname_taxonomy`:

| Command | Arguments | Description | Mandatory arguments |
| ------- | --------- | ----------- | ------------------- |
| proname_taxonomy | --qseqs | Path to query sequences. It should be either 'rep_seqs.qza' (if data was imported into QIIME2 at the previous step) or 'rep_seqs.fasta'. | X |
|  | --qtable | Path to query sequence abundance table. It should be 'rep_table.qza'. This argument is only needed if data was imported into QIIME2 at the previous step. | ~ |
|  | --db | Path to the name of the reference database used by blastn to carry out the taxonomic analysis. The rEGEN-B (rrn operons Extracted from GENomes of Bacteria) database as well as the Silva138 and Greengenes2 databases are already precompiled in the PRONAME environment and are located in the folder '$HOME/miniconda3/envs/proname/db' (Change "proname" in this path if you named your environment differently). The user can provide the name of another blastn database if desired. Note that the database must be formatted to run with the BLAST Command Line Applications (for more info, see: https://www.ncbi.nlm.nih.gov/books/NBK569841/). | X |
|  | --reftax | Path to the taxonomic lineages of the sequences included in the reference database. The taxonomy tsv files associated with the rEGEN-B, Silva138 and Greengenes2 databases are already precompiled in the PRONAME environment and are located in the folder '$HOME/miniconda3/envs/proname/db' (Change "proname" in this path if you named your environment differently). The user can provide another reference database and associated taxonomic lineages if desired. | X |
|  | --evalue | Expectation value (E) threshold to keep hits. [Default: 0.001] |  |
|  | --percid | Percent identity threshold between query and subject sequences to keep hits. [Default: 80] |  |
|  | --qcover | Percent query coverage threshold to keep hits. Note that this option does not correspond to the blastn -qcov_hsp_perc option. Here, --qcover reflects the percent coverage of the whole query sequence, so that the results provided are more consistant with those obtained with the online BLASTn tool. [Default: 80] |  |
|  | --threads | Number of threads to use for the blastn analysis. You can know the number of available threads on your computer by running the command 'nproc --all' [Default: 2] |  |
|  | --metadata | Path to the metadata file. It should be 'sample_metada.tsv'. This argument is only required if data was imported into QIIME2 at the previous step. | ~ |
|  | --assay | Name of your metabarcoding assay, that will appear in the name of taxonomy files produced. | X |
|  | --version | Print the version of the pipeline. |  |
|  | --help | Print the help menu. |  |

The taxa barplot generated can be used to visualize the results of the taxonomic analysis using [QIIME2 View](https://view.qiime2.org/) for instance:



This marks the end of the PRONAME pipeline, which allowed generating high-accuracy consensus sequences starting from raw Nanopore metabarcoding data, and performing the subsequent taxonomic analysis.

## 5. After PRONAME

Since we decided to import the generated files into QIIME2, it is easy to carry on with further analyses using this platform. In this tutorial, we will focus on diversity analysis and differential abundance testing.

First of all, we will add information in the metadata file to define to which sample group each sample belongs:

~~~
awk 'BEGIN {FS=OFS="\t"} \
  $1 ~ /sample-id/ {$3="treatment"} \
  $1 ~ /q2:types/ {$3="categorical"} \
  $1 ~ /sample1|sample2|sample3|sample4|sample5/ {$3="treatment1"} \
  $1 ~ /sample6|sample7|sample8|sample9|sample10/ {$3="treatment2"} \
  1' sample-metadata.tsv > sample-metadata.tsv_tmp && mv sample-metadata.tsv_tmp sample-metadata.tsv
~~~

### 5.1. Diversity analysis

#### 5.1.1. Alpha diversity

We will first compute some alpha diversity metrics:

~~~
qiime diversity alpha \
  --i-table rep_table.qza \
  --p-metric observed_features \
  --o-alpha-diversity alpha_observed_features.qza

qiime diversity alpha \
  --i-table rep_table.qza \
  --p-metric pielou_e \
  --o-alpha-diversity alpha_pielou_evenness.qza

qiime diversity alpha \
  --i-table rep_table.qza \
  --p-metric shannon \
  --o-alpha-diversity alpha_shannon.qza
~~~
Note that it is also possible to carry out alpha rarefaction to select a rarefaction depth before computing diversity metrics, see an example [here](https://docs.qiime2.org/2024.5/tutorials/pd-mice/#alpha-rarefaction-and-selecting-a-rarefaction-depth).

We can now test whether the distribution of sequences is significantly different between both treatment groups:

~~~
qiime diversity alpha-group-significance \
  --i-alpha-diversity alpha_observed_features.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha_observed_features-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity alpha_pielou_evenness.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha_pielou_evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity alpha_shannon.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha_shannon-group-significance.qzv
~~~

#### 5.1.2. Beta diversity

~~~
qiime diversity beta \
  --i-table Taxonomy_table.qza \
  --p-metric braycurtis \
  --o-distance-matrix beta_braycurtis.qza
~~~

~~~
qiime diversity beta-group-significance \
  --i-distance-matrix beta_braycurtis.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column Treatment \
  --o-visualization beta_braycurtis-treatment-significance.qzv \
  --p-pairwise
~~~

### 5.2. Differential abundance

~~~
qiime taxa collapse \
  --i-table Taxonomy_table.qza \
  --i-taxonomy Taxonomy_table_taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table Taxonomy_table_l6.qza

qiime composition ancombc \
  --i-table Taxonomy_table_l6.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-formula 'Treatment' \
  --p-reference-levels 'Treatment::T' \
  --o-differentials ancombc_treatment_l6.qza

qiime composition da-barplot \
  --i-data ancombc_treatment_l6.qza \
  --p-significance-threshold 0.05 \
  --o-visualization da_barplot_treatment_l6.qzv
~~~

  





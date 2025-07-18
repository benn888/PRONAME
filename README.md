![PRONAME_logo](./images/proname_logo.PNG?raw=true "PRONAME logo")

# PRONAME: PROcessing NAnopore MEtabarcoding data

PRONAME is an open-source bioinformatics pipeline that allows processing Nanopore metabarcoding sequencing data. The pipeline is written mainly in bash and is compiled in a Docker image which simply needs to be pulled from Docker Hub to be ready to use. The Docker image includes all developed scripts, dependencies and precompiled reference databases.

The pipeline is divided into four steps: (i) Nanopore sequencing data is first imported into PRONAME to trim adapter and primer sequences (optional) and to visualize raw read length and quality (`proname_import`). (ii) One of the main advantages of the second script of the pipeline (`proname_filter`) is that it allows diffentiating simplex from duplex reads and, thus, take advantage of higher-accuracy duplex reads introduced with the V14 sequencing chemistry. Reads that do not meet length and quality criteria are then filtered out. (iii) The next script of the pipeline (`proname_refine`) performs a read clustering, uses Medaka, i.e. a Nanopore data-dedicated tool, to correct sequencing errors by polishing, and discards chimera sequences. (iv) The last script (`proname_taxonomy`) allows performing the taxonomic analysis of the generated high-accuracy consensus sequences. The pipeline offers the possibility to generate a phyloseq object and to import the generated files into QIIME2 for further analyses (diversity, abundance, etc.), if desired.

The rEGEN-B (rrn operons Extracted from GENomes of Bacteria) database developed in this work is included in the Docker image and is also directly available on [Figshare](https://doi.org/10.6084/m9.figshare.26380702).

Additional information can also be found in our [associated publication](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2024.1483255/full).

# Installation

If you don't have Docker on your computer, you can find instructions to install it [here](https://docs.docker.com/engine/install/).

The PRONAME Docker image is available on [Docker Hub](https://hub.docker.com/repository/docker/benn888/proname/general).

This repository includes two Docker images, each optimized for a specific architecture:

* amd64: For x86_64 (Intel/AMD) processors (i.e. suitable for most Linux & Windows machines)
* arm64: For ARM-based processors (e.g., Apple M1/M2, Raspberry Pi)

These images provide flexibility for users on different hardware platforms.

#### Download
To download an image, please run one of the following commands:

- **Command to pull image for amd64 architecture:**  
  ```bash
   docker pull benn888/proname:v2.1.2-amd64

- **Command to pull image for arm64 architecture:**  
  ```bash
   docker pull benn888/proname:v2.1.2-arm64
Note that, depending on your installation, running Docker commands may require `sudo` privileges.

You can run this command to confirm that the image has successfully been downloaded and is available:

~~~
docker images
~~~

Then, the simplest way to run a new container is to use this command:

~~~
docker run -it --name proname_container benn888/proname:v2.1.2-<arch>
~~~
Where `<arch>` should be replaced by `amd64` or `arm64`.

However, a more effective way to launch a container is to set up a shared volume that mounts a host directory directly in the container. This setup allows access to raw sequencing data in the container and enables direct access to PRONAME results from the host machine:

~~~
docker run -it --rm --name proname_container -v /path/to/host/data:/data benn888/proname:v2.1.2-<arch>
~~~
where `/path/to/host/data` is the path to the directory on your host machine containing the raw sequencing data, and `/data` is the directory in the container where this data will be accessible. Place any files resulting from the PRONAME analysis in `/data` to access them directly from the host machine.

Note that, although we did not encounter any memory issue when testing and using PRONAME, it is good to keep in mind that [fine-tuning Docker's memory usage](https://docs.docker.com/engine/containers/resource_constraints/) may be useful in certain cases.

For Docker Desktop on macOS, in particular, ensure that the setting **Settings > General > Use Rosetta for x86_64/amd64 emulation on Apple Silicon** is unchecked.

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
  --sequencingkit SQK-LSK114 \
  --trimprimers yes \
  --fwdprimer AGRGTTYGATYMTGGCTCAG \
  --revprimer CGACATCGAGGTGCCAAAC
~~~

Here is the complete list of available arguments for `proname_import`:

| Command | Arguments | Description | Mandatory arguments |
| ------- | --------- | ----------- | ------------------- |
| proname_import | --inputpath | Path to the folder containing raw fastq files. To prevent file conflicts and ensure accurate sequence counting, your raw FASTQ files must be stored in a separate directory (e.g., /data/RawData). Do not place them directly in your working directory, as this is where proname_import writes all its output files. Mixing input and output files in the same location can lead to errors and unreliable results. | Yes |
|  | --threads | Number of threads to use for the Guppy adapter-trimming step and/or the Cutadapt primmer-trimming step. You can know the number of available threads on your computer by running the command 'nproc --all' [Default: 2] | No |
|  | --duplex | Indicate whether your sequencing data include duplex reads or not. Duplex reads are high-quality reads that were introduced with the kit 14 chemistry. [Option: "yes" or "no"] | Yes |
|  | --trimadapters | Indicate whether your sequencing data contain adapters that should be trimmed. [Option: "yes" or "no"] | Yes |
|  | --sequencingkit | Name of the ONT sequencing kit used to generate the library(-ies). [Default: "SQK-LSK114"] | No | 
|  | --trimprimers | Indicate whether your sequencing data contain primers that should be trimmed. [Option: "yes" or "no"] | Yes |
|  | --fwdprimer | The sequence of the forward primer used during PCR to amplify DNA. If barcoded primers were used to multiplex samples, please provide here only the target-specific part of the primer in 5'->3' orientation. This argument is required if --trimprimers is set to "yes". | ~ |
|  | --revprimer | The sequence of the reverse primer used during PCR to amplify DNA. If barcoded primers were used to multiplex samples, please provide here only the target-specific part of the primer in 5'->3' orientation. This argument is required if --trimprimers is set to "yes". | ~ |
|  | --nocounting | When this argument is set to "yes", counting of simplex/duplex reads is not performed. [Options: "yes" or "no", Default: "no"] | No |
|  | --plotformat | Format of the scatterplot visualization files produced. It can be either "png" or "html". Since NanoPlot produces empty png plots for an unknown reason, it is only used to generate html visualizations. png plots are produced using the custom script scaleq.py. [Default: "png"] | No |
|  | --noscatterplot | When this argument is set to 'yes', no length vs. quality scatterplot is generated. Since this is a time-consuming step, this possiblity has been made available to increase the flexibility of the pipeline. However, it is strongly discouraged to skip this scatterplot generation. Visual inspection of these plots is crucial for deciding which type of read to work with (duplex and/or simplex) and which length and quality thresholds to apply. [Options: "yes" or "no", Default: "no"] | No |
|  | --version | Print the version of the pipeline. | No |
|  | --verbose | Activate verbose/debug mode (no redirections). | No |
|  | --help | Print the help menu. | No |

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

We can therefore work only with duplex reads and discard simplex reads. Analyzing the `LengthvsQualityScatterPlot_duplex.png` file (or `LengthvsQualityScatterPlot_duplex.html`) helps determine which length and quality thresholds to apply at the following step:

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
| proname_filter | --datatype | Indicate whether you want to work with simplex reads, duplex reads or both. [Option: "simplex", "duplex" or "both"] | Yes |
|  | --filtminlen | Reads with a length below this threshold will be discarded during quality filtering. [Option: integer] | Yes |
|  | --filtmaxlen | Reads with a length above this threshold will be discarded during quality filtering. [Option: integer] | Yes |
|  | --filtminqual | Reads with a quality score below this threshold will be discarded during quality filtering. [Option: integer] | Yes |
|  | --threads | Number of threads to use. You can know the number of available threads on your computer by running the command 'nproc --all' [Default: 2] | No |
|  | --inputpath | Path to the folder containing raw fastq files. This must be the same path than the one provided while running proname_import. To prevent file conflicts and ensure accurate sequence counting, your raw FASTQ files must be stored in a separate directory (e.g., /data/RawData). Do not place them directly in your working directory, as this is where proname_import writes all its output files. Mixing input and output files in the same location can lead to errors and unreliable results. | Yes |
|  | --deletefiles | Delete all non-essential files, i.e. files generated with proname_import that are no more needed for the rest of the analysis through PRONAME. [Option: "yes" or "no", Default: no] | No |
|  | --nocounting | When this argument is set to "yes", counting of HQ simplex/duplex reads is not performed. [Options: "yes" or "no", Default: "no"] | No |
|  | --plotformat | Format of the scatterplot visualization file produced. It can be either "png" or "html". Since NanoPlot produces an empty png plot for an unknown reason, it is only used to generate html visualizations. png plots are produced using the custom script scaleq.py. [Default: "png"] | No |
|  | --noscatterplot | When this argument is set to 'yes', no length vs. quality scatterplot is generated. Since this is a time-consuming step, this possiblity has been made available to increase the flexibility of the pipeline. However, it is strongly discouraged to skip this scatterplot generation. Visual inspection of these plots is crucial for deciding which type of read to work with (duplex and/or simplex) and which length and quality thresholds to apply. [Options: "yes" or "no", Default: "no"] | No |
|  | --version | Print the version of the pipeline. | No |
|  | --verbose | Activate verbose/debug mode (no redirections). | No |
|  | --help | Print the help menu. | No |

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


The new scatterplot `LengthvsQualityScatterPlot_HQ_duplex.png` (or `LengthvsQualityScatterPlot_HQ_duplex.html`) allows visualizing the impact of these filterings on the read quality and length distributions:

![HQ_Duplex_scatterplot](./images/LengthvsQualityScatterPlot_HQ_duplex.png?raw=true "HQ_Duplex_scatterplot")

## 3. proname_refine

The HQ duplex reads will now undergo a serie of processing steps wrapped in the proname_refine script:

* They will be clustered according to a similarity threshold defined by the `--clusterid` argument;
* The centroid sequence will be extracted from each cluster;
* 300 other reads, defined as 'sub-reads', will be randomly extracted from each cluster;
* Each centroid sequence will be polished using its associated sub-reads to generate an error-corrected consensus sequence;
> ⚡ Version 2.0.0 introduced major improvements in polishing performance, with speedups of at least 10× observed on HPC systems.
* Chimera sequences will be removed;
* The generated files can then be imported into QIIME2 if desired, by setting the `--qiime2import` argument to 'yes'.

~~~
proname_refine \
  --clusterid 0.90 \
  --inputpath RawData \
  --medakamodel r1041_e82_400bps_sup_v4.2.0 \
  --chimeradb /opt/db/rEGEN-B/rEGEN-B_sequences.fasta \
  --qiime2import yes
~~~

Here is the complete list of available arguments for `proname_refine`:

| Command | Arguments | Description | Mandatory arguments |
| ------- | --------- | ----------- | ------------------- |
| proname_refine | --clusterid | The percentage of identity at which clustering should be performed. [Option: decimal between 0 and 1] | Yes |
|  | --clusteringmethod | The tool used to cluster sequences. It can be either "vsearch" or "mmseqs2". VSEARCH offers high accuracy and is well-suited for small to medium datasets, but can be slower on large datasets. MMseqs2 is slightly less accurate but dramatically faster, especially on large datasets. [Options: "vsearch" or "mmseqs2", Default: "vsearch"] | No |
|  | --clusterthreads | Number of threads to use for the clustering step. You can know the number of available threads on your computer by running the command 'nproc --all' [Default: 2] | No |
|  | --inputpath | Path to the folder containing raw fastq files. This must be the same path than the one provided while running proname_import and proname_filter. | Yes |
|  | --subsampledreads | Number of subsampled reads that will be aligned against the centroid sequence during polishing. [Default: 300] | No |
|  | --medakabatchsize | Controls memory use. The default value has been set to 100. If Medaka shows out of memory errors, the batch size should be reduced. [Default: 100] | No |
|  | --medakathreads | Number of threads to use for the polishing step. If running Medaka on a CPU, it is recommended to set this value to the maximum number of available threads. PRONAME will automatically split the polishing process into multiple parallel jobs to significantly increase the analysis speed. Each job uses 2 threads for Medaka, with 2 additional threads allocated to each job for parallel overhead. If running on a GPU, this setting has little impact, as the main computation is handled by PyTorch and CUDA. [Default: 4] | No |
|  | --medakamodel | Basecalling model used to generate raw fastq files. This model will be used by medaka to polish data. The list of available models can be found by running 'medaka tools list\_models' | Yes |
|  | --chimeramethod | Specify the chimera detection method: "ref" (reference-based) or "denovo" (de novo detection). [Default: "ref"] | No |
|  | --chimeradb | Path to the reference database to use for the chimera detection. | Yes |
|  | --qiime2import | Indicate whether the generated representative sequences and table must be imported into QIIME2. [Option: "yes" or "no"] | Yes |
|  | --deletefiles | Delete all non-essential files, i.e. files generated with proname_filter that are no more needed for the rest of the analysis through PRONAME. [Option: "yes" or "no", Default: no] | No |
|  | --version | Print the version of the pipeline. | No |
|  | --verbose | Activate verbose/debug mode (no redirections). | No |
|  | --help | Print the help menu. | No |

## 4. proname_taxonomy

The files generated at the previous step gathering all the consensus sequences (`rep-seqs.qza`) and associated frequency table (`rep-table.qza`) are used here to perform the taxonomic analysis and produce a taxonomy file and a taxa barplot:

~~~
proname_taxonomy \
  --qseqs rep_seqs.qza \
  --qtable rep_table.qza \
  --db /opt/db/rEGEN-B/rEGEN-B_sequences.fasta \
  --reftax /opt/db/rEGEN-B/rEGEN-B_taxonomy.tsv \
  --metadata sample_metadata.tsv \
  --assay assay1
~~~

Here is the complete list of available arguments for `proname_taxonomy`:

| Command | Arguments | Description | Mandatory arguments |
| ------- | --------- | ----------- | ------------------- |
| proname_taxonomy | --qseqs | Path to query sequences. It should be either 'rep_seqs.qza' (if data was imported into QIIME2 at the previous step) or 'rep_seqs.fasta'. | Yes |
|  | --qtable | Path to query sequence abundance table. It should be 'rep_table.qza'. This argument is only needed if data was imported into QIIME2 at the previous step. | ~ |
|  | --db | Path to the name of the reference database used by blastn to carry out the taxonomic analysis. The rEGEN-B (rrn operons Extracted from GENomes of Bacteria) database as well as the Silva138 and Greengenes2 databases are already precompiled in the PRONAME environment and are located in the folder '$HOME/miniconda3/envs/proname/db' (Change "proname" in this path if you named your environment differently). The user can provide the name of another blastn database if desired. Note that the database must be formatted to run with the BLAST Command Line Applications (for more info, see [here](https://www.ncbi.nlm.nih.gov/books/NBK569841/)). | Yes |
|  | --reftax | Path to the taxonomic lineages of the sequences included in the reference database. The taxonomy tsv files associated with the rEGEN-B, Silva138 and Greengenes2 databases are already precompiled in the PRONAME environment and are located in the folder '$HOME/miniconda3/envs/proname/db' (Change "proname" in this path if you named your environment differently). The user can provide another reference database and associated taxonomic lineages if desired. | Yes |
|  | --evalue | Expectation value (E) threshold to keep hits. [Default: 0.001] | No |
|  | --percid | Percent identity threshold between query and subject sequences to keep hits. [Default: 80] | No |
|  | --qcover | Percent query coverage threshold to keep hits. Note that this option does not correspond to the blastn -qcov_hsp_perc option. Here, --qcover reflects the percent coverage of the whole query sequence, so that the results provided are more consistant with those obtained with the online BLASTn tool. [Default: 80] | No |
|  | --threads | Number of threads to use for the blastn analysis. You can know the number of available threads on your computer by running the command 'nproc --all' [Default: 2] | No |
|  | --metadata | Path to the metadata file. It should be 'sample_metada.tsv'. This argument is only required if data was imported into QIIME2 at the previous step. | ~ |
|  | --assay | Name of your metabarcoding assay, that will appear in the name of taxonomy files produced. | Yes |
|  | --phyloseq | Specify if a phyloseq object must be generated. [Options: "yes" or "no", Default: "no"] | No |
|  | --version | Print the version of the pipeline. | No |
|  | --verbose | Activate verbose/debug mode (no redirections). | No |
|  | --help | Print the help menu. | No |

The taxa barplot generated can be used to visualize the results of the taxonomic analysis using [QIIME2 View](https://view.qiime2.org/) or the [Dokdo API](https://dokdo.readthedocs.io/en/latest/dokdo_api.html) for instance:

![taxa_barplot](./images/taxa_barplot_tutorial.png?raw=true "taxa_barplot")

This marks the end of the PRONAME pipeline, which allowed generating high-accuracy consensus sequences starting from raw Nanopore metabarcoding data, and performing the subsequent taxonomic analysis.

## 5. After PRONAME

Since we decided to import the generated files into QIIME2, it is easy to carry on with further analyses using this platform. In this tutorial, we will focus on diversity analyses and differential abundance testing.

First of all, we will add information in the metadata file to define to which treatment group each sample belongs:

~~~
awk 'BEGIN {FS=OFS="\t"} \
  $1 ~ /sample-id/ {$3="treatment"} \
  $1 ~ /q2:types/ {$3="categorical"} \
  $1 ~ /sample1|sample2|sample3|sample4|sample5/ {$3="treatment1"} \
  $1 ~ /sample6|sample7|sample8|sample9|sample10/ {$3="treatment2"} \
  1' sample_metadata.tsv > sample_metadata.tsv_tmp && mv sample_metadata.tsv_tmp sample_metadata.tsv
~~~

### 5.1. Diversity analysis

We will generate a tree for phylogenetic diversity analyses:

~~~
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep_seqs.qza \
  --o-alignment aligned_rep_seqs.qza \
  --o-masked-alignment masked_aligned_rep_seqs.qza \
  --o-tree unrooted_tree.qza \
  --o-rooted-tree rooted_tree.qza
~~~

Before computing diversity metrics, we need to choose a rarefaction depth thanks to alpha rarefaction plotting:

~~~
qiime diversity alpha-rarefaction \
  --i-table rep_table.qza \
  --m-metadata-file sample_metadata.tsv \
  --o-visualization alpha_rarefaction_curves.qzv \
  --p-min-depth 1 \
  --p-max-depth 12000
~~~

![alpha_rarefaction](./images/alpha_rarefaction.png?raw=true "alpha_rarefaction")

Based on the previous figure, we will select 11000 as the rarefaction depth for diversity calculations:

~~~
qiime diversity core-metrics-phylogenetic \
  --i-table rep_table.qza \
  --i-phylogeny rooted_tree.qza \
  --m-metadata-file sample_metadata.tsv \
  --p-sampling-depth 11000 \
  --output-dir core-metrics-results
~~~

#### 5.1.1. Alpha diversity

Different alpha diversity metrics have been computed at the previous step, including Faith's phylogenetic diversity, the number of observed features and the Shannon index:

![alpha_diversity](./images/alpha_diversity.png?raw=true "alpha_diversity")

We can now test whether the distribution of sequences is significantly different between both treatment groups:

~~~
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file sample_metadata.tsv \
  --o-visualization alpha_observed_features_group_significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample_metadata.tsv \
  --o-visualization alpha_faith_pd_group_significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file sample_metadata.tsv \
  --o-visualization alpha_shannon_group_significance.qzv
~~~

The Kruskal-Wallis tests performed provided the following p-values:

|     | faith_pd | observed_features | shannon |
| --- | -------- | ----------------- | ------- |
| p-value | 0.602 | 0.249 | **0.028** |

#### 5.1.2. Beta diversity

Different beta diversity metrics were also computed at the above core metric calculation step and visualization files were generated. Here is an example with the unweighted UniFrac PCoA results:

![beta_diversity](./images/beta_diversity.png?raw=true "beta_diversity")

We can now analyze the statistical trends for different metrics using PERMANOVA:

~~~
# Bray Curtis
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization beta_braycurtis_treatment_significance.qzv

# Jaccard
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization beta_jaccard_treatment_significance.qzv

#Weighted UniFrac
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization beta_weighted_unifrac_treatment_significance.qzv

# Unweighted UniFrac
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample_metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization beta_unweighted_unifrac_treatment_significance.qzv
~~~

The PERMANOVA tests performed provided the following p-values:

|     | Bray Curtis | Jaccard | Weighted UniFrac | Unweighted UniFrac |
| --- | ----------- | ------- | ---------------- | ------------------ |
| p-value | 0.151 | **0.009** | **0.006** | 0.057 |

### 5.2. Differential abundance

The last step of this tutorial will test whether individual taxa show a stastistically significant enrichement or depletion in different sample groups. This is done using [ANCOM-BC](https://pubmed.ncbi.nlm.nih.gov/32665548/), that takes into account the sources of biases affecting microbiome data to perform bias correction.

In this example, we will perform the differential abundance test at a the genus taxonomic level. To do so, we need first to collapse the sequences in the `rep_table.qza` at this taxonomic level:

~~~
qiime taxa collapse \
  --i-table rep_table.qza \
  --i-taxonomy taxonomy_blastn.qza \
  --p-level 6 \
  --o-collapsed-table rep_table_l6.qza
~~~

We can now  run the ANCOM-BC analysis, taking the treatment1 as the reference group:

~~~
qiime composition ancombc \
  --i-table rep_table_l6.qza \
  --m-metadata-file sample_metadata.tsv \
  --p-formula 'treatment' \
  --p-reference-levels 'treatment::treatment1' \
  --o-differentials ancombc_treatment_l6.qza

qiime composition da-barplot \
  --i-data ancombc_treatment_l6.qza \
  --p-significance-threshold 0.05 \
  --o-visualization da_barplot_treatment_l6.qzv
~~~

![differential_abundance](./images/differential_abundance.png?raw=true "differential_abundance")


# Reference

Dubois, B., Delitte, M., Lengrand, S., Bragard, C., Legrève, A., & Debode, F. (2024). PRONAME: a user-friendly pipeline to process long-read Nanopore metabarcoding data by generating high-quality consensus sequences. Frontiers in Bioinformatics, 4, 1483255. doi:  https://doi.org/10.3389/fbinf.2024.1483255



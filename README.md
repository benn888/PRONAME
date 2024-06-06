![PRONAME_logo](./images/PRONAME_logo.jpg?raw=true "PRONAME logo")

# PRONAME: PROcessing NAnopore MEtabarcoding data

PRONAME is an open-source bioinformatics pipeline that allows processing Nanopore metabarcoding sequencing data.

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
The best way to go is to type the name of each script followed by "--help" (e.g. `proname_import --help`) to get the list of all arguments and a usage example.

# Tutorial



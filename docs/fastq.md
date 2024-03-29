## Setup

Despite its name, bam2bakR is also able to take fastq files as input, trimming and aligning them before running them through the standard bam2bakR pipeline. The current implementation of this is a bit rough around the edges and inflexible though. For aligning fastq files to create bam files that can be passed to bam2bakR, I will shamelessly suggest my own Snakemake pipeline called [THE_Aligner](https://github.com/isaacvock/THE_Aligner) which implements a number of different aligners, makes use of automated unit testing, is more configurable, and includes QC of fastqs and output bam files that bam2bakR does not currently implement. A new and improved version of the full fastq-to-cB file pipeline implemented in bam2bakR is currently under development at [fastq2EZbakR](https://github.com/isaacvock/fastq2EZbakR), and is almost fully operational (as of 11/29/2023; pipeline is fully implemented, in the process of running several test datasets). Documentation for bam2bakR's fastq processing strategy is included here for posterity's sake.

The steps for making use of bam2bakR's fastq processing functionality are pretty much identical to those for using bam files as input, but I will reproduce those steps here in their entirety for convenience sake. The main difference is that there are now additional parameters in the config to make note of. The steps are:

1. [Install conda (or mamba) on your system](#conda). This is the package manager that bam2bakR uses to make setting up the necessary dependencies a breeze.
1. [Deploy workflow](#deploy) with [Snakedeploy](https://snakedeploy.readthedocs.io/en/latest/index.html)
1. [Edit the config file](#config) (located in config/ directory of deployed/cloned repo) to your liking
1. [Run it!](#run)

The remaining documentation on this page will describe each of these steps in greater detail and point you to additional documentation that might be useful.

### Install conda (or mamba)<a name="conda"></a>
[Conda](https://docs.conda.io/projects/conda/en/latest/index.html) is a package/environment management system. [Mamba](https://mamba.readthedocs.io/en/latest/) is a newer, faster, C++ reimplementation of conda. While often associated with Python package management, lots of software, including all of the TimeLapse pipeline dependencies, can be installed with these package managers. They have pretty much the same syntax and can do the same things, so I highly suggest using Mamba in place of Conda whenever possible. 

One way to install Mamba is to first install Conda following the instructions at [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Then you can call:

``` bash
conda install -n base -c conda-forge mamba
```
to install Mamba.

A second strategy would be to install Mambaforge, which is similar to something called Miniconda but uses Mamba instead of Conda. I will reproduce the instructions to install Mambaforge below, as this is probably the easiest way to get started with the necessary installation of Mamba. These instructions come from the [Snakemake Getting Started tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html), so go to that link if you'd like to see the full original details:

* For Linux users with a 64-bit system, run these two lines of code from the terminal:

``` bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```
* For Mac users with x86_64 architecture: 
``` bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh -o Mambaforge-MacOSX-x86_64.sh
bash Mambaforge-MacOSX-x86_64.sh
```
* And for Mac users with ARM/M1 architecture:
``` bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-arm64.sh -o Mambaforge-MacOSX-arm64.sh
bash Mambaforge-MacOSX-arm64.sh
```

When asked this question:
``` bash
Do you wish the installer to preprend the install location to PATH ...? [yes|no]
```
answer with `yes`. Prepending to PATH means that after closing your current terminal and opening a new one, you can call the `mamba` (or `conda`) command to install software packages and create isolated environments. We'll be using this in the next step.

### Deploy workflow<a name="deploy"></a>

Version 1.0.1 of bam2bakR is now compatible with deployment using the tool [Snakedeploy](https://snakedeploy.readthedocs.io/en/latest/index.html). To get started with Snakedeploy, you first need to create a simple conda environment with Snakemake and Snakedeploy (Python version has to be pinned to < 3.12 due to bugs discussed [here](https://github.com/snakemake/snakemake/issues/2459) and [here](https://github.com/snakemake/snakemake/issues/2480)):


``` bash
mamba create -c conda-forge -c bioconda --name deploy_snakemake snakemake snakedeploy "python<3.12"
```

Next, create a directory that you want to run bam2bakR in (I'll refer to it as `workdir`) and move into it:
``` bash
mkdir workdir
cd workdir
```

Now, activate the `deploy_snakemake` environment and deploy the workflow as follows:

``` bash
conda activate deploy_snakemake
snakedeploy deploy-workflow https://github.com/simonlabcode/bam2bakR.git . --branch main
```

`snakedeploy deploy-workflow https://github.com/simonlabcode/bam2bakR.git` copies the content of the `config` directory in the bam2bakR Github repo into the directoy specified (`.`, which means current directory, i.e., `workdir` in this example). It also creates a directory called `workflow` that contains a singular Snakefile that instructs Snakemake to use the workflow hosted on the main branch (that is what `--branch main` determines) of the bam2bakR Github repo. `--branch main` can also be replaced with `--tag 1.0.2` to ensure that you are consistently using the same version of bam2bakR (version 1.0.2 release).

### Edit the config file<a name="config"></a>

**!!This is the part that differs from the standard bam-file input deployment!!**

In the `config/` directory you will find a file named `config.yaml`. If you open it in a text editor, you will see several parameters which you can alter to your heart's content. The first parameter that you have to set is at the top of the file:

``` yaml
samples:
  WT_1: data/fastq/WT_1
  WT_2: data/fastq/WT_2
  WT_ctl: data/fastq/WT_ctl
  KO_1: data/fastq/KO_1
  KO_2: data/fastq/KO_2
  KO_ctl: data/fastq/KO_ctl
```
`samples` is the list of sample IDs and paths to directory containing .fastq files that you want to process. **NOTE: each directory must contain a single fastq file or pair of fastq files**.  Delete the existing sample names and paths and add yours. The sample names in this example are `WT_1`, `WT_2`, `WT_ctl`, `KO_1`, `KO_2`, and `KO_ctl`. These are the sample names that will show up in the `sample` column of the output cB.csv file. The `:` is necessary to distinguish the sample name from what follows, the path to the relevant bam file. Note, the path is NOT an absolute path, it is relative to the directory that you deployed to (i.e., `workdir` in this example). Thus, in this example, the fastq file-containing directories are located in a directory called `fastq` that is inside of a directory called `data` located in `workdir`. Your data can be wherever you want it to be, but it might be easiest if you put it in a `data` directory inside the bam2bakR directory as in this example. 

As another example, imagine that the `data` directory was in the directory that contains `workdir`, and that there was no `samples` subdirectory inside of `data`. In that case, the paths would look something like this:

``` yaml
samples:
  WT_1: ../data/WT_1
  WT_2: ../data/WT_2
  WT_ctl: ../data/WT_ctl
  KO_1: ../data/KO_1
  KO_2: ../data/KO_2
  KO_ctl: ../data/KO_ctl
```
where `../` means navigate up one directory. 

The next parameter you have to set denotes the sample names of any -s4U control samples (i.e., samples that were not fed s4U or a similar metabolic label):

``` yaml
control_samples: ["WT_ctl", "KO_ctl"]
```

In this case, the samples named WT_ctl and KO_ctl are the -s4U control samples. -s4U controls will be used to call any single nucleotide polymorphisms (SNPs) in your cell line so as to avoid conflating them with T-to-C mutations induced by the nucleotide recoding chemistry. 

The third crucial parmaeter immediately follows:

``` yaml
annotation: data/annotation/GRCh38.gtf
```
This is the path to the GTF annotation file for the genome that reads were mapped to. The same rules apply when it comes to specifying this path.

Finally, the path to the genome fasta file that you used for alignment must also be specified:

``` yaml
genome_fasta: data/genome/GRCh38.fasta
```

The other parameters that can be altered are:

* `strandedness`: whether the first read in a pair (or the only read if single-end) represents the original sequence of the RNA (F), or its reverse complement (R). For example, set this parameter to "F" if your library is an FR paired-end library, and "R" if it is an RF paired-end library.
* `FORMAT`: whether the reads are paired-end (PE) or single-end (SE).
* `mut_tracks`: the type of mutation (e.g., T-to-C mutations) that sequencing tracks will be colored by. If you are most interested in the T-to-C mutational content, then `mut_tracks` should be TC. If G-to-A, then `mut_tracks` should be GA. If both, then `mut_tracks` should be "TC,GA".
* `minqual`: Minimum base quality to call it a mutation.
* `keepcols`: Names of columns to keep in cB.csv output file. See [Output](output.md) for details of columns you can keep.
* `spikename`: If spike-ins are present, this should be a string that is common to all gene_ids for spike-in transcripts in annotation gtf. For example, in Ensembl annotations for Drosophila melanogaster, all gene_ids start with "FBgn". Therefore, if you have Drosophila spike-ins, `spikename` should be "FBgn".
* `normalize`: If True, then scale factor calculated with edgeR is used to normalize sequencing tracks.
* `WSL`: whether you are running this on the Windows subsystem for linux (0 = yes; 1= no)

Configurable parameters that uniquely impact the fastq processing are:

* `HISAT2`: Path to directory containing hisat2 indices. These must be pre-built by the user. See [here](https://daehwankimlab.github.io/hisat2/manual/) for details.
* `HISAT_3N`: Path to directory containing hisat-3n indices. These must be pre-built by the user. See [here](https://daehwankimlab.github.io/hisat2/hisat-3n/) for details.
* `STAR_index`: Path to directory containing STAR indices. These can be created by the pipeline if `build_star` is set to True. They will be created at the directory specified by `STAR_index`.
* `use_hisat3n`: If True, HISAT-3N will be used for alignment.
* `use_star`: If True, and `use_hisat3n` is False, STAR will be used for alignment. If both `use_hisat3n` and `use_star` are false, then HISAT2 will be used.
* `build_star`: If True, STAR indices will be built automatically at path specified in `STAR_index` if not provided by user.
* `hisat3n_path`: Path to HISAT-3N executable; HISAT-3N cannot be installed automatically like all other dependencies as it is not installable via conda.
* `chr_tag`: If True, chr is added to chromosome names during alignment (HISAT2 and HISAT-3N only). Useful when aligner index is number-based but annotation file has "chr" in its seqnames.
* `Yale`: If True, HISAT-3N can be loaded as a module available on the McCleary cluster. Otherwise, HISAT-3N must be installed manaully by the user and the path to the exectuable provided.
* `flattened`: This parameter should be kept to False unless you know what you are doing.
* `adapter`: Code snippet that is passed to [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) for adapter trimming.
* `cutadapt_extra`: Extra, optional parameters passed to Cutadapt.
* `star_extra`: Extra, optional parameters passed to STAR for alignment.
* `hisat2_extra`: Extra, optional parameters passed to HISAT2 for alignment.



 
Edit the values in the config file as necessary and move on to the last step.

### Run it!<a name="run"></a>

Once steps 1-3 are complete, bam2bakR can be run from the directory you deployed the workflow to as follows:

``` bash
snakemake --cores all --use-conda
```
There are **A LOT** of adjustable parameters that you can play with when running a Snakemake pipeline. I would point you to the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) 
for the details on everything you can change when running the pipeline.


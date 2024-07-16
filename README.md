# Snakemake workflow: RNA sequencing with pytximport

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0.0-brightgreen.svg)](https://snakemake.bitbucket.io)

End-to-end bulk RNA sequencing analysis from paired-end FASTQ files to differential gene expression and enrichment analyses, consisting of the following steps:

1. Quality control
    - FastQC metrics of all read files
    - fastp trimming and quality control of read files
     - MultiQC report of FastQC and fastp metrics
2. Reference file preparation
    - Download of a FASTA genome from Ensembl
    - Download of a FASTA transcriptome from Ensembl
    - Download of a .gtf annotation file from Ensembl
    - Construction of a gene id to transcript id and gene name mapping table from the annotation file
3. Salmon index building for alignment-free quantification
    - Decoy preparation based on the genome and construction of a gentrome
    - Index building based on the gentrome with salmon
4. Salmon quantification-free alignment with bootstraps
5. Gene-level summarization with pytximport
6. Differential gene expression analysis with PyDESeq2
7. Gene set enrichment analysis with decoupler-py

This workflow is based on the Snakemake cookiecutter template, the kallisto-sleuth example workflow from Johannes Köster and many available Snakemake wrappers maintained by different authors.

## Authors

* Malte Kuehl (@maltekuehl)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI.

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.tsv` to specify your sample setup.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

### Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.

### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 7: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:complextissue/snakemake-bulk-rna-sequencing-workflow.git` or `git remote add -f upstream https://github.com/complextissue/snakemake-bulk-rna-sequencing-workflow.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 8: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.

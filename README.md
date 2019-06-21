# EM
###Run EM (Expectation Maximization) on AWS by Python Snakemake


(I assumed that you already have an AWS instance and a Github account.)



####1. Download this git repository to your AWS instance
 - Click on "Clone or download" on the top right side of this page  and then "copy" button
 - Enter below on AWS instance command line:
```
git clone {address that you have copied}
```

####2. Install Miniconda.
Conda is a package management system and environment management system. Click [here](https://www.biostars.org/p/335903/) for a brief introduction.
To install Miniconda, on command line, enter:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

The address might have changed. If so, go to the conda website download page, and copy the address of Python3.7, Linux, 64-bit version.

After the installation is complete, you probably need to `source ~/.bash_profile` to activate `conda` command.

####3. Install Snakemake and all other dependencies
Snakemake is a workflow management system to create reproducible and scalable data analyses. The EM workflow runs scripts written in perl and java, so we also depend on these two softwares. It's most likely that your AWS already installed perl and java, but just for version control purposes, we install the specific versions of perl and java. The versions in the `EM/env/environment.yml` were tested OK in Jun 2019.
To install the dependencies, in the EM folder, enter:

```
conda env create -f env/environment.yml
```

It may take some time to complete the installation.
After the installation is complete, activate the "EM" environment by:
```
conda activate EM
```

####4. Edit the configuration file.
The configuration file (EM/config/example\_config.yaml) specifies input directory, output directory, and what ethnic groups you want to run on. The input directory is the path to your GL\_greedy output, which contains all ethnic groups' glid and pull files. You don't need to "mkdir" the output directory, as it will be created by itself if it does not exist yet. In the `example_config.yaml` file, all 26 ethnic groups were listed but commented out, except JAPI. To include other ethnic groups, just remove the `#` in front of them.

####5. Execute the workflow.
To do a "dry" run (which only prints out the jobs to be run), in the EM directory, enter:

```
snakemake --cores 2 --configfile config/example_config.yaml -n
```

This dry-runs the workflow written in `Snakefile` file, and  `-n` stands for dry run.
If the print-outs are in green, it means there's no issues in the workflow. Otherwise, the print-outs are in red.
Omit `-n` for a real run.

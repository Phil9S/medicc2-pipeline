# medicc2-pipeline

## Authors

* Philip Smith (@phil9s)

## Pipeline setup

### Step 1 Clone the repo

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system.

```
git clone https://github.com/Phil9S/swgs-absolutecn.git
cd swgs-absolutecn/
```

### Step 2 Install conda

Run the following to install conda whilst following the on-screen instructions.
- When asked to run `conda init` and initialise conda please respond with 'yes'

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -p $HOME/miniconda/
source ~/.bashrc
rm Miniconda3-latest-Linux-x86_64.sh
```

See installing [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for more information.

#### For those with Conda already installed

For systems where conda is already available the following requirements need to be met:
- conda must be available on the PATH
- conda version `4.8.3' or greater
- the location of the installation folder is required

Check the installed version of conda using the following:
```
conda -V
```
*If this command does not work then conda is also not available on the PATH*

Find your installation directory using the following:
```
whereis conda | sed 's%condabin/conda%%'
```

### Step 3 Installing additional dependencies

From within the repository directory, run the install_env.sh script to generate a conda environment and install custom packages:
```
./install_env.sh $HOME/miniconda/

```

If you used a previously installed conda build please use the conda or miniconda installation directory when running this section instead of '$HOME/miniconda/' to correctly initialise the conda environment.

```
conda activate medicc2
```

### Step 4 Preparing the input files

#### Copy number data

The `medicc2-pipeline` accepts three input formats in order to perform multi-sample/multi-patient medicc2 analysis. These input types are `segment tables`, `QDNAseq` data objects with class `QDNAseqCopyNumbers`, or an existing folder of medicc2-suitable input files. As medicc2 requires both a single file per patient and for each sample to contain a consistent number of segments, pre-processing is applied input types which do not conform to these requirements (`segment tables` and `QDNAseq`).

#### metadata

A metadata tab-seperated file is required to correctly split samples into patient-specific medicc2 inputs found in the `segment table` and `QDNAseq` input types. This will generate one input file per sample in the metadata file containing all samples associated with that patient. The table below demonstrates the metadata schema required. The pipeline currently makes no checks that this file is valid or that sample names match between copy number data and metadata inputs. Additionally, the headers must be as written below as no method is currently implemented to alter or change the metadata table headers.

|PATIENT_ID|SAMPLE_ID|
|----------|---------|
|PAT1      |SAM1     |
|PAT1      |SAM2     |

An example meta.tsv is included in this repository.

### Step 5 Running medicc2-pipeline



# medicc2-pipeline

## Description

*Summary*

Reconstruct phylogenetic trees based on minimum-event distance (MED) for multiple samples with variable copy number input formats using the [medicc2](https://bitbucket.org/schwarzlab/medicc2/src/master/) algorithm.

## Pipeline setup
### Clone the repo
```
git clone https://github.com/Phil9S/medicc2-pipeline.git
cd medicc2-pipeline/
```
### Install conda
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
- conda version `4.8.2' or greater
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

### Installing additional dependencies

From within the repository directory, run the install_env.sh script to generate a conda environment and install custom packages:
```
./install_env.sh $HOME/miniconda/
```

If you used a previously installed conda build please use the conda or miniconda installation directory when running this section instead of '$HOME/miniconda/' to correctly initialise the conda environment.

```
conda activate medicc2
```

## Running pre-processing
### Preparing the input files

#### copy number data

The `medicc2-pipeline` accepts three input formats in order to perform multi-sample/multi-patient medicc2 analysis. These input types are segment tables, [QDNAseq](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html) data objects with class `QDNAseqCopyNumbers`, or an existing folder of suitable input files required by `medicc2`. As `medicc2` requires both a single file per patient and for each sample to contain a consistent number of segments, pre-processing is applied input types which do not conform to these requirements (segment tables and QDNAseq).

##### segment table
Segment tables should be either total or allelic-specific inputs as described below:

##### Total
|chromosome|start|end |segVal|sample|
|----------|-----|----|------|------|
|chr1      |1    |1000|1     |SAM1  |
|chr1      |1    |2000|2     |SAM2  |

##### Allele-specific
|chromosome|start|end  |segValA|segValB|sample|
|----------|-----|-----|-------|-------|------|
|chr1      |1    |1000 |1      |3      |SAM1  |
|chr1      |1    |2000 |2      |2      |SAM2  |

Example copy number input files (total and allele-specific) are provided (`resources/segment_table_example.tsv` & `resources/segment_table_allele_specific_example.tsv`, respectively).

##### QDNAseqCopyNumbers
An Rda file containing a [QDNAseq](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html) object of class `QDNAseqCopyNumber`. This should have been saved using the R function `saveRDS`. An example QDNAseq file `resources/qdnaseqmod_example_file.rds` is provided in this repository.

##### Medicc2
A folder containing a list of files (one per patient) containing the required medicc2 format. Example MEDICC2 input files (total and allele-specific) are provided `resources/medicc_input_total/` & `resources/medicc2_input_as`, respectively)..

#### metadata
A metadata tab-seperated file is required to correctly split samples into patient-specific medicc2 inputs found in the `segment table` and `QDNAseq` input types. This will generate one input file per sample in the metadata file containing all samples associated with that patient. The table below demonstrates the metadata schema required. The pipeline currently makes no checks that this file is valid or that sample names match between copy number data and metadata inputs. Additionally, the headers must be as written below as no method is currently implemented to alter or change the metadata table headers.

|PATIENT_ID|SAMPLE_ID|
|----------|---------|
|PAT1      |SAM1     |
|PAT1      |SAM2     |

An example file `resources/mapping_file_example.tsv` is included in this repository.

## Running medicc2 pipeline

The medicc2 pipeline can be run on the command line using the following;
```
./run_medicc2.sh -h
```
Which will list the options and help information.

### Example runs

#### segment table (allele-specific)
```
./run_medicc2.sh -t segment -i resources/segment_table_allele_specfic_example.tsv -o resources/example_run/ -m resources/mapping_file_example.tsv
```
#### segment table (total)
```
./run_medicc2.sh -t segment -i resources/segment_table_example.tsv -o resources/example_run/ -m resources/mapping_file_example.tsv -wn
```
#### qdnaseq
```
./run_medicc2.sh -t qdnaseq -i resources/qdnaseqmod_example_file.rds -o resources/example_run/ -m resources/mapping_file_example.tsv
```
#### medicc (allele-specific)
```
./run_medicc2.sh -t medicc -i resources/medicc_input_as/ -o resources/example_run/
```
#### medicc (total)
```
./run_medicc2.sh -t medicc -i resources/medicc_input_total/ -o resources/example_run/ -wn
```
## Output
The pipeline generates multiple outputs, the output folder contains the following;

- `medicc2_output` contains the data and plots generated by MEDICC2 during minimum event distance and phylogeny computation. Additionally, a file called `all.trees.new` which contains all the generated trees concatonated together (useful for downstream analysis).
- `input_files` contains the MEDICC2 conforming input files for each 'patient' in the submitted analysis, one file per patient grouping.
- `medicc_out.rds` is an R RDS object containing the same data as contained in the `input_files` folder.
- `norm_segs.rds` is an R RDS object containin the segment normalised data prior to reformatting to match the MEDICC2 input requirements.
- `run_medicc2.sh_log.txt` is a log file for the medicc2-pipeline run and contains all stdout and stderr. 

### Authors

* Philip Smith (@phil9s)

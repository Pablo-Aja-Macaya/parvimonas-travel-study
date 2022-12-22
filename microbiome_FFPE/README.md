# Microbiome analyses (FFPE)
16S microbiome study of "Formalin-Fixed Paraffin Embedded" samples (FFPE).
- [Requirements](#requirements)
- [Download reads and metadata](#download-reads-and-metadata)
- [QIIME2 processing](#qiime2-processing)
    - [Sidle tool](#sidle-tool)
- [Post-processing](#post-processing)
    - [Relative abundance in patient 89](#relative-abundance-in-patient-89)
    - [R Libraries Installation Instructions](#r-libraries-installation-instructions)

## Requirements
In this project, there will be required to install a base Conda environment (where another necessary tools will be installed) as well as an R interpreter to run the scripts. The table below shows the shortcuts to the sections with the guidelines for installing/running the scripts/tools:

 - [Initialize Conda](#initialize-conda)
 - [Create and install QIIME2 environment](#create-and-install-qiime2-environment)
 - [Install Sidle](#install-sidle)

### Initialize Conda:
To use Conda environment, it is mandatory inicialize in terminal command-line or inside bash script (in the last case, after the bash header): 

```sh
#!/bin/bash
conda init bash
conda activate <insert-conda-environment-name>
```

### Create and install QIIME2 environment 

```sh
# Download .yml file
wget https://data.qiime2.org/distro/core/qiime2-2021.4-py38-linux-conda.yml
conda env create -n qiime2-2021.4 --file qiime2-2021.4-py38-linux-conda.yml
rm qiime2-2021.4-py38-linux-conda.yml # OPTIONAL CLEANUP
```
NOTE: 
- Also, an available file with environment plugins and packages configuration can be found at `microbiome-FFPE/qiime2/envs/qiime2.yml`).


[In case of complications during Conda or QIIME2 version installation, go to](https://docs.qiime2.org/2021.4/install/native/)

### Install Sidle
Sidle can be installed after activate your conda environment as explained previously.
```sh
# Activate conda environment
conda activate qiime2-2022.4
# Sidle requirements
conda install dask
conda install -c conda-forge -c bioconda -c qiime2 -c defaults xmltodict
# Install rescript
pip install git+https://github.com/bokulich-lab/RESCRIPt.git
# Install Sidle
pip install git+https://github.com/jwdebelius/q2-sidle
qiime dev refresh-cache
```
NOTE: 
- Another way to install Rescript is explained at `microbiome_non-FFPE/README.rmd`.
- We have used a Conda environment, but another environments, as mamba, can be used.

## Download reads and metadata

Download reads from bioproject **PRJNA911189**. Metadata can be downloaded either from the bioproject (Using both biosample and SRA metadata), but is also available in this repository at `microbiome_FFPE/data/ffpe_metadata.tsv`.

Either way, the expected metadata for the following scripts consists in a .tsv file with the columns SampleID, id_tissue, subject_type, CCRbiome, location_en, histologic_grade_en, tnm_staging, sample_nature, sequencing_run and  control_applies_to. 
- `Metadata Columns`:
    - `SampleID`: Sample identification code.
    - `id_tissue`: Type of sample (tumor, adenoma, normal, liver-met or liver-healthy). 
    - `subject_type`: The subject id code without the tissue type (CCR89).
    - `subject`: The sample classification (ccr or control).
    - `CCRbiome`: The number associated to the subject_type code.
    - `Location`: Sample location inside gastrointestinal system. 
    - `Histological_grade`: A system to classify the differences between cancerous cells apearance.
    - `tnm_staging`: One of the systems used to stage bowel (colon and rectal) cancer (informs about the tumor size as well as how far it has spread).
    - `sample_nature`: Sample classification (ffpe or control). This will let to difference between samples types (ffpe or not ffpe in metadata with more patients).
    - `sequecing_run`: sequencing run id code (it will be used to remove negative controls contaminants from those samples with which they share sequencing run id). 
    - `control_applies_to`: 


## QIIME2 processing

The full QIIME2 workflow code can be found at the file pipeline_FFPE_Samples__QIIME2_Conda.sh, located at `microbiome_FFPE/qiime2/` path. Inside it, there would be a preprocessing before using the sidle plugin (a great tool to work with problematics samples as paraffined ones). 

### Sidle tool
Paraffined samples has been analyzed with SMURF (Short MUltiple Regions Framework) python implementation in QIIME2 Framework, called q2-sidle. This software combines sequencing results of any number of amplified 16S rRNA regions to obtain a reliable microbial profiling. [To know more, go to](https://q2-sidle.readthedocs.io/en/latest/index.html).
    
### Output Files: 
When it is done multiple output files will be available (but we must know that some outputs files turn into as input files in another steps of the pipeline). The following schema is a brief outline of all the files obtained during the execution of the pipeline:

[Sidle tool](#sidle-tool)
- `regional-databases`:
    - `silva-silvaversion-Rx-trimmedlength.qza`
- `asvs`: 
    - `cutadapt`:
        - `cutadapt-Rx.qza`
        - `cutadapt-Rx.qzv`
    - `dada2`:
        - `dada2-table.qza`
        - `rep-seqs-dada2.qza`

[Post-processing](#post-processing)
- `post-dada2-trimming`:
    - `rep-seqs-dada2-trimmedlength.qza`
    - `table-dada2-trimmedlength.qza`
- `alignment`:
  - `alignment-Rx-database.qza`
- `reconstruction`:
    - `reconstruct-database`:
      - `map-file-[...].qza`
      - `summary-file-[...].qza`
    - `reconstruct-counts`:
      - `reconstructed-table-[...].qza`
- `taxonomy`:  
  - `reconstructed-taxonomy-table.qza`
      
In the dada2 folder as well as in the reconstruction folder, there could be more output files. In the case of dada2, there is a denoising file, which gives a summary of the filtering processing step in this tool. Also, in reconstruction step, there will be a map and summary files in the first step 
NOTES:
- The x in Rx represents the number of the amplified region. As we have amplified five regions, so x go from 1 to 5. 
- The [...] in the file names means that the file name is given by the person who runs the code.  

It is optional to export results to biom or tsv format. In this case, we used a R library, QIIME2R, which let to import .qza files without exported it to other formats. But, another option is export from qza format to tsv or biom and run the rmarkdown file. To see how rmarkdown would look like, go to `microbiome_non-FFPE/post-processing/ccr89.Rmd` and check until control contaminants removing. The above section will not change. 

NOTE: During the pipeline_FFPE_Samples__QIIME2_Conda.sh, there are database files required throughout the full pipeline. The sidle documentation page indicates how to prepare the chosen database (SILVA, GreenGenes, etc) to the different regions of the study samples.
    
## Post-processing
### Relative abundance in patient 89
A rmarkdown file, located at `microbiome_non-FFPE/post-processing/CCR89_FFPE_Replicates.Rmd`, can be executed to obtain the data results shows in html format.
We have used Rstudio to run R scripts.(See Readme available at `microbiome_non-FFPE/README.rmd`). 
Note: For results replication purpose, the metadata used in the rmarkdown file should be the one called metadataccr89_FFPE.tsv. This one is a simplification of the metadata called ffpe_metadata.tsv (both of them are located at: `microbiome_FFPE/data/`)

### R Libraries Installation Instructions    
In case of any problems with the libraries required in  the R post-processing script, installation instructions for the libraries can be found in the R script libraries_installation_instructions.R, located at `microbiome_FFPE/post-processing/libraries_installation_instructions.R`. 







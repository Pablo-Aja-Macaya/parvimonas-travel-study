# Microbiome analyses (FFPE)
16S microbiome study of "Formalin-Fixed Paraffin Embedded" samples (FFPE).
- [Requirements](#requirements)
- [Download reads and metadata](#download-reads-and-metadata)
- [QIIME2 processing](#qiime2-processing)
    -[Sidle tool](#sidle-tool)
- [Post-processing](#post-processing)
    - [Relative abundance in patient 89](#relative-abundance-in-patient-89)
    - [R Libraries Installation Instructions](#r-libraries-installation-instructions)

## Requirements
Each step will have different requirements, but a base Conda environment, (where QIIME2 is installed), will be used to run the scripts. To go quicker access to installation instructions sections:
 - [Initialize Conda](#install-conda)
 - [Create and install QIIME2 environment](#create-and-install-qiime2-environment)
 - [Install Sidle](#install-sidle)

### Initialize Conda:
To inicialize Conda environment (it should be included inside bash script, after the bash header):

```sh
#!/bin/bash
conda init bash
conda activate <insert-conda-environment-name>
```

### Create and install QIIME2 environment (also available at `microbiome-FFPE/qiime2/envs/qiime2.yml`):
```sh
# Download .yml file
wget https://data.qiime2.org/distro/core/qiime2-2021.4-py38-linux-conda.yml
conda env create -n qiime2-2021.4 --file qiime2-2021.4-py38-linux-conda.yml
rm qiime2-2021.4-py38-linux-conda.yml # OPTIONAL CLEANUP
```
[In case of complications with conda or qiime2 version installation, go to](https://docs.qiime2.org/2021.4/install/native/)

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
Another way to install Rescript is explained at microbiome_non-FFPE/README.rmd.



## Download reads and metadata

Download reads from bioproject **PRJNA911189**. Metadata can be downloaded either from the bioproject (Using both biosample and SRA metadata), but is also available in this repository at `microbiome_FFPE/data/ffpe_metadata.tsv`.

Either way, the expected metadata for the following scripts consists in a .tsv file with the columns SampleID, id_tissue, subject_type, CCRbiome, location_en, histologic_grade_en, tnm_staging, sample_nature, sequencing_run and  control_applies_to. Where id_tissue_sp specify the type of sample (tumor, adenoma, normal, liver-met or liver-healthy). 
As for the columns subject_type is the subject classification (ccr or control), CCRbiome gets the subject id code without the sample type as well as location, histological_grade and tnm_staging shows the sample location, the histological grade and the tnm_staging.
Finally,the columns sample_nature, control_id_seq_run and control_applies_to are related to sample classification (ffpe or control), to sequencing run code and negative control identification ("control-applies-to*" will be used to remove contamination of specific samples/tissues in their respective sequencing runs), respectively.

## QIIME2 processing
### Sidle tool
Paraffined samples has been analyzed with SMURF (Short MUltiple Regions Framework) python implementation in QIIME2 Framework, called q2-sidle. This software combines sequencing results of any number of amplified 16S rRNA regions to obtain a reliable microbial profiling. [To know more, go to]:(https://q2-sidle.readthedocs.io/en/latest/index.html).

A bash script is located at `microbiome_FFPE/qiime2/<insert-bash-script-name>.sh`.    
    
### Output Files: 
When it is done multiple outputs will be available, which will be the ones used in:

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

    
## Post-processing
### Relative abundance in patient 89
A rmarkdown, located at `microbiome_non-FFPE/post-processing/CCR89_FFPE_Replicates.Rmd`, can be executed to obtain the data results shows in html format.
We have used Rstudio to run R scripts, but there are another environments used, as mamba.(See Readme available at `microbiome_non-FFPE/README.rmd`). 
Note: For results replication purpose, the metadata used in the rmarkdown file should be the one called metadataccr89_FFPE.tsv. This one is a simplification of the metadata called ffpe_metadata.tsv (both of them are located at: `microbiome_FFPE/data/`)

### R Libraries Installation Instructions    
In case of any problems with the libraries required in  the R post-processing script, installation instructions for the libraries can be found in the R script libraries_installation_instructions.R, located at `microbiome_FFPE/post-processing/libraries_installation_instructions.R`. 







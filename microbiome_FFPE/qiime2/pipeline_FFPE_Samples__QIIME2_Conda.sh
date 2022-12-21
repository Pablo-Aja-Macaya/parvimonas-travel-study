#!/bin/bash

###---- PIPELINE from FFPE SAMPLES 

# READ PREPARATION TUTORIAL: https://q2-sidle.readthedocs.io/en/latest/read_preparation.html

# ---- Conda bash activation ---- 

conda init bash 

conda activate qiime2-2021.4 # Qiime version 


set -e # Exit if error



# ---- Import ----

# For any questions, check: https://docs.qiime2.org/2022.8/tutorials/importing/


qiime tools  import \
 --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path muestrasccr89jul22/ \
 --input-format CasavaOneEightSingleLanePerSampleDirFmt \
 --output-path pacienteccr89_SMURF.qza


# ---- CUTADAPT ----

# In case of doubt, read: https://docs.qiime2.org/2022.2/plugins/available/cutadapt/trim-paired/ 

# R1

qiime cutadapt trim-paired         --i-demultiplexed-sequences pacienteccr89_SMURF.qza         --p-front-f 'TGGCGAACGGGTGAGTAA'         --p-front-r   'CCGTGTCTCAGTCCCARTG'   --p-match-read-wildcards         --p-discard-untrimmed         --p-match-adapter-wildcards         --p-cores 8         --o-trimmed-sequences CCR-smurf-con-controles-neg-PE-paciente-ccr89_R1.qza   --verbose

##  visualization 
qiime demux summarize --i-data CCR-smurf-con-controles-neg-PE-paciente-ccr89_R1.qza --o-visualization CCR-smurf-con-controles-neg-PE-paciente-ccr89_R1.qzv



# R2

qiime cutadapt trim-paired         --i-demultiplexed-sequences  pacienteccr89_SMURF.qza        --p-front-f 'ACTCCTACGGGAGGCAGC'         --p-front-r 'GTATTACCGCGGCTGCTG'   --p-match-read-wildcards         --p-discard-untrimmed         --p-match-adapter-wildcards         --p-cores 8         --o-trimmed-sequences CCR-smurf-con-controles-neg-PE-paciente-ccr89_R2.qza   --verbose

## visualization 

qiime demux summarize --i-data CCR-smurf-con-controles-neg-PE-paciente-ccr89_R2.qza --o-visualization CCR-smurf-con-controles-neg-PE-paciente-ccr89_R2.qzv


# R3 

qiime cutadapt trim-paired         --i-demultiplexed-sequences   pacienteccr89_SMURF.qza       --p-front-f 'GTGTAGCGGTGRAATGCG'         --p-front-r 'CCCGTCAATTCMTTTGAGTT'   --p-match-read-wildcards         --p-discard-untrimmed         --p-match-adapter-wildcards         --p-cores 8         --o-trimmed-sequences CCR-smurf-con-controles-neg-PE-paciente-ccr89_R3.qza   --verbose

## visualization

qiime demux summarize --i-data CCR-smurf-con-controles-neg-PE-paciente-ccr89_R3.qza --o-visualization CCR-smurf-con-controles-neg-PE-paciente-ccr89_R3.qzv

# R4

qiime cutadapt trim-paired         --i-demultiplexed-sequences  pacienteccr89_SMURF.qza        --p-front-f 'GGAGCATGTGGWTTAATTCGA'         --p-front-r 'CGTTGCGGGACTTAACCC'   --p-match-read-wildcards         --p-discard-untrimmed         --p-match-adapter-wildcards         --p-cores 8         --o-trimmed-sequences CCR-smurf-con-controles-neg-PE-paciente-ccr89_R4.qza   --verbose


## visualization
qiime demux summarize --i-data CCR-smurf-con-controles-neg-PE-paciente-ccr89_R4.qza --o-visualization CCR-smurf-con-controles-neg-PE-paciente-ccr89_R4.qzv


# R5 

qiime cutadapt trim-paired   --i-demultiplexed-sequences pacienteccr89_SMURF.qza        --p-front-f 'GGAGGAAGGTGGGGATGAC'         --p-front-r 'AAGGCCCGGGAACGTATT'   --p-match-read-wildcards         --p-discard-untrimmed         --p-match-adapter-wildcards         --p-cores 8         --o-trimmed-sequences CCR-smurf-con-controles-neg-PE-paciente-ccr89_R5.qza   --verbose

## visualization

qiime demux summarize --i-data CCR-smurf-con-controles-neg-PE-paciente-ccr89_R5.qza --o-visualization CCR-smurf-con-controles-neg-PE-paciente-ccr89_R5.qzv



# ---- DADA2 ----
dada2_CCR89_withfilters="dada2_CCR89_withfilters"
[ -d "$dada2_CCR89_withfilters" ]  || mkdir "$dada2_CCR89_withfilters" # Create folder to save dada2 output files

# R1

qiime dada2 denoise-paired --i-demultiplexed-seqs CCR-smurf-con-controles-neg-PE-paciente-ccr89_R1.qza --p-trim-left-f 0 --p-trim-left-r 11 --p-trunc-len-f 134 --p-trunc-len-r 155 --o-table "$dada2_CCR89_withfilters" /table-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros.qza --o-representative-sequences "$dada2_CCR89_withfilters" /rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros.qza --o-denoising-stats "$dada2_CCR89_withfilters" /stats-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros.qza  --verbose 
qiime metadata tabulate --m-input-file "$dada2_CCR89_withfilters"/stats-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros.qza  --o-visualization "$dada2_CCR89_withfilters"/stats-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros.qzv

# R2

qiime dada2 denoise-paired --i-demultiplexed-seqs CCR-smurf-con-controles-neg-PE-paciente-ccr89_R2.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 129 --p-trunc-len-r 133 --o-table "$dada2_CCR89_withfilters" /table-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros.qza --o-representative-sequences "$dada2_CCR89_withfilters" /rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros.qza --o-denoising-stats "$dada2_CCR89_withfilters" /stats-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros.qza  --verbose 

qiime metadata tabulate --m-input-file "$dada2_CCR89_withfilters"/stats-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros.qza  --o-visualization "$dada2_CCR89_withfilters"/stats-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros.qzv


# R3
qiime dada2 denoise-paired --i-demultiplexed-seqs CCR-smurf-con-controles-neg-PE-paciente-ccr89_R3.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 133 --p-trunc-len-r 131 --o-table "$dada2_CCR89_withfilters" /table-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros.qza --o-representative-sequences "$dada2_CCR89_withfilters" /rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros.qza --o-denoising-stats "$dada2_CCR89_withfilters" /stats-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros.qza  --verbose 

qiime metadata tabulate --m-input-file "$dada2_CCR89_withfilters"/stats-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros.qza  --o-visualization "$dada2_CCR89_withfilters"/stats-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros.qzv 

# R4 

qiime dada2 denoise-paired --i-demultiplexed-seqs CCR-smurf-con-controles-neg-PE-paciente-ccr89_R4.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 125 --p-trunc-len-r 132 --o-table "$dada2_CCR89_withfilters" /table-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros.qza --o-representative-sequences "$dada2_CCR89_withfilters" /rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros.qza --o-denoising-stats "$dada2_CCR89_withfilters" /stats-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros.qza  --verbose 

qiime metadata tabulate --m-input-file "$dada2_CCR89_withfilters"/stats-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros.qza  --o-visualization "$dada2_CCR89_withfilters"/stats-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros.qzv 

# R5 

qiime dada2 denoise-paired --i-demultiplexed-seqs CCR-smurf-con-controles-neg-PE-paciente-ccr89_R5.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 131 --p-trunc-len-r 132 --o-table "$dada2_CCR89_withfilters" /table-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros.qza --o-representative-sequences "$dada2_CCR89_withfilters" /rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros.qza --o-denoising-stats "$dada2_CCR89_withfilters" /stats-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros.qza  --verbose 

qiime metadata tabulate --m-input-file "$dada2_CCR89_withfilters"/stats-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros.qza  --o-visualization "$dada2_CCR89_withfilters"/stats-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros.qzv


# ---- POST-DADA2 TRIMMING ----

# TRIMMING TO 100 NTS

# R1

qiime sidle trim-dada2-posthoc \
 --i-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros.qza \
 --i-representative-sequences "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros.qza \
 --p-trim-length 100 \
 --o-trimmed-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros-100nt.qza \
 --o-trimmed-representative-sequences "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros-100nt.qza 

qiime feature-table summarize --i-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros-100nt.qza --o-visualization "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros-100nt.qzv 

qiime feature-table tabulate-seqs --i-data "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros-100nt.qza --o-visualization "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros-100nt.qzv 


# R2 

qiime sidle trim-dada2-posthoc \
 --i-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros.qza \
 --i-representative-sequences "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros.qza \
 --p-trim-length 100 \
 --o-trimmed-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros-100nt.qza \
 --o-trimmed-representative-sequences "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros-100nt.qza 

qiime feature-table summarize --i-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros-100nt.qza --o-visualization "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros-100nt.qzv 
qiime feature-table tabulate-seqs --i-data "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros-100nt.qza --o-visualization "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros-100nt.qzv 


# R3 

qiime sidle trim-dada2-posthoc \
 --i-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros.qza \
 --i-representative-sequences "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros.qza \
 --p-trim-length 100 \
 --o-trimmed-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros-100nt.qza \
 --o-trimmed-representative-sequences "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros-100nt.qza 

qiime feature-table summarize --i-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros-100nt.qza --o-visualization "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros-100nt.qzv 

qiime feature-table tabulate-seqs --i-data "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros-100nt.qza --o-visualization "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros-100nt.qzv 



# R4 

qiime sidle trim-dada2-posthoc \
 --i-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros.qza \
 --i-representative-sequences "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros.qza \
 --p-trim-length 100 \
 --o-trimmed-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros-100nt.qza \
 --o-trimmed-representative-sequences "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros-100nt.qza 

qiime feature-table summarize --i-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros-100nt.qza --o-visualization "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros-100nt.qzv 

qiime feature-table tabulate-seqs --i-data "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros-100nt.qza --o-visualization "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros-100nt.qzv 


# R5 

qiime sidle trim-dada2-posthoc \
 --i-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros.qza \
 --i-representative-sequences "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros.qza \
 --p-trim-length 100 \
 --o-trimmed-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros-100nt.qza \
 --o-trimmed-representative-sequences "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros-100nt.qza 

qiime feature-table summarize --i-table "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros-100nt.qza --o-visualization "$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros-100nt.qzv 

qiime feature-table tabulate-seqs --i-data "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros-100nt.qza --o-visualization "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros-100nt.qzv 



# ---- ALIGNMENT ----
# Create folder where saved output alignment files
alignment_smurfCCR89_silva138_dbagosto2021="alignment_smurfCCR89_silva138_dbagosto2021"
[ -d "$alignment_smurfCCR89_silva138_dbagosto2021" ]  || mkdir "$alignment_smurfCCR89_silva138_dbagosto2021" # Create folder to save alignment output files


for ((i=1; i<=5; i++)); do

      	echo "Región $i"
        echo "Alignment step of regional rep-seqs with regional database kmers files."
	      echo "DADA2 files and database kmers files goes from R1 to R5 and database kmers files."

        qiime sidle align-regional-kmers \
          --i-kmers /home/elsa/SMURF/silva_sidle_agosto2021/R$i-sidle-silva-db-100nt-kmers.qza  \
          --i-rep-seq "$dada2_CCR89_withfilters"/rep-seqs-CCR-smurf-con-controles-neg-PE-CCR89_R$i-confiltros-100nt.qza \
          --p-region Region_R$i \
          --p-n-workers  6  \
          --o-regional-alignment "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R$i-silva138-dbagosto2021-R$i-100nt.qza \
	        --verbose

        qiime metadata tabulate  \
          --m-input-file "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R$i-silva138-dbagosto2021-R$i-100nt.qza  \
          --o-visualization "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R$i-silva138-dbagosto2021-R$i-100nt.qzv


done


# ---- DATABASE RECONSTRUCTION ----

echo "Database Reconstruction Step  Silva138"

#######################################
# ---- Memory Usage Error Exit ----
#######################################

# if --p-block-size 1000 does not work, probe with 100. 
# --p-block-size 100

reconstruction_memoryrefactor_silva138julio22="reconstruction_memoryrefactor_silva138julio22"
[ -d "$reconstruction_memoryrefactor_silva138julio22" ]  || mkdir "$reconstruction_memoryrefactor_silva138julio22" # Create folder to save alignment output files


qiime sidle reconstruct-database \
  --p-region Region_R0 \
    --i-kmer-map /home/elsa/SMURF/silva_sidle_agosto21_abril22/R0-sidle-silva-db-100nt-map.qza   \
    --i-regional-alignment  "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R1-silva138-dbagosto2021-R1-100nt.qza  \
  --p-region Region_R1\
    --i-kmer-map  /home/elsa/SMURF/silva_sidle_agosto21_abril22/R1-sidle-silva-db-100nt-map.qza  \
    --i-regional-alignment  "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R2-silva138-dbagosto2021-R2-100nt.qza   \
  --p-region Region_R2 \
    --i-kmer-map /home/elsa/SMURF/silva_sidle_agosto21_abril22/R2-sidle-silva-db-100nt-map.qza   \
    --i-regional-alignment  "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R3-silva138-dbagosto2021-R3-100nt.qza  \
  --p-region Region_R3  \
    --i-kmer-map  /home/elsa/SMURF/silva_sidle_agosto21_abril22/R3-sidle-silva-db-100nt-map.qza  \
    --i-regional-alignment  "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R4-silva138-dbagosto2021-R4-100nt.qza    \
  --p-region Region_R4  \
    --i-kmer-map  /home/elsa/SMURF/silva_sidle_agosto21_abril22/R4-sidle-silva-db-100nt-map.qza  \
    --i-regional-alignment  "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R5-silva138-dbagosto2021-R5-100nt.qza   \
  --p-block-size  1000  \
  --p-n-workers  6  \
  --o-database-map "$reconstruction_memoryrefactor_silva138julio22"/map-silva138julio22-memoryrefactor.qza \
  --o-database-summary "$reconstruction_memoryrefactor_silva138julio22"/summary-silva138julio22-memoryrefactor.qza \
  --verbose


echo "Database Reconstruction Step Silva138"

qiime metadata tabulate \
  --m-input-file "$reconstruction_memoryrefactor_silva138julio22"/summary-silva138julio22-memoryrefactor.qza  \
  --o-visualization "$reconstruction_memoryrefactor_silva138julio22"/summary-silva138julio22-memoryrefactor.qzv

# ---- TABLE RECONSTRUCTION ----

echo "Counts Reconstruction Step Silva138"

qiime sidle reconstruct-counts \
  --p-region Region_R1  \
   --i-regional-alignment "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R1-silva138-dbagosto2021-R1-100nt.qza \
   --i-regional-table /home/elsa/CCR89julio22/SMURF/"$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R1-confiltros-100nt.qza \
  --p-region Region_R2 \
   --i-regional-alignment "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R2-silva138-dbagosto2021-R2-100nt.qza  \
   --i-regional-table /home/elsa/CCR89julio22/SMURF/"$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R2-confiltros-100nt.qza  \
  --p-region Region_R3 \
    --i-regional-alignment "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R3-silva138-dbagosto2021-R3-100nt.qza  \
    --i-regional-table  /home/elsa/CCR89julio22/SMURF/"$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R3-confiltros-100nt.qza \ 
  --p-region Region_R4 \
    --i-regional-alignment "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R4-silva138-dbagosto2021-R4-100nt.qza  \
    --i-regional-table /home/elsa/CCR89julio22/SMURF/"$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R4-confiltros-100nt.qza \
  --p-region Region_R5 \
    --i-regional-alignment "$alignment_smurfCCR89_silva138_dbagosto2021"/align-map-R5-silva138-dbagosto2021-R5-100nt.qza  \
    --i-regional-table /home/elsa/CCR89julio22/SMURF/"$dada2_CCR89_withfilters"/table-CCR-smurf-con-controles-neg-PE-CCR89_R5-confiltros-100nt.qza  \
 --i-database-map "$reconstruction_memoryrefactor_silva138julio22"/map-smurfCCR89-julio22.qza \
 --i-database-summary "$reconstruction_memoryrefactor_silva138julio22"/summary-smurfCCR89-julio22.qza  \
 --p-block-size 1000 \
 --p-debug   \
 --p-n-workers 6 \
 --o-reconstructed-table "$reconstruction_memoryrefactor_silva138julio22"/reconstructed-table-CCR-smurf-con-controles-neg-PE-CCR89_julio22.qza  \
 --verbose


echo "Table Reconstruction Step Silva138"

# ---- TAXONOMY  ----

echo "Extracting taxonomy from reconstructed table"

qiime sidle reconstruct-taxonomy \
 --i-reconstruction-map "$reconstruction_memoryrefactor_silva138julio22"/map-smurfCCR89-julio22.qza  \
 --i-taxonomy /home/elsa/SMURF/silva_sidle_agosto2021/silva-138-99-tax.qza \
 --p-database 'silva' \
 --p-define-missing 'inherit' \
 --o-reconstructed-taxonomy "$reconstruction_memoryrefactor_silva138julio22"/reconstructed-taxonomy-CCR-smurf-con-controles-neg-PE-CCR89_julio22.qza \
 --verbose

qiime metadata tabulate \
 --m-input-file "$reconstruction_memoryrefactor_silva138julio22"/reconstructed-taxonomy-CCR-smurf-con-controles-neg-PE-CCR89_julio22.qza  \
 --o-visualization  "$reconstruction_memoryrefactor_silva138julio22"/reconstructed-taxonomy-CCR-smurf-con-controles-neg-PE-CCR89_julio22.qzv


# With silva128, you can construct a phylogenetic tree


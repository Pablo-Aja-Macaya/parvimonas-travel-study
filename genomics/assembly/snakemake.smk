import sys
sys.path.append(".")

# Create DAG tree visualization
# snakemake -s snakefile_illumina.smk --rulegraph | dot -Tpng > dag.png

###############################################
# ------- Libraries/Functions set up -------- #
###############################################

import glob
import json
import os
import pandas as pd
import numpy as np
import re
from subprocess import call, run
import pathlib
import datetime
from colorama import Fore, Style
from snakemake.logging import logger
from pathlib import Path

###############################################
# ------------ Config variables ------------- #
###############################################
"""
- `BASE_FOLDER`: Output folder (Full path)
- `SHORT_READS_INPUT`: Input folder with subfolders containing short reads (Full path). Structure inside the folder must be like `path/sr_reads/{sample_id}/{sample_id}_R{1,2}.fastq.gz`.
- `LONG_READ_INPUT`: Input folder with subfolders containing long reads (Full path). Structure inside the folder must be like `path/lr_reads/{sample_id}/{sample_id}.fastq.gz`.
- `TRIMMOMATIC_CLIP`: Trimmomatic clip in FASTA format file like "NexteraPE-PE.fa" (Full path)
- `PHIX_REF`: PhiX reference genome in FASTA format (Full path)
"""

# Output folder for pipeline
BASE_FOLDER = os.path.abspath("output/assembly")

# Reads folder
SHORT_READS_INPUT = os.path.abspath("/home/usuario/Proyectos/Results/Kelly/Parvimonas/pmicra_merged__unicycler_hybrid/reads_qc/01_prepared")
LONG_READ_INPUT = os.path.abspath("/home/usuario/Proyectos/Results/Kelly/Parvimonas/pmicra_merged__unicycler_hybrid/00_long_clean_reads")

# Others
TRIMMOMATIC_CLIP = os.path.abspath("/home/usuario/Proyectos/Innova/data/illumina_clips/NexteraPE-PE.fa")
PHIX_REF = os.path.abspath("")

# Environments
CONDA_ENVS = {
    # If these strings finish in .yml/.yaml snakemake will install the environment
    # If they do not, they will be considered as names of already installed environments
    "main": "envs/parvimonas_assembly.yml",
}

###############################################
# ----------------- Input ------------------- #
###############################################

# ---- Sample wildcards ----
wc_tmp = glob_wildcards(Path(SHORT_READS_INPUT, "{xxx}", "{barcode}_R1.fastq.gz"))
barcodes = wc_tmp.barcode

# Check that each barcodes has exactly 2 fastq files
# If it raises an error one barcode is duplicated and files will need to be renamed
for bc in set(barcodes):
    if barcodes.count(bc) != 1:
        raise Exception(f"Unexpected")

barcodes = set(barcodes)

# ---- Set wildcards dict ----
wildcards_dict = {"barcode": barcodes}

# ---- Log to user ----
logger.info("")
logger.info(" Pipeline main variables ".center(70, "~"))
logger.info(f"{Fore.LIGHTBLACK_EX}- Selected barcodes: {wildcards_dict['barcode']}{Style.RESET_ALL}")
logger.info(f"{Fore.LIGHTBLACK_EX}- Short reads input folder: {SHORT_READS_INPUT}{Style.RESET_ALL}")
logger.info(f"{Fore.LIGHTBLACK_EX}- Long reads input folder: {LONG_READ_INPUT}{Style.RESET_ALL}")
logger.info(f"{Fore.LIGHTBLACK_EX}- Output folder: {BASE_FOLDER}{Style.RESET_ALL}")

###############################################
# -------------- Folder set up -------------- #
###############################################

# Cleaning short reads
SHORT_READS_PIPELINE_OUTPUT = f"{BASE_FOLDER}/short_reads_qc"
PREPARED_SHORT_READS_OUTPUT = f"{SHORT_READS_PIPELINE_OUTPUT}/01_prepared"
FASTQC_BEFORE_OUTPUT = f"{SHORT_READS_PIPELINE_OUTPUT}/02_fastqc_before"
BBDUK_OUTPUT = f"{SHORT_READS_PIPELINE_OUTPUT}/03_bbduk"
CLUMP_OUTPUT = f"{SHORT_READS_PIPELINE_OUTPUT}/04_clumpify"
TRIMMOMATIC_OUTPUT = f"{SHORT_READS_PIPELINE_OUTPUT}/05_trimmomatic"
FASTP_OUTPUT = f"{SHORT_READS_PIPELINE_OUTPUT}/05_fastp"
BMTAGGER_OUTPUT = f"{SHORT_READS_PIPELINE_OUTPUT}/06_bmtagger"
SHORT_CLEAN_READS_OUTPUT = f"{SHORT_READS_PIPELINE_OUTPUT}/06_clean_reads"
FASTQC_AFTER_OUTPUT = f"{SHORT_READS_PIPELINE_OUTPUT}/07_fastqc_after"
MULTIQC_OUTPUT = f"{SHORT_READS_PIPELINE_OUTPUT}/08_multiqc"

# Cleaning long reads
LONG_READS_PIPELINE_OUTPUT = f"{BASE_FOLDER}/long_reads_qc"
PREPARED_LONG_READS_OUTPUT = f"{LONG_READS_PIPELINE_OUTPUT}/01_prepared"
NANOPLOT_QC_BEFORE_OUTPUT = f"{LONG_READS_PIPELINE_OUTPUT}/02_nanoplot_before"
PORECHOP_OUTPUT = f"{LONG_READS_PIPELINE_OUTPUT}/03_porechop"
FILTLONG_OUTPUT = f"{LONG_READS_PIPELINE_OUTPUT}/04_filtlong"
NANOPLOT_QC_AFTER_OUTPUT = f"{LONG_READS_PIPELINE_OUTPUT}/05_nanoplot_after"
LONG_CLEAN_READS_OUTPUT = f"{LONG_READS_PIPELINE_OUTPUT}/06_clean_reads"

# Assembly
ASSEMBLY_PIPELINE_OUTPUT = f"{BASE_FOLDER}/assembly"
ASSEMBLY_OUTPUT = f"{ASSEMBLY_PIPELINE_OUTPUT}/01_assembly"
QUAST_OUTPUT = f"{ASSEMBLY_PIPELINE_OUTPUT}/02_quast"
COVERAGE_OUTPUT = f"{ASSEMBLY_PIPELINE_OUTPUT}/03_coverage"
FILTERED_ASSEMBLY_OUTPUT = f"{ASSEMBLY_PIPELINE_OUTPUT}/04_filtered_assembly"
FINAL_ASSEMBLY_OUTPUT = f"{ASSEMBLY_PIPELINE_OUTPUT}/05_final_assembly"
FILTERED_QUAST_OUTPUT = f"{ASSEMBLY_PIPELINE_OUTPUT}/06_filtered_quast"
CHECKM_OUTPUT = f"{ASSEMBLY_PIPELINE_OUTPUT}/07_checkm"

###############################################
# ---------------- Output ------------------- #
###############################################


rule all:
    input:
        # Clean reads
        expand(PREPARED_SHORT_READS_OUTPUT + "/{barcode}", **wildcards_dict),
        expand(PREPARED_LONG_READS_OUTPUT + "/{barcode}", **wildcards_dict),

        # QC plots for Illumina reads
        expand(FASTQC_BEFORE_OUTPUT + "/{barcode}", **wildcards_dict),
        expand(FASTQC_AFTER_OUTPUT + "/{barcode}", **wildcards_dict),
        MULTIQC_OUTPUT + "/after",

        # QC plots for ONT reads
        expand(NANOPLOT_QC_BEFORE_OUTPUT + "/{barcode}", **wildcards_dict),
        expand(NANOPLOT_QC_AFTER_OUTPUT + "/{barcode}", **wildcards_dict),

        # Assembly
        expand(ASSEMBLY_OUTPUT + "/{barcode}/{barcode}_assembly.fasta", **wildcards_dict),



###############################################
# -------- Short reads processing ----------- #
###############################################
rule prepare_short_reads:
    input:
        reads_1 = SHORT_READS_INPUT + '/{barcode}/{barcode}_R1.fastq.gz',
        reads_2 = SHORT_READS_INPUT + '/{barcode}/{barcode}_R2.fastq.gz',
    output:
        folder=directory(PREPARED_SHORT_READS_OUTPUT + "/{barcode}"),
        reads_1=PREPARED_SHORT_READS_OUTPUT + "/{barcode}/{barcode}_R1.fastq.gz",
        reads_2=PREPARED_SHORT_READS_OUTPUT + "/{barcode}/{barcode}_R2.fastq.gz",
    shell:
        """
        mkdir -p {output.folder}
        cp {input.reads_1} {output.reads_1}
        cp {input.reads_2} {output.reads_2}
        """

rule fastqc_before:
    """FastQC quality control"""
    input:
        rules.prepare_short_reads.output.reads_1,
        rules.prepare_short_reads.output.reads_2,
    output:
        folder=directory(FASTQC_BEFORE_OUTPUT + "/{barcode}"),
    shell:
        """
        mkdir -p {output}
        fastqc {input} -o {output.folder}
        """

rule bbduk:
    """
    Extract lambda phage PhiX
    """
    input:
        reads_1=rules.prepare_short_reads.output.reads_1,
        reads_2=rules.prepare_short_reads.output.reads_2,
    output:
        folder=directory(BBDUK_OUTPUT + "/{barcode}"),
        reads_1=temp(BBDUK_OUTPUT + "/{barcode}/{barcode}_R1.fastq.gz"),
        reads_2=temp(BBDUK_OUTPUT + "/{barcode}/{barcode}_R2.fastq.gz"),
        stats=BBDUK_OUTPUT + "/{barcode}/stats.txt",
    params:
        phix_ref=PHIX_REF,
    conda:
        CONDA_ENVS["main"]
    shell:
        """
        mkdir -p {output.folder}
        bbduk.sh in={input.reads_1} in2={input.reads_2} out={output.reads_1} out2={output.reads_2} ref={params.phix_ref} stats={output.stats}
        """


rule clumpify:
    """
    Eliminate duplicates / "compress" reads
    """
    input:
        reads_1=rules.bbduk.output.reads_1,
        reads_2=rules.bbduk.output.reads_2,
    output:
        folder=directory(CLUMP_OUTPUT + "/{barcode}"),
        reads_1=temp(CLUMP_OUTPUT + "/{barcode}/{barcode}_R1.fastq.gz"),
        reads_2=temp(CLUMP_OUTPUT + "/{barcode}/{barcode}_R2.fastq.gz"),
    conda:
        CONDA_ENVS["main"]
    shell:
        """
        mkdir -p {output.folder}
        clumpify.sh -Xmx1g dedupe=t in={input.reads_1} in2={input.reads_2} out={output.reads_1} out2={output.reads_2} reorder
        """


rule trimmomatic:
    """
    Trim and quality control reads
    """
    input:
        reads_1=rules.clumpify.output.reads_1,
        reads_2=rules.clumpify.output.reads_2,
    output:
        folder=directory(TRIMMOMATIC_OUTPUT + "/{barcode}"),
        reads_1=temp(TRIMMOMATIC_OUTPUT + "/{barcode}/{barcode}_R1.fastq.gz"),
        reads_2=temp(TRIMMOMATIC_OUTPUT + "/{barcode}/{barcode}_R2.fastq.gz"),
    log: TRIMMOMATIC_OUTPUT + "/{barcode}/log.txt"
    params:
        clip=TRIMMOMATIC_CLIP,
        window_size=4,
        window_qual=15,
        minlen=50,
        headcrop=15,
        leading=10,
        trailing=10,
    conda:
        CONDA_ENVS["main"]
    shell:
        """
        mkdir -p {output.folder}
        trimmomatic PE \
            -phred33 {input.reads_1} {input.reads_2} \
            {output.reads_1} {output.folder}/{wildcards.barcode}_forward_unpaired.fq.gz \
            {output.reads_2} {output.folder}/{wildcards.barcode}_reverse_unpaired.fq.gz \
            ILLUMINACLIP:{params.clip}:2:30:10 \
            LEADING:{params.leading} TRAILING:{params.trailing} \
            HEADCROP:{params.headcrop} \
            SLIDINGWINDOW:{params.window_size}:{params.window_qual} \
            MINLEN:{params.minlen} 2>&1 | tee {log}

        rm {output.folder}/*_unpaired.fq.gz 
        """


rule short_clean_reads:
    input:
        reads_1=rules.trimmomatic.output.reads_1,
        reads_2=rules.trimmomatic.output.reads_2,
    output:
        reads_1=SHORT_CLEAN_READS_OUTPUT + "/{barcode}/{barcode}_R1.fastq.gz",
        reads_2=SHORT_CLEAN_READS_OUTPUT + "/{barcode}/{barcode}_R2.fastq.gz",
        folder=directory(SHORT_CLEAN_READS_OUTPUT + "/{barcode}"),
    shell:
        """
        cp {input.reads_1} {output.reads_1}
        cp {input.reads_2} {output.reads_2}
        """

rule fastqc_after:
    """FastQC quality control"""
    input:
        rules.short_clean_reads.output.reads_1,
        rules.short_clean_reads.output.reads_2,
    output:
        folder=directory(FASTQC_AFTER_OUTPUT + "/{barcode}"),
    shell:
        """
        mkdir -p {output}
        fastqc {input} -o {output.folder}
        """


rule multiqc:
    """Create multiqc reports"""
    input:
        expand(FASTQC_AFTER_OUTPUT + "/{barcode}", **wildcards_dict),
        expand(FASTQC_BEFORE_OUTPUT + "/{barcode}", **wildcards_dict),
    output:
        folder=directory(MULTIQC_OUTPUT),
        multi_before=directory(MULTIQC_OUTPUT + "/before"),
        multi_after=directory(MULTIQC_OUTPUT + "/after"),
    params:
        fastqc_before=FASTQC_BEFORE_OUTPUT,
        fastqc_after=FASTQC_AFTER_OUTPUT,
    conda:
        CONDA_ENVS["main"]
    shell:
        """
        mkdir -p {output.folder}
        multiqc {params.fastqc_before} -o {output.multi_before}
        multiqc {params.fastqc_after} -o {output.multi_after}
        """




###############################################
# --------- Long reads processing ----------- #
###############################################
rule prepare_long_reads:
    input:
        reads= LONG_READ_INPUT + '/{barcode}/{barcode}.fastq.gz',
    output:
        folder=directory(PREPARED_LONG_READS_OUTPUT + "/{barcode}"),
        reads=PREPARED_LONG_READS_OUTPUT + "/{barcode}/{barcode}.fastq.gz",
    shell:
        """
        mkdir -p {output.folder}
        cp {input.reads} {output.reads}
        """

rule nanoplot_before:
    input:
        reads=rules.prepare_long_reads.output.reads,
    output:
        folder=directory(NANOPLOT_QC_BEFORE_OUTPUT + "/{barcode}"),
    threads: 1
    conda:
        CONDA_ENVS["main"]
    shell:
        "NanoPlot -t {threads} --fastq {input.reads} -o {output.folder}"

rule porechop:
    input:
        reads=rules.prepare_long_reads.output.reads,
    output:
        folder=directory(PORECHOP_OUTPUT + "/{barcode}"),
        reads=temp(PORECHOP_OUTPUT + "/{barcode}/{barcode}.fastq.gz"),
    log:
        PORECHOP_OUTPUT + "/{barcode}/porechop.log",
    threads: 30
    conda:
        CONDA_ENVS["main"]
    shell:
        "porechop -i {input.reads} -o {output.reads} -t {threads} 2>&1 | tee {log}"


rule filtlong:
    """
    Quality control long reads
    """
    input:
        reads=rules.porechop.output.reads,
    output:
        folder=directory(FILTLONG_OUTPUT + "/{barcode}"),
        reads=temp(FILTLONG_OUTPUT + "/{barcode}/{barcode}.fastq.gz"),
    conda:
        CONDA_ENVS["main"]
    params:
        min_length = 550, # In order to remove 501bp reads (they appear when there are a lot of adapters in the run)
        keep_percent = 95,
    shell:
        """
        mkdir -p {output.folder}
        filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} {input.reads} | gzip > {output.reads}
        """


rule long_clean_reads:
    input:
        reads=rules.filtlong.output.reads,
    output:
        reads=LONG_CLEAN_READS_OUTPUT + "/{barcode}/{barcode}.fastq.gz",
        folder=directory(LONG_CLEAN_READS_OUTPUT + "/{barcode}"),
    shell:
        """
        cp {input.reads} {output.reads}
        """


rule nanoplot_after:
    input:
        reads=rules.long_clean_reads.output.reads,
    output:
        directory(NANOPLOT_QC_AFTER_OUTPUT + "/{barcode}"),
    threads: 1
    conda:
        CONDA_ENVS["main"]
    shell:
        "NanoPlot -t {threads} --fastq {input.reads} -o {output}"




###############################################
# -------------- Assembly ------------------- #
###############################################

rule assembly:
    """
    Assembly using Unicycler
    """
    input:
        reads_1 = rules.short_clean_reads.output.reads_1,
        reads_2 = rules.short_clean_reads.output.reads_2,
        long_reads = rules.long_clean_reads.output.reads,
    output:
        folder=directory(ASSEMBLY_OUTPUT + "/{barcode}"),
        assembly=ASSEMBLY_OUTPUT + "/{barcode}/{barcode}_assembly.fasta",
        graph=ASSEMBLY_OUTPUT + "/{barcode}/{barcode}_assembly.gfa",
    params:
        sr_mode = "normal", # {conservative,normal}
        hybrid_mode = "normal", 
    threads: 12
    conda:
        CONDA_ENVS["main"]
    shell:
        """
        mkdir -p {output.folder}
        if [[ -e "{input.long_reads}" ]]
        then
            echo "Long reads are present for {wildcards.barcode}, will do hybrid assembly"
            unicycler -1 {input.reads_1} -2 {input.reads_2} -l {input.long_reads}  -o {output.folder} -t {threads} --mode {params.hybrid_mode}
        else
            echo "Long reads are not present for {wildcards.barcode}, will do short-read assembly"
            unicycler -1 {input.reads_1} -2 {input.reads_2} -o {output.folder} -t {threads} --mode {params.sr_mode}
        fi

        mv {output.folder}/assembly.fasta {output.assembly}
        mv {output.folder}/assembly.gfa {output.graph}
        """




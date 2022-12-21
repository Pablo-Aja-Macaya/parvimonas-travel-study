import sys
sys.path.append(".")

# Create DAG tree visualization
# snakemake -s snakefile_illumina.smk --rulegraph | dot -Tsvg > dag.svg

###############################################
# ------- Libraries/Functions set up -------- #
###############################################

import glob
import os
import shutil
import pandas as pd
import re
from subprocess import call, run
import pathlib

from Bio import SeqIO

# ---- Tools and databases ----

# ---- Functions ----
from utils.assembly_utils import rename_contigs
from utils.snakemake_common import replace_proxy_rule_output
from utils.process_pirate import new_or_missing_pirate
from utils.pgap_to_prokka import pgap_gff_to_prokka_gff, extract_ncl_seqs_from_gbk, clean_pgap_fasta
from utils.prepare_prokka_for_synima import prepare_prokka_for_synima

###############################################
# ------------ Config variables ------------- #
###############################################

# Output folder for pipeline
BASE_FOLDER = os.path.abspath('/home/usuario/Proyectos/Results/Kelly/Parvimonas/pgap_pipeline/pipeline')

# Inputs
GENOME_INPUT_FOLDER = os.path.abspath("/home/usuario/Proyectos/Results/Kelly/Parvimonas/pgap_pipeline/data/input_for_pipeline")

# Reference
# unkown_GCA_900637905.1_52756_G01_genomic.fna

# pgap
PGAP_PATH = os.path.abspath('~/Proyectos/programas/pgap/pgap.py')
SYNIMA_PATH = os.path.abspath("/home/usuario/Proyectos/programas/Synima")

# Environment
CONDA_ENV = os.path.abspath("envs/parvimonas.yml")

# Dont touch
INPUT_METADATA = {
    "GCF_016552165.1":{'locus_tag':'SAMN15967530','genus_species':'Parvimonas parva'},
    "GCA_000154405.1":{'locus_tag':'SAMN00627099','genus_species':'Parvimonas micra'},
    "GCA_000493795.1":{'locus_tag':'SAMN02471847','genus_species':'Parvimonas micra'},
    "GCA_000800295.1":{'locus_tag':'SAMN03135876','genus_species':'Parvimonas micra'},
    "GCA_003454775.1":{'locus_tag':'SAMN09781070','genus_species':'Parvimonas micra'},
    "GCA_003938725.1":{'locus_tag':'SAMN10163193','genus_species':'Parvimonas micra'},
    "GCA_008633155.1":{'locus_tag':'SAMD00149547','genus_species':'Parvimonas micra'},
    "GCA_900637905.1":{'locus_tag':'SAMEA48415918','genus_species':'Parvimonas micra'},
    "GCA_902373675.1":{'locus_tag':'SAMEA5850804','genus_species':'Parvimonas micra'},
    "GCA_023147965.1":{'locus_tag':'SAMN19899812','genus_species':'Parvimonas micra'},
    "GCA_023148065.1":{'locus_tag':'SAMN19899808','genus_species':'Parvimonas micra'},
    "GCA_023148105.1":{'locus_tag':'SAMN19899807','genus_species':'Parvimonas micra'},
    "GCA_023148225.1":{'locus_tag':'SAMN19899803','genus_species':'Parvimonas micra'},
    "PM102KC_G_1":{'locus_tag':'102KCG1','genus_species':'Parvimonas micra'},
    "PM114KC_G_1":{'locus_tag':'14KCG1','genus_species':'Parvimonas micra'},
    "PM79KC_AC_1":{'locus_tag':'79KCAC1','genus_species':'Parvimonas micra'},
    "PM79KC_AC_2":{'locus_tag':'79KCAC2','genus_species':'Parvimonas micra'},
    "PM79KC_AC_3":{'locus_tag':'79KCAC3','genus_species':'Parvimonas micra'},
    "PM79KC_AC_4":{'locus_tag':'79KCAC4','genus_species':'Parvimonas micra'},
    "PM79KC_G_1":{'locus_tag':'79KCG1','genus_species':'Parvimonas micra'},
    "PM89KC_AC_1":{'locus_tag':'89KCAC1','genus_species':'Parvimonas micra'},
    "PM89KC_G_1":{'locus_tag':'89KCG1','genus_species':'Parvimonas micra'},
    "PM89KC_G_2":{'locus_tag':'89KCG2','genus_species':'Parvimonas micra'},
    "PM94KC_G_1":{'locus_tag':'94KCG1','genus_species':'Parvimonas micra'},
}


###############################################
# -------------- Folder set up -------------- #
###############################################

# ---- Input preparation ----
INPUT_PREPARATION = f'{BASE_FOLDER}/01_input_preparation'
TRANSFORMED_ASSEMBLY_OUTPUT = f'{INPUT_PREPARATION}/rename_contigs'
CIRCLATOR_OUTPUT = f'{INPUT_PREPARATION}/circlator'

# ---- PGAP ----
PGAP_FOLDER = f'{BASE_FOLDER}/02_pgap'
PGAP_PREP_OUTPUT = f'{PGAP_FOLDER}/input'
PGAP_OUTPUT = f'{PGAP_FOLDER}/output'

PGAP_TO_PROKKA_OUTPUT = f'{BASE_FOLDER}/03_converted_to_prokka'


# ---- Comparative analysis ----
COMPARATIVE_PIPELINE_OUTPUT = f'{BASE_FOLDER}/comparative_analysis'

PIRATE_OUTPUT = f'{COMPARATIVE_PIPELINE_OUTPUT}/03_pirate'
PIRATE_INPUT_FOLDER = f'{PIRATE_OUTPUT}/input_folder'
PIRATE_OUTPUT_FOLDER = f'{PIRATE_OUTPUT}/output_folder'
NEW_OR_MISSING_OUTPUT = f'{PIRATE_OUTPUT}/new_or_missing_pirate'
PIRATE_RAXML_OUT = f'{PIRATE_OUTPUT}/raxml'

# ---- Synteny ----
SYNTENY_PIPELINE_OUTPUT = f'{BASE_FOLDER}/synteny'
TRANSFORMED_PROKKA_OUTPUT = f'{SYNTENY_PIPELINE_OUTPUT}/01_transformed_prokka'
SYNIMA_OUTPUT = f'{SYNTENY_PIPELINE_OUTPUT}/02_synima'


###############################################
# ----------------- Input ------------------- #
###############################################
barcodes = glob_wildcards(GENOME_INPUT_FOLDER + '/{barcode}.fasta')[0]

wildcards_dict = {
    'barcode':barcodes,
}


###############################################
# ---------------- Output ------------------- #
###############################################
rule all:
    input:
        expand(PGAP_TO_PROKKA_OUTPUT + '/{barcode}/{barcode}.fna', zip, **wildcards_dict),
        PIRATE_OUTPUT_FOLDER + '/pangenome_alignment.fasta',
        NEW_OR_MISSING_OUTPUT + '/.touch.done',
        PIRATE_RAXML_OUT + '/RAxML_bestTree.png',
        # ---- Synima ----
        expand(TRANSFORMED_PROKKA_OUTPUT + '/{barcode}/.finished.touch', zip, **wildcards_dict),
        SYNIMA_OUTPUT + '/.finished_blast.touch',
        SYNIMA_OUTPUT + '/.finished_dag.touch',


###############################################
# ------------ Prepare input ---------------- #
###############################################

rule rename_assembly:
    input:
        assembly = GENOME_INPUT_FOLDER + '/{barcode}.fasta'
    output:
        folder = directory(TRANSFORMED_ASSEMBLY_OUTPUT + '/{barcode}'),
        assembly = TRANSFORMED_ASSEMBLY_OUTPUT + '/{barcode}/{barcode}.fasta'
    run:
        rename_contigs(input.assembly, output.assembly, wildcards)


rule conditional_circlator:
    input:
        assembly = rules.rename_assembly.output.assembly
    output:
        folder = directory(CIRCLATOR_OUTPUT + '/{barcode}'),
        assembly = CIRCLATOR_OUTPUT + '/{barcode}/{barcode}.fasta',
        check = touch(CIRCLATOR_OUTPUT + '/{barcode}/.done.touch'),
    conda:
        CONDA_ENV
    threads: 6
    shell:
        """
        # WARNING: Only valid for genomes that are expected to have only one chromosome, with no plasmids (it counts number of contigs)
        seq_count=$(grep -c "^>" {input.assembly})
        if [ $seq_count -eq 1 ]
            then
                echo "Only one sequence in {input.assembly}"
                echo "Running circlator"
                circlator fixstart {input.assembly} {output.assembly}
                mv {output.assembly}.fasta {output.assembly} # circlator puts .fasta at the end so it is renamed
                rm {output.assembly}.*promer*
                
        elif [ $seq_count -gt 1 ]
            then
                echo "More than one sequence in {input.assembly}"
                echo "Setting linear mode (not modifying file)"
                cat {input.assembly} > {output.assembly}
        fi        
        """

###############################################
# ----------------- PGAP -------------------- #
###############################################

def get_locus_tag(wildcards, input_metadata=INPUT_METADATA):
    return input_metadata[wildcards.barcode]['locus_tag']

def get_genus_species(wildcards, input_metadata=INPUT_METADATA):
    return input_metadata[wildcards.barcode]['genus_species']

rule pgap_preparation:
    input:
        assembly = rules.conditional_circlator.output.assembly,
        check = rules.conditional_circlator.output.check,
    output:
        folder = directory(PGAP_PREP_OUTPUT + '/{barcode}'),
        assembly = PGAP_PREP_OUTPUT + '/{barcode}/{barcode}.fna',
        input_yaml = PGAP_PREP_OUTPUT + '/{barcode}/input.yaml',
        submol_yaml = PGAP_PREP_OUTPUT + '/{barcode}/submol.yaml',
    params:
        locus_tag = get_locus_tag,
        genus = get_genus_species,
    shell:
        """
        seq_count=$(grep -c "^>" {input.assembly})
        if [ $seq_count -eq 1 ]
            then
                echo "Only one sequence in {input.assembly}"
                echo "Setting circular mode"
                completness_setting="circular"
                
        elif [ $seq_count -gt 1 ]
            then
                echo "More than one sequence in {input.assembly}"
                echo "Setting linear mode"
                completness_setting="linear"
        fi    
        
        cp {input.assembly} {output.assembly}

        # Create input.yaml
        > {output.input_yaml}
        echo "fasta:" >> {output.input_yaml}
        echo "    class: File" >> {output.input_yaml}
        echo "    location: {wildcards.barcode}.fna" >> {output.input_yaml}
        echo "submol:" >> {output.input_yaml}
        echo "    class: File" >> {output.input_yaml}
        echo "    location: submol.yaml" >> {output.input_yaml}

        # Create submol.yaml
        > {output.submol_yaml}
        echo "topology: ${{completness_setting}}" >> {output.submol_yaml}
        echo "location: 'chromosome'" >> {output.submol_yaml}
        echo "organism:" >> {output.submol_yaml}
        echo "    genus_species: '{params.genus}'" >> {output.submol_yaml}
        echo "    strain: '{wildcards.barcode}'" >> {output.submol_yaml}     
        echo "locus_tag_prefix: '{params.locus_tag}'" >> {output.submol_yaml}  

        echo "Required PGAP files ready. Manually check their atributes"

        """


rule pgap:
    input:
        input_yaml = rules.pgap_preparation.output.input_yaml
    output:
        folder = directory(PGAP_OUTPUT + '/{barcode}'),
        assembly = PGAP_OUTPUT + '/{barcode}/{barcode}.fna',
        assembly_gbk = PGAP_OUTPUT + '/{barcode}/{barcode}.gbk',
        assembly_gff = PGAP_OUTPUT + '/{barcode}/{barcode}.gff',
    params:
        pgap = PGAP_PATH
    threads: 20
    resources:
        mem_mb=25000
    shell:
        """
        # Run pgap
        rm {output.folder} -r # remove folder because snakemake creates it first and pgap does not overwrite it
        {params.pgap} --report-usage-false --no-internet --cpus {threads} --memory {resources.mem_mb}m -o {output.folder} {input.input_yaml}
        
        # Rename with barcodes
        mv {output.folder}/{wildcards.barcode}.fna {output.folder}/{wildcards.barcode}_org.fasta
        mv {output.folder}/annot.fna {output.assembly}
        mv {output.folder}/annot.gbk {output.assembly_gbk}
        mv {output.folder}/annot.gff {output.assembly_gff}
        mv {output.folder}/annot.faa {output.folder}/{wildcards.barcode}.faa
        mv {output.folder}/annot.sqn {output.folder}/{wildcards.barcode}.sqn
        """


rule pgap_to_prokka:
    input:
        folder = rules.pgap.output.folder,
        assembly = rules.pgap.output.assembly,
        assembly_gff = rules.pgap.output.assembly_gff,
        assembly_gbk = rules.pgap.output.assembly_gbk,
    output:
        folder = directory(PGAP_TO_PROKKA_OUTPUT + '/{barcode}'),
        assembly = PGAP_TO_PROKKA_OUTPUT + '/{barcode}/{barcode}.fna',
        assembly_gff = PGAP_TO_PROKKA_OUTPUT + '/{barcode}/{barcode}.gff',
        assembly_gbk = PGAP_TO_PROKKA_OUTPUT + '/{barcode}/{barcode}.gbk',
        assembly_ffn = PGAP_TO_PROKKA_OUTPUT + '/{barcode}/{barcode}.ffn',
    run:
        # ---- GFF transformation ----
        _ = pgap_gff_to_prokka_gff(wildcards.barcode, input.assembly, input.assembly_gff, output.assembly_gff)
        
        # ---- File with all genes in nucleotide form ----
        _ = extract_ncl_seqs_from_gbk(input.assembly_gbk, output.assembly_ffn)

        # ---- Genome and proteins ----
        _ = clean_pgap_fasta(f'{input.folder}/{wildcards.barcode}.fna', f'{output.folder}/{wildcards.barcode}.fna', version='fna')
        _ = clean_pgap_fasta(f'{input.folder}/{wildcards.barcode}.faa', f'{output.folder}/{wildcards.barcode}.faa')

        # ---- Unchanged ----
        shell("cp {input.assembly_gbk} {output.assembly_gbk}")
        shell("cp {input.folder}/*.sqn {output.folder}")


###############################################
# ---------------- Synima ------------------- #
###############################################

rule modify_prokka:
    input:
        folder = rules.pgap_to_prokka.output.folder
    output:
        folder = directory(TRANSFORMED_PROKKA_OUTPUT + '/{barcode}'),
        check = touch(TRANSFORMED_PROKKA_OUTPUT + '/{barcode}/.finished.touch')
    run:
        mod_files = prepare_prokka_for_synima(wildcards.barcode, input.folder, output.folder)


rule synima_blast_all:
    input:
        prokka_folders = lambda wc: replace_proxy_rule_output(wc, rules.modify_prokka.output.folder, wildcards_dict),
    output:
        folder = directory(SYNIMA_OUTPUT),
        check = touch(SYNIMA_OUTPUT + '/.finished_blast.touch')
    conda:
        CONDA_ENV
    params:
        synima = SYNIMA_PATH
    shell:
        """
        # Copy prokka folders to the synima folder
        cp -r {input.prokka_folders} {output.folder}

        # Move there
        cd {output.folder}

        # Create the config file
        echo "//" > Repo_spec.txt
        for folder in {input.prokka_folders}
        do
            barcode=$(basename $folder)
            echo "Genome" $barcode
            echo "Annotation" $barcode
            echo "//"
        done >> Repo_spec.txt

        perl {params.synima}/util/Create_full_repo_sequence_databases.pl -r ./Repo_spec.txt -f gene
        perl {params.synima}/util/Blast_grid_all_vs_all.pl -r ./Repo_spec.txt
        """

rule synima_dag:
    input:
        check = rules.synima_blast_all.output.check
    output:
        check = touch(SYNIMA_OUTPUT + '/.finished_dag.touch')
    conda:
        CONDA_ENV
    params:
        folder = SYNIMA_OUTPUT,
        synima = SYNIMA_PATH,
    shell:
        """
        cd {params.folder}
        perl {params.synima}/util/Blast_all_vs_all_repo_to_OrthoMCL.pl -r ./Repo_spec.txt
        perl {params.synima}/util/Orthologs_to_summary.pl -o all_orthomcl.out
        perl {params.synima}/util/DAGchainer_from_gene_clusters.pl -r ./Repo_spec.txt -c GENE_CLUSTERS_SUMMARIES.OMCL/GENE_CLUSTERS_SUMMARIES.clusters
        """

rule synima_plot:
    input:
        check = rules.synima_dag.output.check
    output:
        check = SYNIMA_OUTPUT + '/.finished_plot.touch'
    conda:
        CONDA_ENV
    params:
        folder = SYNIMA_OUTPUT,
        synima = SYNIMA_PATH,
    shell:
        """
        cd {params.folder}
        perl {params.synima}/SynIma.pl -a Repo_spec.txt.dagchainer.aligncoords -b Repo_spec.txt.dagchainer.aligncoords.spans -g c
        """

# -- Manual synima --
# perl ${synima}/util/Blast_all_vs_all_repo_to_OrthoMCL.pl -r ./Repo_spec.txt
# perl ${synima}/util/Orthologs_to_summary.pl -o all_orthomcl.out
# perl ${synima}/util/DAGchainer_from_gene_clusters.pl -r ./Repo_spec.txt -c GENE_CLUSTERS_SUMMARIES.OMCL/GENE_CLUSTERS_SUMMARIES.clusters
# perl ${synima}/SynIma.pl -a Repo_spec.txt.dagchainer.aligncoords -b Repo_spec.txt.dagchainer.aligncoords.spans -g c
# cat SynIma-output__tumor/config.txt.Rscript | grep -v "^text(" > SynIma-output__tumor/config_no_text.txt.Rscript
# Rscript SynIma-output/config.txt.Rscript 

# -- Manual Sibelia --
# Sibelia -s fine -m 5000 --outdir xxx__vs__xxx input/xxxx.fna input/xxx.fna


# -- Manual mummer2circos --
# mummer2circos -r input/queries/PM89KC_G_1.fna -q input/queries/PM89KC_G_2.fna

###############################################
# ---------------- PIRATE ------------------- #
###############################################

rule copy_prokka_gff_to_pirate_folder:
    input:
        lambda wc: replace_proxy_rule_output(wc, rules.pgap_to_prokka.output.assembly_gff, wildcards_dict),
    output:
        folder = directory(PIRATE_INPUT_FOLDER)
    shell:
        """
        mkdir -p {output.folder}
        cp {input} {output.folder}
        """

rule pirate:
    input:
        rules.copy_prokka_gff_to_pirate_folder.output.folder
    output:
        folder = directory(PIRATE_OUTPUT_FOLDER),
        pan_aln = PIRATE_OUTPUT_FOLDER + '/pangenome_alignment.fasta',
        gene_families_renamed = PIRATE_OUTPUT_FOLDER + '/PIRATE.gene_families.ordered.renamed.tsv'
    threads: 60
    conda:
        CONDA_ENV
    shell:
        """
        # Run PIRATE
        # PIRATE -i {input} -o {output.folder} -t {threads} -s "40,50,60,70,80,90,95,96,97,98,99,100" -k "--cd-low 100 --e 1E-12 --diamond" -a -r
        PIRATE -i {input} -o {output.folder} -t {threads} -s "70,80,90,95,96,97,98,99,100" -k "--cd-low 100 --e 1E-12 --diamond" -a -r

        # Transform PIRATE output to contain previous locus tag
        scripts/pirate__subsample_outputs.pl -g {output.folder}/modified_gffs -i {output.folder}/PIRATE.gene_families.ordered.tsv -o {output.folder}/PIRATE.gene_families.ordered.renamed.tsv --field "prev_locus"

        """


rule new_or_missing_pirate:
    input:
        folder = rules.pirate.output.folder,
        gene_families_renamed = rules.pirate.output.gene_families_renamed
    output:
        folder = directory(NEW_OR_MISSING_OUTPUT),
        check = touch(NEW_OR_MISSING_OUTPUT + '/.touch.done')
    run:
        query_groups = { # Names must not contain dots because PIRATE replaces them by underscore
            'group_1':['PM89KC_G_1','PM89KC_G_2'],
            'group_2':['PM79KC_AC_1','PM79KC_AC_2','PM79KC_AC_3','PM79KC_AC_4'],
            'group_3':['PM89KC_AC_1'],
            'group_4':['PM94KC_G_1'],
            'group_5':['PM102KC_G_1'],
            'group_6':['PM114KC_G_1'],
            'group_7':['GCA_000154405_1'],
            'group_8':['GCA_000493795_1'],
            'group_9':['GCA_902373675_1'],
            'group_10':['GCA_003938725_1'],
            'group_11':['GCA_008633155_1'],
            'group_12':['GCA_000800295_1'],
            'group_13':['GCA_003454775_1'],
            'group_14':['GCF_016552165_1'],
            'group_15':['GCA_900637905_1'],
            'group_16': ['GCA_023147965_1','GCA_023148065_1','GCA_023148105_1','GCA_023148225_1']
        }
        new_or_missing_pirate(input.gene_families_renamed, output.folder, query_groups)



###############################################
# ---------------- RAxML -------------------- #
###############################################

rule fasta2phylip:
    input:
        'alignment.fasta'
    output:
        'alignment.phy'
    shell:
        'bash scripts/convert_fasta2phylip.sh {input} > {output}'

rule raxml:
    input:
        'alignment.phy'
    output:
        folder = directory('path' + '/raxml'),
        res = 'path' + '/raxml/RAxML_bestTree.tree',
        check = touch('path' + '/raxml/raxml.done'),
    threads: 60
    params:
        bootstrap_perms = 1000, #--bootstop-perms={params.bootstrap_perms}
        model = 'GTRCAT',
        wkdir = os.path.dirname(os.path.realpath(__file__)),
    shell:
        'mkdir -p {output.folder}'
        ' && '
        'raxmlHPC-PTHREADS-AVX -s {input} -w {params.wkdir}/{output.folder} -m {params.model} -T {threads} -n tree -p 1 -N 1000 -p 12345 -x 12345 -f a'

rule draw_raxml_tree:
    # Draw tree with R
    input:
        tree="tree.newick",
        check="raxml.done",
    output:
        tree="tree.png",
    shell:
        "Rscript scripts/draw_tree.R {input.tree} {output.tree}"


# ---- Pirate x RAxML ----
use rule fasta2phylip as pirate_fasta2phylip with:
    input:
        rules.pirate.output.pan_aln
    output:
        PIRATE_RAXML_OUT + '/alignment.phy'

use rule raxml as pirate_raxml with:
    input:
        rules.pirate_fasta2phylip.output
    output:
        folder = directory(PIRATE_RAXML_OUT + '/raxml'),
        res = PIRATE_RAXML_OUT + '/raxml/RAxML_bestTree.tree',
        check = touch(PIRATE_RAXML_OUT + '/raxml/raxml.done'),

use rule draw_raxml_tree as draw_pirate_raxml with:
    input:
        tree = rules.pirate_raxml.output.res,
        check = rules.pirate_raxml.output.check,
    output:
        tree = PIRATE_RAXML_OUT + '/RAxML_bestTree.png'



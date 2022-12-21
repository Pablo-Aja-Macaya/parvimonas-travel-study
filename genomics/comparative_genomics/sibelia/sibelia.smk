
import glob
import os
import itertools
import pandas as pd

###############################################
# ------------ Config variables ------------- #
###############################################


# Output folder
BASE_DIR = os.path.abspath("/home/usuario/Proyectos/Results/Kelly/Parvimonas/pgap_pipeline/manual/sibelia")

# Input folder with genomes ending in .fna
INPUT_FOLDER = os.path.abspath(f"{BASE_DIR}/input")

# CIRCOS config file
config_file = f"config/default_circos.conf"

# Parameters for Sibelia
window_sizes = [50,100,500,1000,2000,3000,5000]

# Environment
CONDA_ENV = "envs/sibelia.yaml"


###############################################
# -------------- Folder set up -------------- #
###############################################
SIBELIA_OUTPUT = os.path.abspath(f"{BASE_DIR}/output/sibelia")
CLEAN_CIRCOS_OUTPUT = os.path.abspath(f"{BASE_DIR}/output/clean_circos")
CIRCOS_OUTPUT = os.path.abspath(f"{BASE_DIR}/output/circos")


###############################################
# ----------------- Input ------------------- #
###############################################
barcodes = [os.path.basename(i).replace('.fna','') for i in glob.glob(f"{INPUT_FOLDER}/*.fna")]

def get_output(wc, barcodes, window_sizes, folder):
    l = []
    for w in window_sizes:
        used_combinations = []  
        for bc1, bc2 in itertools.combinations(barcodes, 2):
            if set([bc1,bc2]) not in used_combinations and bc1 != bc2:
                l.append(folder + f'/{w}/{bc1}__vs__{bc2}/{bc1}__vs__{bc2}__circos.png')
            else:
                continue
    return l

###############################################
# ---------------- Output ------------------- #
###############################################
rule all:
    input:
        lambda wc: get_output(wc, barcodes, window_sizes, CIRCOS_OUTPUT)

rule sibelia:
    input:
        assembly_1 = INPUT_FOLDER + '/{bc1}.fna',
        assembly_2 = INPUT_FOLDER + '/{bc2}.fna',
    output:
        folder = directory(SIBELIA_OUTPUT + '/{window_size}/{bc1}__vs__{bc2}'),
        circos_folder = directory(SIBELIA_OUTPUT + '/{window_size}/{bc1}__vs__{bc2}/circos'),
    threads: 5
    conda: CONDA_ENV
    shell:
        """        
        mkdir -p {output.folder}
        Sibelia -s loose -m {wildcards.window_size} --outdir {output.folder} {input.assembly_1} {input.assembly_2} 
        """

rule clean_reciprocal_hits:
    """
    Removes reciprocal hits in circos.segdup.txt (links will only unite to the other chromosome)
    """
    input:
        circos_folder = rules.sibelia.output.circos_folder
    output:
        folder = directory(CLEAN_CIRCOS_OUTPUT + '/{window_size}/{bc1}__vs__{bc2}'),
        circos_folder = directory(CLEAN_CIRCOS_OUTPUT + '/{window_size}/{bc1}__vs__{bc2}/circos'),
    run:
        try:
            shell("cp -r {input.circos_folder} {output.folder}")
            df = pd.read_csv(f'{output.circos_folder}/circos.segdup.txt', delim_whitespace=True, header=None)
            df = df.drop_duplicates([0,1], keep=False)
            df.to_csv(f'{output.circos_folder}/circos.segdup.txt', index=None, sep=' ', header=None)
        except Exception as e:
            raise Exception(e)


rule circos:
    input:
        circos_folder = rules.clean_reciprocal_hits.output.circos_folder
    output:
        folder = directory(CIRCOS_OUTPUT + '/{window_size}/{bc1}__vs__{bc2}'),
        circos_folder = directory(CIRCOS_OUTPUT + '/{window_size}/{bc1}__vs__{bc2}/circos'),
        image = CIRCOS_OUTPUT + '/{window_size}/{bc1}__vs__{bc2}/{bc1}__vs__{bc2}__circos.png',
    conda: CONDA_ENV
    params:
        config_file = config_file
    shell:
        """
        cp -r {input.circos_folder} {output.folder}
        cat {params.config_file} > {output.circos_folder}/circos.conf
        cd {output.circos_folder} && circos && cd {output.folder}
        mv {output.circos_folder}/circos.png {output.image}
        """


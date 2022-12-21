from colorama import Fore, Style
import pandas as pd
from subprocess import run
import re
import glob
import os

import snakemake
from snakemake.io import expand


def replace_proxy_rule_output(
    wildcards, expandable_structure, wildcards_dict, zip_expansion="yes"
):
    """
    Takes path structure (fakepath/{barcode}/{barcode}_assembly.fasta) and expands the ProxyRule
    Inputs:
        wildcards = snakemake wildcard object produced in rule (it is not used in the function)
        expandable_structure =  "fakepath/{wildcard1}/{wildcard2}_assembly.fasta"
        wildcards_dict = {'wildcard1':['A','B','C], 'wildcard2':['D','E','F']}
    """
    if zip_expansion == "yes":
        l = expand(expandable_structure, zip, **wildcards_dict)
    else:
        l = expand(expandable_structure, **wildcards_dict)

    if l:
        return l
    else:
        raise Exception(
            f"Error: Structure {expandable_structure} does not contain any known wildcards."
        )


# def run_snakemake(snakefile: str, parameters: dict, ppln_name: str):
#     snakemake.snakemake(
#         snakefile,
#         printshellcmds=False,
#         forceall=True,
#         force_incomplete=True,
#         # workdir=config[KEY_TEMPDIR],
#         config=config,
#         cores=parameters["cores"],
#         lock=False,
#         quiet=True,
#         # log_handler=logger.log_handler,
#     )


def run_snakemake(snakefile: str, parameters: dict, ppln_name: str):
    print(f"{Fore.YELLOW}Ititializing {ppln_name} pipeline... {Style.RESET_ALL}\n")

    # Init command
    cmd = [f"snakemake -s {snakefile}"]

    # Add previz if required
    if parameters["previsualize"] == "yes":
        cmd.append("-n --quiet")

    # Add cores
    cmd.append(f"-c{parameters['cores']}")

    # Use conda (can be optional)
    if parameters["use_conda"]:
        cmd.append("--use-conda")

    # Params to df
    config_df = pd.DataFrame.from_dict(parameters.items())
    config_df = config_df.set_index(0).T
    print(f"{Fore.YELLOW}User supplied params:{Style.RESET_ALL}")
    print(f"{Fore.LIGHTBLACK_EX}{config_df}{Style.RESET_ALL}\n")

    # Create folder and params file
    user_config_file = f"{parameters['output_folder']}/user_config_file.tsv"
    _ = run(f"mkdir -p {parameters['output_folder']}", check=True, shell=True)
    config_df.to_csv(user_config_file, sep="\t", index=None)

    # TODO: make this optional
    cmd.append(
        "--rerun-triggers input mtime"
    )  # default: ['mtime', 'params', 'input', 'software-env', 'code']

    # Add config file path to command
    cmd.append("--config")
    cmd.append(f"user_config_file={user_config_file}")

    # Create workflow image (HAS TO BE THE LAST ARGUMENT)
    draw_wf = parameters.get("draw_wf")
    if draw_wf:
        draw_wf_file = draw_wf
        print(
            f"{Fore.YELLOW}Creating workflow image (rulegraph) in {draw_wf_file}...{Style.RESET_ALL}"
        )
        print(
            f"{Fore.LIGHTBLACK_EX}Note: If image is created then steps do not show in table{Style.RESET_ALL}"
        )
        cmd.append(f"--rulegraph | dot -Tpdf > {draw_wf_file}")

    # Execute
    cmd = " ".join(cmd)
    try:
        cmd_res = run(cmd, check=True, shell=True)
    except Exception as e:
        print(f"{Fore.RED}ERROR: {e}{Style.RESET_ALL}")

    print(f"\n{Fore.MAGENTA}Used command: {cmd}{Style.RESET_ALL}\n")


def discern_reads(reads: str, read_type: str) -> dict:
    """
    Given a file path for reads and read_type, discern the main ID of the sample and which type it is,
    it then replaces the ID of the sample with the string "{barcode}", preparing it for snakemake
    Input:
    - reads: file path to reads in fastq.gz
    - read_type: one of the following ["illumina_pe_reads","illumina_se_reads","ont_reads"]
    """
    structure_dict = {
        "illumina_pe_reads": {
            "pe_basespace": {  # Paired-end reads from basespace
                "re": r"/(.*)_.*ds\..*/\1_S\d+_L\d+_R[12]_\d+.fastq.gz",
                "smk_format": "{folder_id}/{barcode}_S{trash1}_R{orientation}_{trash2}.fastq.gz",
            },
            "pe_any": {  # Paired-end reads from other place (takes whole id)
                "re": r"/(.*)/\1_R[12].fastq.gz",
                "smk_format": "{folder_id}/{barcode}_R{orientation}.fastq.gz",
            },
        },
        "illumina_se_reads": {
            "se_reads": {  # Single-end reads (?)
                "re": r"/(.*)/\1.fastq.gz",
                "smk_format": "{folder_id}/{barcode}.fastq.gz",
            },
        },
        "ont_reads": {
            "ont_guppy": {  # Reads coming from Guppy
                "re": r"/(barcode\d+[A-Z]*)/\1.fastq.gz",
                "smk_format": "{folder_id}/{barcode}.fastq.gz",
            },
            "ont_porechop": {  # Reads coming from Porechop
                "re": r"/(bc\d+)[A-Z]*/\1.fastq.gz",
                "smk_format": "{folder_id}/{barcode}.fastq.gz",
            },
            "ont_any": {  # Reads from other place (takes the whole id)
                "re": r"/(.*)/\1.fastq.gz",
                "smk_format": "{folder_id}/{barcode}.fastq.gz",
            },
        },
    }
    res_dict = {}
    for structure_name, structure in structure_dict[read_type].items():
        structure_re = structure["re"]
        res = re.search(structure_re, reads)
        if res:
            # print(f"FOUND | ID: {res.group(1)} | STR: {structure_name}")
            barcode = res.group(1)
            smk_format = structure["smk_format"]
            res_dict = {
                "barcode": barcode,
                "found": True,
                "structure_name": structure_name,
                "smk_format": smk_format,
            }
            return res_dict
        else:
            pass
            # print(f"NOT FOUND | STR: {structure_name} | READS: {reads}")

    return {"barcode": None, "found": False}
    # raise Exception(f"Could not identify structure of file {reads}")


def check_reads_substitution_structure(parent_folder: str, read_type: str) -> str:
    """
    Find substitution structure for reads
    The input is the parent_folder, which contains multiple folders, and
    each folder can have 1 or 2 files with reads ending in ".fastq.gz".
    These files can be from illumina or ONT, so this is discerned in discern_reads()
    For each file, a structure is returned (if it follows any)
    If all structures are the same, the function passes and returns the final common structure,
    if not it raises an error
    The final common structure can then be used in snakemake.io.glob_wildcards()

    Input:
    - parent_folder: path to folder containing folders with reads
    - read_type: "ont_reads" | "illumina_pe_reads" | "illumina_se_reads"
    Output:
    - common_structure: "{folder_id}/{barcode}.fastq.gz"
    - common_structure_name: ['ont_guppy']
    Possible exceptions:
    - No reads in subfolder
    - Could not find matching structure
    - Not all files follow same substitution structure
    - No available files
    """
    reads_folders = glob.glob(f"{parent_folder}/*/")
    smk_formats_found = []
    structure_names_found = []
    for folder in reads_folders:
        reads_files = glob.glob(f"{folder}/*.fastq.gz")
        if len(reads_files) == 0:
            raise Exception(f"ERROR: No files in folder: {folder}")
        else:
            for reads in reads_files:
                res = discern_reads(reads, read_type)
                if res["found"]:
                    # print(f"{reads} | {res['barcode']} | {res['smk_format']}")
                    # print(res)
                    smk_formats_found.append(res["smk_format"])
                    structure_names_found.append(res["structure_name"])
                else:
                    raise Exception(
                        f"ERROR: Could not find structure using {read_type} preset: {reads}"
                    )

    smk_formats_found = set(smk_formats_found)
    structure_names_found = set(structure_names_found)
    if len(smk_formats_found) > 1:
        raise Exception(f"ERROR: Not all files follow the same substitution structure")
    elif len(smk_formats_found) == 0:
        raise Exception(f"ERROR: No input files available in {parent_folder}")
    else:
        common_structure = list(smk_formats_found)[0]
        common_structure_name = list(structure_names_found)
        # print(
        #     f"All files in {parent_folder} follow same structure: Structure='{common_structure}' | Structure names={common_structure_name}"
        # )
        return common_structure, common_structure_name

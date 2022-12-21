import glob
import pandas as pd
from pathlib import Path
from subprocess import run

############################################################
# ---------------------- Params -------------------------- #
############################################################
# ---- Input ----
# VFDB fasta and diamond database of this fasta
vfdb_fasta = "/home/usuario/Proyectos/Innova/data/db/virulence/VFDB/diamond_db/VFDB_A_pro/VFDB_setA_pro.fas"
database = "/home/usuario/Proyectos/Innova/data/db/virulence/VFDB/diamond_db/VFDB_A_pro/database.dmnd"

# Folder with subfolders of structure "in_folder/{sample_id}/{sample_id}.faa"
in_folder = "/home/usuario/Proyectos/Results/Kelly/Parvimonas/pgap_pipeline/pipeline/03_converted_to_prokka"

# ---- Output ----
# Output folder
out_folder = "/home/usuario/Proyectos/Results/tests/virulence"

# ---- Params ----
outfmt_fields = "qseqid qtitle sseqid stitle pident length mismatch gapopen full_qseq full_sseq qstart qend qlen sstart send slen evalue bitscore"
outfmt = 6
max_target_seqs = 1
threads = 60
max_hsps = 1
min_identity = 25

############################################################
# ------------------- Run DIAMOND ------------------------ #
############################################################
cmd = f"mkdir -p {out_folder}"
print(cmd)
run(cmd, shell=True, check=True)

for input_folder in glob.glob(f"{in_folder}/*"):
    bc = Path(input_folder).name
    query = f"{input_folder}/{bc}.faa"
    output_file = f"{out_folder}/{bc}/vfdb_dmnd.tsv"

    cmd = f"mkdir -p {out_folder}/{bc}"
    print(cmd)
    run(cmd, shell=True, check=True)

    cmd = [
        "diamond blastp",
        f"--in {vfdb_fasta}",
        f"--db {database}",
        f"--query {query}",
        f"--outfmt {outfmt} {outfmt_fields}",
        f"--header -o {output_file} --max-target-seqs {max_target_seqs}",
        f"--max-hsps {max_hsps} --threads {threads} --salltitles",
        f"--quiet --more-sensitive",
    ]
    cmd = " ".join(cmd)
    print(cmd)
    run(cmd, shell=True, check=True)


############################################################
# ----------------- Process DIAMOND ---------------------- #
############################################################
def process_diamond_vfdb(
    diamond_res, outfmt_fields, output_file=None, min_identity=70, min_length_pct=60
):
    """
    Process results from Diamond on VFDB database
    """
    diamond_res = diamond_res
    diamond_fields = outfmt_fields.split()

    with open(diamond_res, "rt") as handle:
        for line in handle:
            if "Fields:" in line:
                line = line.strip()
                line = line.replace("Fields: ", "").replace("# ", "")
                diamond_cols_long_names = line.split(", ")
                break

    diamond_field_correspondence = {
        short: long for short, long in zip(diamond_fields, diamond_cols_long_names)
    }

    df = pd.read_csv(diamond_res, sep="\t", skiprows=3, header=None)
    df.columns = diamond_fields

    # ---- Add VFDB annotation to result ----
    records = df["stitle"].tolist()

    dicc_list = []
    for record in records:
        try:
            # Get ID and remove from record
            id = record[: record.find(" ")]
            record = record.replace(id + " ", "")

            # Get gene name and remove from record
            gene_name = record[: record.find(" ")]
            record = record.replace(gene_name + " ", "")
            gene_name = gene_name.replace("(", "").replace(")", "")

            # Get full description
            description = record[: record.find("[") - 1]
            record = record.replace(description + " ", "")

            # Get second ID and organism
            product, organism = [
                i.replace("[", "").replace("]", "") for i in record.split("] [")
            ]

            dicc = {
                "sseqid": id,
                "gene_name": gene_name,
                "description": description,
                "product": product,
                "organism": organism,
            }
            dicc_list.append(dicc)

        except Exception:
            print(f"Error processing: {record}")

    subject_desc_df = pd.DataFrame(dicc_list)

    # Merging of diamond result with VFDB information
    merged_df = df.merge(subject_desc_df, on="sseqid")
    merged_df = merged_df.drop_duplicates()

    # Filter
    merged_df = merged_df[merged_df["pident"] >= min_identity]
    merged_df = merged_df[
        merged_df["length"] / merged_df["slen"] * 100 >= min_length_pct
    ]

    # Order
    merged_df = merged_df.sort_values(by=["pident", "bitscore"], ascending=False)

    # Create length of alignment / length of subject gene field
    merged_df["length/slen"] = (
        merged_df["length"].astype(str) + "/" + merged_df["slen"].astype(str)
    )

    # ---- Select output columns ----
    merged_df = merged_df[
        [
            "qseqid",
            "sseqid",
            "gene_name",
            "product",
            "description",
            "organism",
            "pident",
            "length/slen",
            "mismatch",
            "gapopen",
            "evalue",
            "bitscore",
        ]
    ]  # ,'full_qseq','full_sseq'

    # ---- Rename columns ----
    diamond_field_correspondence["gene_name"] = "Gene"
    diamond_field_correspondence["product"] = "Product"
    diamond_field_correspondence["description"] = "Description"
    diamond_field_correspondence["organism"] = "Organism"
    diamond_field_correspondence["length/slen"] = "Alignment length / Gene length"
    diamond_field_correspondence["mismatch"] = "Mismatch"
    diamond_field_correspondence["qseqid"] = "Query"
    diamond_field_correspondence["sseqid"] = "Subject"
    # diamond_field_correspondence['qtitle'] = 'Query description'

    final_df = merged_df.rename(columns=diamond_field_correspondence)

    if output_file:
        final_df.to_csv(output_file, sep="\t", index=None)

    return final_df


for folder in glob.glob(f"{out_folder}/*"):
    try:
        diamond_res = f"{folder}/vfdb_dmnd.tsv"
        _ = process_diamond_vfdb(
            diamond_res,
            outfmt_fields,
            output_file=f"{folder}/vfdb_dmnd_processed.tsv",
            min_identity=min_identity,
        )
    except:
        pass
